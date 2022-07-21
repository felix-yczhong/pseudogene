import subprocess
import tempfile
import pathlib
import re
import functools
import itertools

import multiprocessing
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

import smaca.constants as C
import pseudogene
from pseudogene.bam import PseudoBam
from pseudogene.utils import *
from pseudogene.worker import *

def read_input_RawArray(file_path, ref, true_gene_region, pseudo_gene_region, true_gene_pos, pseudo_gene_pos, true_genes, pseudo_genes, idx):
    bam_file = PseudoBam(file_path)

    global control_gene_cov, gene_region_cov, true_table, pseudo_table
    control_gene_cov[idx, ...] = bam_file.get_cov_ranges(C.POSITIONS[ref]["GENES"])
    gene_region_cov[idx, ...] = bam_file.get_cov_ranges_(true_gene_region, pseudo_gene_region)

    offset = 0
    for true_gene, pos in true_gene_pos.items():
        for j, pos_info in enumerate(pos):
            true_table[idx, offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(pos_info[:2])).squeeze()
        offset += len(pos)

    offset = 0
    for pseudo_gene, pos in pseudo_gene_pos.items():
        for j, pos_info in enumerate(pos):
            pseudo_table[idx, offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(pos_info[:2])).squeeze()
        offset += len(pos)

# calculator for a certain gene
class PseudoGeneCalculator():
    def read_config(self, config, gene_group, ref):
        self.ref = ref
        self.gene_group = gene_group
        self.config = config
        self.true_gene = config[gene_group]['true_gene']
        self.pseudo_gene = config[gene_group]['pseudo_gene']
        self.num_gene = len(self.true_gene) + len(self.pseudo_gene)

        self.true_exception = {true_gene: config[self.gene_group][self.ref][true_gene]['exception'] for true_gene in self.true_gene}
        self.pseudo_exception = {pseudo_gene: config[self.gene_group][self.ref][pseudo_gene]['exception'] for pseudo_gene in self.pseudo_gene}

    def get_gene_region(self):
        """
        read fasta files created by query_alignment() to get true gene(s) and pseudo gene(s) region information, specifically
        * for true gene(s): exon number, chromosome, start, end
        * for pseudo gene(s): chromosome, start, end

        :param self:
        :return true_gene_info: a list of (exon number, chromosome, start, end) tuple in the same order of self.true_gene
        :return pseudo_gene_info: a list of (chromosome, start, end) tuple in the same order of self.pseudo_gene
        """
        def read_fasta(true_gene, ref):
            input_path = self.config['data']['fasta']['true_seq'].format(gene=true_gene, genome_ver=ref)
            with open(input_path) as fasta:
                header = fasta.readline()
            ma = prog.match(header)
            if ma is None:
                raise ValueError(f"Cannot recognize {input_path} header")
            chr, start, end = ma.group('chr', 'start', 'end')
            start, end = convert_coordiante_fasta_to_sam(int(start), int(end))
            return chr, start, end
        
        pattern = ">(?P<gene_name>[\w]+[\d]*)-(?P<exon_num>[\d]*)-(?P<chr>chr[XY\d]+)-(?P<start>[\d]+)-(?P<end>[\d]+)"
        prog = re.compile(pattern)

        true_gene_region = [read_fasta(true_gene, self.ref) for true_gene in self.aligned_true_gene]
        pseudo_gene_region = [read_fasta(pseudo_gene, self.ref) for pseudo_gene in self.aligned_pseudo_gene]
        return true_gene_region, pseudo_gene_region

    def construct_alignment_table(self):
        """
        read bam files created by query_alignment() to get true gene(s) and pseudo gene(s) alignment information and contruct alignment table

        :param self:
        :return alignment_table: a pandas dataframe where each row represents how true gene(s) alignment position(s) align to pseudo gene(s) alignment position(s)
        :return true_gene_pos: a dictionary with each true gene name as key and a list of their alignment points as value
        :return pseudo_gene_pos: a dictionary with each pseudo gene name as key and a list of their alignment points as value
        """
        sam_path = self.config["data"]["sam"]["aligned_seq"].format(gene_group=self.gene_group, genome_ver=self.ref)
        samfile = pysam.AlignmentFile(sam_path, 'rb')

        points = set()

        # get all alignment pos we are interested
        for read in samfile.fetch():
            for read_offset, ref_offset, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                ref_gene_name, ref_exon_num, ref_chr, ref_start, _ = analyze_seq_name(read.reference_name)
                ref_base = ref_base.upper()
                ref_pos = ref_start + ref_offset

                read_gene_name, read_exon_num, read_chr, read_start, _ = analyze_seq_name(read.query_name) # with chr as prefix
                read_seq = read.get_forward_sequence()
                read_base = BASE_MATCH[read_seq[::-1][read_offset]] if read.is_reverse else read_seq[read_offset]
                read_base = read_base.upper()

                # check if hard-clipped at the begining of the read
                clipped_num = 0
                for (cigar_op, op_length) in read.cigartuples:
                    if cigar_op == 5: # 5 = hard-clipping
                        clipped_num += op_length
                    break
                read_pos = read_start + read_offset + clipped_num

                is_true_exception = (t := self.true_exception[ref_gene_name]) is not None and (ref_chr, ref_pos) in t.keys()
                is_pseudo_exception = (p := self.pseudo_exception[read_gene_name]) is not None and (read_chr, read_pos) in p.keys()
                if read_base != ref_base or is_true_exception or is_pseudo_exception:
                    point = (ref_gene_name, ref_chr, ref_pos)
                    points.add(point)
 
        # record info of the interesting alignment pos
        # we are doing this in two parts becauase of multiple reads alignment

        true_gene_pos = defaultdict(lambda: set())
        pseudo_gene_pos = defaultdict(lambda: set())
        alignment_relation = defaultdict(lambda: dict())
        for read in samfile.fetch():
            for read_offset, ref_offset, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                ref_gene_name, ref_exon_num, ref_chr, ref_start, _ = analyze_seq_name(read.reference_name)
                ref_base = ref_base.upper()
                ref_pos = ref_start + ref_offset

                read_gene_name, read_exon_num, read_chr, read_start, _ = analyze_seq_name(read.query_name) # with chr as prefix
                read_seq = read.get_forward_sequence()
                read_base = BASE_MATCH[read_seq[::-1][read_offset]] if read.is_reverse else read_seq[read_offset]
                read_base = read_base.upper()

                # check if hard-clipped at the begining of the read
                clipped_num = 0
                for (cigar_op, op_length) in read.cigartuples:
                    if cigar_op == 5: # 5 = hard-clipping
                        clipped_num += op_length
                    break
                read_pos = read_start + read_offset + clipped_num

                if (ref_gene_name, ref_chr, ref_pos) in points:
                    is_true_exception = (t := self.true_exception[ref_gene_name]) is not None and (ref_chr, ref_pos) in t.keys()
                    is_pseudo_exception = (p := self.pseudo_exception[read_gene_name]) is not None and (read_chr, read_pos) in p.keys()

                    if is_true_exception:
                        read_alt = base.upper() if (base:= t[(ref_chr, ref_pos)]) != '' else ref_base
                    elif is_pseudo_exception:
                        read_alt = base.upper() if (base:= p[(read_chr, read_pos)]) != '' else read_base
                    else:
                        read_alt = read_base

                    true_gene_pos[ref_gene_name].add((ref_chr, ref_pos, ref_base))
                    pseudo_gene_pos[read_gene_name].add((read_chr, read_pos, read_alt))

                    alignment_relation[(ref_chr, ref_pos, read_alt)][(read_chr, read_pos)] = [ref_gene_name, ref_exon_num, ref_base, \
                                                                                            read_gene_name, read_base, is_true_exception or is_pseudo_exception]

        true_gene_pos = {key: list(value) for key, value in true_gene_pos.items()}
        pseudo_gene_pos = {key: list(value) for key, value in pseudo_gene_pos.items()}

        # early return for parallelization / faster execution?
        aa_change_table = self.query_amino_acid_change(alignment_relation)
        for row in aa_change_table.itertuples():
            for key in alignment_relation[(row._1, row.POS, row.ALT)].keys():
                alignment_relation[(row._1, row.POS, row.ALT)][key].append(row.aa_change)
                
        # construct alignment table
        # will have problem during join if there are more than 1 true genes
        pos_col = ['exon', 'chr', 'pos', 'ref', 'alt', 'aa_change']
        ref_gene_names = set()
        read_gene_names = set()
        rows = defaultdict(lambda: defaultdict(lambda: list()))
        for (ref_chr, ref_pos, read_alt), values in alignment_relation.items():
            for (read_chr, read_pos), value in values.items():
                if len(value) == 6: # if pseudogene has the same ref base as true gene, fill aa change as "No Change"
                    value.append("No Change")
                ref_gene_name, ref_exon_num, ref_base, read_gene_name, read_base, is_exception, aa_change = value
                row = [is_exception, ref_exon_num, ref_chr, ref_pos, ref_base, read_chr, read_pos, read_base, read_alt, aa_change]
                rows[ref_gene_name][read_gene_name].append(row)
                ref_gene_names.add(ref_gene_name)
                read_gene_names.add(read_gene_name)

        tables = list()
        exception_col_names = list(itertools.product(['ratio'], ['exception'], ['is_exception']))
        for ref_gene_name, values in rows.items():
            true_col_names = list(itertools.product(['True'], [ref_gene_name], pos_col[:-2]))
            for read_gene_name, rs in values.items():
                pseudo_col_names = list(itertools.product(['Pseudo'], [read_gene_name], pos_col[1:]))
                col_names = pd.MultiIndex.from_tuples(exception_col_names + true_col_names + pseudo_col_names)
                table = pd.DataFrame(rs, columns=col_names)
                tables.append(table)
            alignment_table = functools.reduce(lambda left, right: left.merge(right, how='outer', \
                                                left_on=exception_col_names + true_col_names, 
                                                right_on=exception_col_names + true_col_names), tables)
        return alignment_table, true_gene_pos, pseudo_gene_pos, list(ref_gene_names), list(read_gene_names)

    def query_amino_acid_change(self, alignment_relation):
        vcf = create_vcf(alignment_relation, self.gene_group)
        with tempfile.NamedTemporaryFile(mode='r+', dir=self.tmp_dir) as vcf_file, \
             tempfile.NamedTemporaryFile(mode='r+', dir=self.tmp_dir) as nirvana_output:
            vcf_file.write(VCF_HEADER)
            vcf.to_csv(vcf_file, sep='\t', index=False, mode='a')
            vcf_file.seek(0)

            tool_path = self.config['tools']['nirvana_path']

            run_nirvana(tool_path, vcf_file.name, nirvana_output.name, self.ref)

            if len(self.true_gene) > 1:
                raise NotImplementedError('More than 1 true genes is not supported.')
            
            NM_num = self.config[self.gene_group][self.ref][self.true_gene[0]]["NM_number"]
            o = pathlib.Path(os.path.abspath(nirvana_output.name))
            aa_change = get_aa_change_from_nirvana(NM_num, vcf, o)
            aa_change_table = pd.concat([vcf, aa_change], axis=1)
            aa_change_table = aa_change_table[['#CHROM', 'POS', 'REF', 'ALT', 'aa_change']]
        return aa_change_table

    def __init__(self, test_file_dict, control_file_dict, ref, config, gene_group, output_path, ncpus=os.cpu_count(), tmp_dir=None, debug=False):
        """
        1. read config
        2. read sam & build alignment table
        3. initialize shared memeory array
        """
        self.read_config(config, gene_group, ref)
        self.output_path = output_path
        self.test_file_dict = test_file_dict
        self.control_file_dict = control_file_dict
        self.n_bam = len(self.test_file_dict) + len(self.control_file_dict)
        self.tmp_dir = tmp_dir
        self.ncpus = ncpus
        self.debug = debug

        if len(self.true_gene) > 1:
            raise NotImplementedError('algorithm for more than one true gene is not implemented!')

        # --------------- serialized
        self.alignment_table, self.true_gene_pos, self.pseudo_gene_pos, self.aligned_true_gene, self.aligned_pseudo_gene = self.construct_alignment_table()
        self.true_gene_region, self.pseudo_gene_region = self.get_gene_region()
        self.aligned_num_gene = len(self.aligned_true_gene) + len(self.aligned_pseudo_gene)
        
        self.control_gene_cov = np.zeros((self.n_bam, len(C.POSITIONS[self.ref]["GENES"])))
        self.gene_region_cov = np.zeros((self.n_bam, self.aligned_num_gene))
        
        self.true_gene_pos_len = sum(len(pos) for pos in self.true_gene_pos.values())
        self.pseudo_gene_pos_len = sum(len(pos) for pos in self.pseudo_gene_pos.values())
        
        self.true_table = np.zeros((self.n_bam, self.true_gene_pos_len, 4), dtype=np.int64)
        self.pseudo_table = np.zeros((self.n_bam, self.pseudo_gene_pos_len, 4), dtype=np.int64)

        # --------------- Raw Array
        # self.true_gene_region, self.pseudo_gene_region = self.get_gene_region()
        # self.alignment_table, self.true_gene_pos, self.pseudo_gene_pos = self.construct_alignment_table()

        # import ctypes
        # global control_gene_cov, gene_region_cov
        # control_gene_cov = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * len(C.POSITIONS[self.ref]["GENES"]))).reshape((self.n_bam, len(C.POSITIONS[self.ref]["GENES"])))
        # gene_region_cov = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * self.aligned_num_gene)).reshape((self.n_bam, self.aligned_num_gene))
        
        # self.true_gene_pos_len = sum(len(pos) for pos in self.true_gene_pos.values())
        # self.pseudo_gene_pos_len = sum(len(pos) for pos in self.pseudo_gene_pos.values())

        # global true_table, pseudo_table
        # true_table = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * self.true_gene_pos_len * 4)).reshape((self.n_bam, self.true_gene_pos_len, 4))
        # pseudo_table = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * self.pseudo_gene_pos_len * 4)).reshape((self.n_bam, self.pseudo_gene_pos_len, 4))

    def build_raw_table(self, idx):
        pos_col = ['exon', 'chr', 'pos', 'ref', 'alt', 'aa_change']
        ref_sub_pos_col = ['chr', 'pos', 'ref']
        alt_sub_pos_col = ['chr', 'pos', 'alt']

        alignment_table = self.alignment_table.copy()

        offset = 0
        for true_gene, pos in self.true_gene_pos.items():
            join_index = list(pd.MultiIndex.from_product([['True'], [true_gene], ref_sub_pos_col]))
            pos_info = pd.DataFrame(pos, columns=ref_sub_pos_col)
            cov_info = pd.DataFrame(self.true_table[idx, offset: offset + len(pos), :], columns=BASES)
            to_merged = pd.concat([pos_info, cov_info], axis=1, keys=[('True', true_gene), ('True', true_gene)])
            alignment_table = alignment_table.merge(to_merged, how='left', left_on=join_index, right_on=join_index)
            offset += len(pos)

        offset = 0
        for pseudo_gene, pos in self.pseudo_gene_pos.items():
            join_index = list(pd.MultiIndex.from_product([['Pseudo'], [pseudo_gene], alt_sub_pos_col]))
            pos_info = pd.DataFrame(pos, columns=alt_sub_pos_col)
            cov_info = pd.DataFrame(self.pseudo_table[idx, offset: offset + len(pos), :], columns=BASES)
            to_merged = pd.concat([pos_info, cov_info], axis=1, keys=[('Pseudo', pseudo_gene), ('Pseudo', pseudo_gene)])
            alignment_table = alignment_table.merge(to_merged, how='left', left_on=join_index, right_on=join_index)
            offset += len(pos)

        # re-order table columns
        ridx = [('ratio', 'exception', 'is_exception')] + \
                list(pd.MultiIndex.from_product([['True'], self.aligned_true_gene, pos_col[:-2] + BASES])) + \
                list(pd.MultiIndex.from_product([['Pseudo'], self.aligned_pseudo_gene, pos_col[1:] + BASES]))
        alignment_table = alignment_table[ridx]

        alignment_table.sort_values(by=[('True', self.aligned_true_gene[0], 'pos')], inplace=True, ignore_index=True)

        return alignment_table

    def calculate_sum(self, raw_table):
        true_idx = [pd.MultiIndex.from_product([['True'], [true_gene], BASES]) for true_gene in self.aligned_true_gene]
        pseudo_idx = [pd.MultiIndex.from_product([['Pseudo'], [pseudo_gene], BASES]) for pseudo_gene in self.aligned_pseudo_gene]

        sum_each_base = sum(np.nan_to_num(raw_table[col].to_numpy()) for col in true_idx + pseudo_idx)
        aligned_num_gene = sum((~np.isnan(raw_table[col].to_numpy())).astype(int)[:, 0].squeeze() for col in true_idx + pseudo_idx)

        sum_df = pd.DataFrame(sum_each_base, columns=pd.MultiIndex.from_product([["Sum"], ["bases"], BASES]))
        return sum_df, aligned_num_gene

    def calculate_ratio(self, raw_table, sum_df, scale_factor, aligned_num_gene):
        """
        1. get the refs per alignment point (from raw_table)
        2. get the bases sum per alignment point (from sum_df)
        3. calculate ratio (infer ratio)
        """
        def calculate_refs(raw_table):
            def get_refs(row):
                refs = defaultdict(lambda: list())
                for true_gene in self.aligned_true_gene:
                    ref = row[('True', true_gene, 'ref')]
                    refs[ref].append(true_gene)

                for pseudo_gene in self.aligned_pseudo_gene:
                    ref = row[('Pseudo', pseudo_gene, 'ref')]
                    refs[ref].append(pseudo_gene)

                if row[('ratio', 'exception', 'is_exception')]:
                    for pseudo_gene in self.aligned_seudo_gene:
                        alt = row[('Pseudo', pseudo_gene, 'alt')]
                        refs[alt] # make sure base in the key of refs, but not adding new elements to its value

                for key, value in refs.items():
                    if len(value) == 0:
                        refs[key].append('other')

                return pd.Series({' & '.join(gene_names): ref for ref, gene_names in refs.items()})
            refs = raw_table.apply(get_refs, axis=1)
            refs.columns = pd.MultiIndex.from_product([['ratio'], ['ref'], refs.columns])
            return refs

        def calculate_bases(refs, sum_df):
            vfunc = np.vectorize(get_base)

            ref_index = vfunc(refs.to_numpy()).astype(int)
            sum_arr = sum_df.to_numpy()
            extended_sum_arr = np.concatenate([sum_arr, np.full((len(sum_arr), 1), np.nan)], axis=1)

            bases_arr = np.take_along_axis(extended_sum_arr, ref_index, axis=1)
            bases = pd.DataFrame(bases_arr, columns=pd.MultiIndex.from_product([['ratio'], ['bases'], [ref[-1] for ref in refs.columns]]))
            return bases

        def calculate_sub_ratio(bases, scale_factor, num_genes):
            ratio = (bases / np.nansum(bases, axis=1)[:, np.newaxis] * scale_factor * num_genes[:, np.newaxis]).round(0)
            ratio.columns = pd.MultiIndex.from_product([['ratio'], ['Copy Number'], [base[-1] for base in bases.columns]])
            return ratio

        refs = calculate_refs(raw_table)
        bases = calculate_bases(refs, sum_df)
        ratios = calculate_sub_ratio(bases, scale_factor, aligned_num_gene)

        ratio_table = pd.concat([refs, bases, ratios], axis=1)
        return ratio_table

    def read_input(self, idx, file_path):
        bam_file = PseudoBam(file_path)
        self.true_table[idx, ...] = np.zeros((self.true_gene_pos_len, 4), dtype=np.int64)
        self.pseudo_table[idx, ...] = np.zeros((self.pseudo_gene_pos_len, 4), dtype=np.int64)

        self.control_gene_cov[idx, ...] = bam_file.get_cov_ranges(C.POSITIONS[self.ref]["GENES"])
        self.gene_region_cov[idx, ...] = bam_file.get_cov_ranges_(self.true_gene_region, self.pseudo_gene_region)

        offset = 0
        for true_gene, pos in self.true_gene_pos.items():
            for j, pos_info in enumerate(pos):
                self.true_table[idx, offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(pos_info[:2])).squeeze()
            offset += len(pos)

        offset = 0
        for pseudo_gene, pos in self.pseudo_gene_pos.items():
            for j, pos_info in enumerate(pos):
                self.pseudo_table[idx, offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(pos_info[:2])).squeeze()
            offset += len(pos)  

    def calculation_worker(self, idx):
        """ 
        format to pandas table (raw_summary, i.e. without sum & ratio) 

        """
        raw_table = self.build_raw_table(idx)
        sum_df, aligned_num_gene = self.calculate_sum(raw_table)
        ratio_df = self.calculate_ratio(raw_table, sum_df, self.scale_factor[idx], aligned_num_gene)
        raw_table_without_exception = raw_table.drop(columns=('ratio', 'exception', 'is_exception'))
        is_exception = raw_table[('ratio', 'exception', 'is_exception')]
        summary = pd.concat([is_exception, ratio_df, sum_df, raw_table_without_exception], axis=1)

        f = tempfile.NamedTemporaryFile(dir=self.tmp_dir, delete=False)
        summary.to_csv(f, sep='\t', index=False)
        return f.name

    def calculate(self):
        """
        1. read bam, get true/pseudo gene region (mean) coverage -> calculate scale factor
        2. read bam, get control gene (mean) coverage -> calculate scale factor
        3. read bam, get true/pseudo specific region coverage -> merge with alignment table -> make sum -> make ratio -> output
        """

        # --------------- read inputs

        # --------------- read inputs --------------- serialized
        # for idx, file_path in enumerate(itertools.chain(self.test_file_dict.values(), self.control_file_dict.values())):
        #     self.read_input(idx, file_path)

        # --------------- read inputs --------------- multithreaded
        # with ThreadPoolExecutor(max_workers=self.ncpus) as executor:
        #     cs = [executor.submit(self.read_input, idx=idx, file_path=file_path) for idx, file_path in enumerate(itertools.chain(self.test_file_dict.values(), self.control_file_dict.values()))]
        #     for future in concurrent.futures.as_completed(cs):
        #         future.result()

        # --------------- read inputs --------------- multiprocessing -> need return
        # with ThreadPoolExecutor(max_workers=self.ncpus) as executor:
        #     cs = [executor.submit(self.read_input, idx=idx, file_path=file_path) for idx, file_path in enumerate(itertools.chain(self.test_file_dict.values(), self.control_file_dict.values()))]
        #     for future in concurrent.futures.as_completed(cs):
        #         future.result()

        # --------------- read inputs --------------- multiprocessing - Copy On Write
        with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
            futures = {
                executor.submit(read_input_CopyOnWrite, 
                                file_path=file_path,
                                ref=self.ref, 
                                true_gene_region=self.true_gene_region,
                                pseudo_gene_region=self.pseudo_gene_region, 
                                true_gene_pos=self.true_gene_pos, 
                                pseudo_gene_pos=self.pseudo_gene_pos, 
                                true_genes=self.aligned_true_gene, 
                                pseudo_genes=self.aligned_pseudo_gene,
                                true_gene_pos_len=self.true_gene_pos_len,
                                pseudo_gene_pos_len=self.pseudo_gene_pos_len): idx
                for idx, file_path in enumerate(itertools.chain(self.test_file_dict.values(), self.control_file_dict.values()))
            }
            for future in concurrent.futures.as_completed(futures):
                idx = futures[future]
                self.control_gene_cov[idx, ...], self.gene_region_cov[idx, ...], self.true_table[idx, ...], self.pseudo_table[idx, ...] = future.result()

        # --------------- read inputs --------------- multiprocessing - Raw Array
        # with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
        #     futures = {
        #         executor.submit(read_input_RawArray, 
        #                         file_path=file_path,
        #                         ref=self.ref, 
        #                         true_gene_region=self.true_gene_region,
        #                         pseudo_gene_region=self.pseudo_gene_region, 
        #                         true_gene_pos=self.true_gene_pos, 
        #                         pseudo_gene_pos=self.pseudo_gene_pos, 
        #                         true_genes=self.true_gene, 
        #                         pseudo_genes=self.pseudo_gene,
        #                         idx=idx): idx
        #         for idx, file_path in enumerate(itertools.chain(self.test_file_dict.values(), self.control_file_dict.values()))
        #     }
        #     for future in concurrent.futures.as_completed(futures):
        #         idx = futures[future]
        #         future.result()

        # Others
        self.cov_ratio = self.gene_region_cov.sum(axis=1).reshape((self.n_bam, 1)) / self.control_gene_cov
        self.mean_cov_ratio = self.cov_ratio.sum(axis=0) / self.n_bam
        self.scale_factor = (self.cov_ratio / self.mean_cov_ratio).sum(axis=1) / len(C.POSITIONS[self.ref]["GENES"])

        # Raw Array
        # global control_gene_cov, gene_region_cov, true_table, pseudo_table
        # self.cov_ratio = gene_region_cov.sum(axis=1).reshape((self.n_bam, 1)) / control_gene_cov
        # self.mean_cov_ratio = self.cov_ratio.sum(axis=0) / self.n_bam
        # self.scale_factor = (self.cov_ratio / self.mean_cov_ratio).sum(axis=1) / len(C.POSITIONS[self.ref]["GENES"])

        # --------------- calculation --------------- serialized
        # fs = {
        #     case_name: self.calculation_worker(idx) for idx, case_name in enumerate(self.test_file_dict.keys())
        # }

        # --------------- calculation --------------- multithreaded
        # with ThreadPoolExecutor(max_workers=self.ncpus) as executor:           
        #     futures = {
        #         executor.submit(self.calculation_worker, idx=idx): case_name for idx, case_name in enumerate(self.test_file_dict.keys())
        #     }

        #     fs = {
        #         futures[future]: future.result() for future in concurrent.futures.as_completed(futures)
        #     }

        # --------------- calculation --------------- multiprocessing - Copy On Write
        # fs = {
        #     case_name: calculation_worker_CopyOnWrite(idx, self.alignment_table, \
        #                                     self.true_gene_pos, self.pseudo_gene_pos, \
        #                                     self.true_table, self.pseudo_table, \
        #                                     self.true_gene, self.pseudo_gene, \
        #                                     self.scale_factor, \
        #                                     self.aligned_num_gene, self.tmp_dir) for idx, case_name in enumerate(self.test_file_dict.keys()) 
        # }
        with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
            futures = {
                executor.submit(calculation_worker_CopyOnWrite, 
                                idx, self.alignment_table, \
                                self.true_gene_pos, self.pseudo_gene_pos, \
                                self.true_table, self.pseudo_table, \
                                self.aligned_true_gene, self.aligned_pseudo_gene, \
                                self.scale_factor, self.tmp_dir): case_name 
                for idx, case_name in enumerate(self.test_file_dict.keys())
            }
            fs = {
                futures[future]: future.result() for future in concurrent.futures.as_completed(futures)
            }

        # --------------- calculation --------------- multiprocessing - Raw Array
        # with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
        #     futures = {
        #         self.executor.submit(self.calculation_worker, idx=idx): case_name for idx, case_name in enumerate(self.test_file_dict.keys())
        #     }

        #     fs = {
        #         futures[future]: future.result() for future in concurrent.futures.as_completed(futures)
        #     }
        #     return fs, self.scale_factor

        if self.debug:
            test_df = pd.DataFrame(self.test_file_dict, index=['file path']).transpose()
            control_df = pd.DataFrame(self.control_file_dict, index=['file path']).transpose()
            index_df = pd.concat([test_df, control_df], axis=0, keys=['test', 'control'])
            index_df.rename_axis(index=['group', 'case name'], inplace=True)
            index_df.reset_index(inplace=True)
     
            # write out scale factors
            scale_factor_file_path = self.output_path / self.config["debug_output"].format(case_name=f"{self.gene_group}_scale_factor")
            scale_factor_df = pd.DataFrame(self.scale_factor, columns=['scale factor'])
            df = pd.concat([index_df, scale_factor_df], axis=1)
            df.set_index('case name', inplace=True)
            df.to_csv(scale_factor_file_path, sep='\t', )

            # write out control genes average coverage
            control_cov_file_path = self.output_path / self.config["debug_output"].format(case_name=f"{self.gene_group}_control_gene_cov")
            control_cov_df = pd.DataFrame(self.control_gene_cov, columns=C.POSITIONS[self.ref]["GENES"].keys())
            df = pd.concat([index_df, control_cov_df], axis=1)
            df.set_index('case name', inplace=True)
            df.to_csv(control_cov_file_path, sep='\t')

            # write out true genes and pseudogenes average coverage
            gene_cov_file_path = self.output_path / self.config["debug_output"].format(case_name=f"{self.gene_group}_gene_cov")
            gene_cov_df = pd.DataFrame(self.gene_region_cov, columns=self.aligned_true_gene + self.aligned_pseudo_gene)
            df = pd.concat([index_df, gene_cov_df], axis=1)
            df.set_index('case name', inplace=True)
            df.to_csv(gene_cov_file_path, sep='\t')

        return fs
