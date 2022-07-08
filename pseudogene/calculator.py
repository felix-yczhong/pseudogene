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

        true_gene_region = [read_fasta(true_gene, self.ref) for true_gene in self.true_gene]
        pseudo_gene_region = [read_fasta(pseudo_gene, self.ref) for pseudo_gene in self.pseudo_gene]
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
                read_seq = read.get_forward_sequence()

                read_gene_name, read_exon_num, read_chr, read_start, _ = analyze_seq_name(read.query_name) # with chr as prefix
                read_base = BASE_MATCH[read_seq[::-1][read_offset]] if read.is_reverse else read_seq[read_offset]

                read_pos = read_start + read_offset
                ref_pos = ref_start + ref_offset
                # if ref_pos == 155235203:
                #     print(ref_start, ref_offset)
                #     print(read_start, read_offset)

                is_true_exception = (t := self.true_exception[ref_gene_name]) is not None and (ref_chr, ref_pos) in t.keys()
                is_pseudo_exception = (p := self.pseudo_exception[read_gene_name]) is not None and (read_chr, read_pos) in p.keys()
                if read_base.upper() != ref_base.upper() or is_true_exception or is_pseudo_exception:
                    point = (ref_gene_name, ref_chr, ref_pos)
                    points.add(point)
 
        # record info of the interesting alignment pos
        # we are doing this in two parts becauase of multiple reads alignment
        true_gene_pos = defaultdict(lambda: list())
        pseudo_gene_pos = defaultdict(lambda: list())
        alignment_relation = dict()
        for read in samfile.fetch():
            for read_offset, ref_offset, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                ref_gene_name, ref_exon_num, ref_chr, ref_start, _ = analyze_seq_name(read.reference_name)
                read_seq = read.get_forward_sequence()

                read_gene_name, read_exon_num, read_chr, read_start, _ = analyze_seq_name(read.query_name) # with chr as prefix
                read_base = BASE_MATCH[read_seq[::-1][read_offset]] if read.is_reverse else read_seq[read_offset]

                read_pos = read_start + read_offset
                ref_pos = ref_start + ref_offset

                if (ref_gene_name, ref_chr, ref_pos) in points:
                    is_true_exception = (t := self.true_exception[ref_gene_name]) is not None and (ref_chr, ref_pos) in t.keys()
                    is_pseudo_exception = (p := self.pseudo_exception[read_gene_name]) is not None and (read_chr, read_pos) in p.keys()

                    # if is_true_exception:
                    #     ref_alt = a.upper() if (a:= t[(ref_chr, ref_pos)]) != '' else ref_base.upper()
                    # else:
                    #     ref_alt = ref_base.upper()
                    # if is_pseudo_exception:
                    #     read_alt = a.upper() if (a:= p[(read_chr, read_pos)]) != '' else read_base.upper()
                    # else:
                    #     read_alt = read_base.upper()
                    ref_alt = ref_base.upper()
                    if is_true_exception:
                        read_alt = a.upper() if (a:= t[(ref_chr, ref_pos)]) != '' else ref_base.upper()
                    elif is_pseudo_exception:
                        read_alt = a.upper() if (a:= p[(read_chr, read_pos)]) != '' else read_base.upper()
                    else:
                        read_alt = read_base.upper()

                    true_gene_pos[ref_gene_name].append((ref_chr, ref_pos, ref_alt))
                    pseudo_gene_pos[read_gene_name].append((read_chr, read_pos, read_alt))

                    alignment_relation[(ref_chr, ref_pos, read_alt)] = [ref_gene_name, ref_exon_num, ref_base.upper(), ref_alt, \
                                                                        read_gene_name, read_pos, read_base.upper()]

        # early return for parallelization / faster execution?
        aa_change_table = self.query_amino_acid_change(alignment_relation)
        for row in aa_change_table.itertuples():
            alignment_relation[(row._1, row.POS, row.ALT)].append(row.aa_change)

        # construct alignment table
        # will have problem during join if there are more than 1 true genes
        pos_col = ['exon', 'chr', 'pos', 'ref', 'alt', 'aa_change']
        rows = defaultdict(lambda: defaultdict(lambda: list()))
        for (ref_chr, ref_pos, read_alt), \
            (ref_gene_name, ref_exon_num, ref_base, ref_alt, \
            read_gene_name, read_pos, read_base, aa_change) in alignment_relation.items():

            row = [ref_exon_num, ref_chr, ref_pos, ref_base, ref_alt, read_chr, read_pos, read_base, read_alt, aa_change]
            rows[ref_gene_name][read_gene_name].append(row)

        tables = list()
        for ref_gene_name, v in rows.items():
            for read_gene_name, r in v.items():
                true_col_names = itertools.product(['True'], [ref_gene_name], pos_col[:-1])
                pseudo_col_names = itertools.product(['Pseudo'], [read_gene_name], pos_col[1:])
                col_names = pd.MultiIndex.from_tuples(itertools.chain(true_col_names, pseudo_col_names))
                table = pd.DataFrame(r, columns=col_names)
                tables.append(table)
        alignment_table = functools.reduce(lambda left, right: left.merge(right, how='outer', left_on=true_col_names, right_on=true_col_names), tables)
        return alignment_table, true_gene_pos, pseudo_gene_pos

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

    def __init__(self, file_dict, ref, config, gene_group, ncpus=os.cpu_count(), tmp_dir=None, debug=False):
        """
        1. read config
        2. read sam & build alignment table
        3. initialize shared memeory array
        """
        self.read_config(config, gene_group, ref)
        self.file_dict = file_dict
        self.n_bam = len(self.file_dict)
        self.tmp_dir = tmp_dir
        self.debug = debug

        if len(self.true_gene) != 1:
            raise NotImplementedError('algorithm for more than one true gene is not implemented!')
        if len(self.pseudo_gene) != 1:
            raise NotImplementedError('algorithm for more than one pseudo gene is not fully implemented!')

        # self.executor = ProcessPoolExecutor(max_workers=ncpus)

        # --------------- serialized
        self.true_gene_region, self.pseudo_gene_region = self.get_gene_region()
        self.alignment_table, self.true_gene_pos, self.pseudo_gene_pos = self.construct_alignment_table()
        
        self.control_gene_cov = np.zeros((self.n_bam, len(C.POSITIONS[self.ref]["GENES"])))
        self.gene_region_cov = np.zeros((self.n_bam, self.num_gene))
        
        self.true_gene_pos_len = sum(len(pos) for pos in self.true_gene_pos.values())
        self.pseudo_gene_pos_len = sum(len(pos) for pos in self.pseudo_gene_pos.values())
        
        self.true_table = np.zeros((self.n_bam, self.true_gene_pos_len, 4), dtype=np.int64)
        self.pseudo_table = np.zeros((self.n_bam, self.pseudo_gene_pos_len, 4), dtype=np.int64)

        # --------------- Raw Array
        # self.true_gene_region, self.pseudo_gene_region = self.get_gene_region()
        # self.alignment_table, self.true_gene_pos, self.pseudo_gene_pos = self.build_alignment_table()

        # import ctypes
        # global control_gene_cov
        # control_gene_cov = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * len(C.POSITIONS[self.ref]["GENES"]))).reshape((self.n_bam, len(C.POSITIONS[self.ref]["GENES"])))
        # global gene_region_cov
        # gene_region_cov = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * self.num_gene)).reshape((self.n_bam, self.num_gene))
        
        # self.true_gene_pos_len = sum(len(pos) for pos in self.true_gene_pos.values())
        # self.pseudo_gene_pos_len = sum(len(pos) for pos in self.pseudo_gene_pos.values())

        # global true_table
        # true_table = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * self.true_gene_pos_len * 4)).reshape((self.n_bam, l1, 4))
        # global pseudo_table
        # pseudo_table = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * self.pseudo_gene_pos_len * 4)).reshape((self.n_bam, l2, 4))

    def build_raw_table(self, idx):
        pos_col = ['exon', 'chr', 'pos', 'ref', 'alt', 'aa_change']
        sub_pos_col = ['chr', 'pos', 'alt']

        alignment_table = self.alignment_table.copy()

        offset = 0
        for true_gene, pos in self.true_gene_pos.items():
            join_index = list(pd.MultiIndex.from_product([['True'], [true_gene], sub_pos_col]))
            pos_info = pd.DataFrame(pos, columns=sub_pos_col)
            cov_info = pd.DataFrame(self.true_table[idx, offset: offset + len(pos), :], columns=BASES)
            to_merged = pd.concat([pos_info, cov_info], axis=1, keys=[('True', true_gene), ('True', true_gene)])
            alignment_table = alignment_table.merge(to_merged, how='left', left_on=join_index, right_on=join_index)
            offset += len(pos)

        offset = 0
        for pseudo_gene, pos in self.pseudo_gene_pos.items():
            join_index = list(pd.MultiIndex.from_product([['Pseudo'], [pseudo_gene], sub_pos_col]))
            pos_info = pd.DataFrame(pos, columns=sub_pos_col)
            cov_info = pd.DataFrame(self.pseudo_table[idx, offset: offset + len(pos), :], columns=BASES)
            to_merged = pd.concat([pos_info, cov_info], axis=1, keys=[('Pseudo', pseudo_gene), ('Pseudo', pseudo_gene)])
            alignment_table = alignment_table.merge(to_merged, how='left', left_on=join_index, right_on=join_index)
            offset += len(pos)

        # re-organize table index
        ridx = list(pd.MultiIndex.from_product([['True'], self.true_gene, pos_col[:-1] + BASES])) + \
               list(pd.MultiIndex.from_product([['Pseudo'], self.pseudo_gene, pos_col[1:] + BASES]))
        alignment_table = alignment_table[ridx]
        alignment_table.sort_values(by=[('True', self.true_gene[0], 'pos')], inplace=True)

        return alignment_table

    def calculate_sum(self, raw_table):
        true_idx = [pd.MultiIndex.from_product([['True'], [true_gene], BASES]) for true_gene in self.true_gene]
        pseudo_idx = [pd.MultiIndex.from_product([['Pseudo'], [pseudo_gene], BASES]) for pseudo_gene in self.pseudo_gene]

        sum_each_base = sum(raw_table[col].to_numpy() for col in true_idx + pseudo_idx)
        sum_df = pd.DataFrame(sum_each_base, columns=pd.MultiIndex.from_product([["Sum"], ["bases"], BASES]))
        return sum_df

    def calculate_ratio(self, raw_table, sum_df, scale_factor):
        """
        1. get the refs per alignment point (from raw_table)
        2. get the bases sum per alignment point (from sum_df)
        3. calculate ratio (infer ratio)
        """
        def calculate_refs(raw_table, true_exception, pseudo_exception):
            def get_refs(row):
                refs = defaultdict(lambda: list())
                for true_gene in self.true_gene:
                    chr, pos = row[('True', true_gene, 'chr')], row[('True', true_gene, 'pos')]
                    if (e := true_exception[true_gene]) is not None and (base := e.get((chr, pos))) is not None:
                        refs[base] # make sure base in the key of refs, but not adding new elements to its value
                    ref = row[('True', true_gene, 'ref')]
                    refs[ref].append(true_gene)
                for pseudo_gene in self.pseudo_gene:
                    chr, pos = row[('Pseudo', pseudo_gene, 'chr')], row[('Pseudo', pseudo_gene, 'pos')]
                    if (e := pseudo_exception[pseudo_gene]) is not None and (base := e.get((chr, pos))) is not None:
                        refs[base] # make sure base in the key of refs, but not adding new elements to its value
                    ref = row[('Pseudo', pseudo_gene, 'ref')]
                    refs[ref].append(pseudo_gene)

                for key, value in refs.items():
                    if len(value) == 0:
                        refs[key].append('other')

                return pd.Series({' & '.join(gene_names): ref for ref, gene_names in refs.items()})
            refs = raw_table.apply(get_refs, axis=1)
            refs.columns = pd.MultiIndex.from_product([['ratio'], ['ref'], refs.columns])
            return refs

        def calculate_bases(refs, sum_df):
            def get_bases(row):
                refs = row[('ratio', 'ref')].tolist()
                bases = [row[('Sum', 'bases', ref)] if isinstance(ref, str) else None for ref in refs]
                return pd.Series(bases)
            bases = pd.concat([refs, sum_df], axis=1).apply(get_bases, axis=1)
            bases.columns = pd.MultiIndex.from_product([['ratio'], ['bases'], [ref[-1] for ref in refs.columns]])
            return bases

        def calculate_sub_ratio(bases, scale_factor, num_genes):
            def get_ratio(row):
                bases = row[('ratio', 'bases')].tolist()
                sum_base = sum(base for base in bases if not np.isnan(base))
                if sum_base == 0:
                    return pd.Series(bases)

                scaled_ratio = [(base / sum_base) * scale_factor * num_genes * PLOIDY for base in bases] # num_genes need to change in multiple true/pseudogene scenario
                # min_ele = min([base for base in scaled_ratio if not np.isnan(base)])
                # multiplier = round(min_ele) / min_ele if min_ele != 0 and round(min_ele) != 0 else 1 # what if 2:0 (min_ele != 0 but rounds to 0)
                # predicted_ratio = [ele if np.isnan(ele) else round(ele * multiplier) for ele in scaled_ratio]
                predicted_ratio = [ratio if np.isnan(ratio) else round(ratio) for ratio in scaled_ratio]
                
                return pd.Series(predicted_ratio)
            ratio = bases.apply(get_ratio, axis=1)
            ratio.columns = pd.MultiIndex.from_product([['ratio'], ['Copy Number'], [base[-1] for base in bases.columns]])
            return ratio

        refs = calculate_refs(raw_table, self.true_exception, self.pseudo_exception)
        bases = calculate_bases(refs, sum_df)
        ratios = calculate_sub_ratio(bases, scale_factor, self.num_gene)

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
        sum_df = self.calculate_sum(raw_table)
        ratio_df = self.calculate_ratio(raw_table, sum_df, self.scale_factor[idx])
        summary = pd.concat([ratio_df, sum_df, raw_table], axis=1)

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
        for idx, file_path in enumerate(self.file_dict.values()):
            self.read_input(idx, file_path)

        # --------------- read inputs --------------- multithreaded
        # with ThreadPoolExecutor(max_worker=self.ncpus) as executor:
        #     cs = [executor.submit(self.read_input, idx=idx, file_path=file_path) for idx, file_path in enumerate(self.file_dict.values())]
        #     for future in concurrent.futures.as_completed(cs):
        #         future.result()

        # --------------- read inputs --------------- multiprocessing - Copy On Write
        # with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
        #     futures = {
        #         executor.submit(read_input_CopyOnWrite, 
        #                         file_path=file_path,
        #                         ref=self.ref, 
        #                         true_gene_region=self.true_gene_region,
        #                         pseudo_gene_region=self.pseudo_gene_region, 
        #                         true_gene_pos=self.true_gene_pos, 
        #                         pseudo_gene_pos=self.pseudo_gene_pos, 
        #                         true_genes=self.true_gene, 
        #                         pseudo_genes=self.pseudo_gene): idx
        #         for idx, file_path in enumerate(self.file_dict.values())
        #     }

        #     for future in concurrent.futures.as_completed(futures):
        #         idx = futures[future]
        #         self.control_gene_cov[idx, ...], self.gene_region_cov[idx, ...], self.true_table[idx, ...], self.pseudo_table[idx, ...] = future.result()

        # --------------- read inputs --------------- multiprocessing - Raw Array
        # with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
        #     futures = [
        #         executor.submit(read_input_RawArray, 
        #                         file_path=file_path,
        #                         ref=self.ref, 
        #                         true_gene_region=self.true_gene_region,
        #                         pseudo_gene_region=self.pseudo_gene_region, 
        #                         true_gene_pos=self.true_gene_pos, 
        #                         pseudo_gene_pos=self.pseudo_gene_pos, 
        #                         true_genes=self.true_gene, 
        #                         pseudo_genes=self.pseudo_gene,
        #                         idx=idx)
        #         for idx, file_path in enumerate(self.file_dict.values())
        #     ]
        #     for future in concurrent.futures.as_completed(futures):
        #         futures[future]

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
        fs = {
            case_name: self.calculation_worker(idx) for idx, case_name in enumerate(self.file_dict.keys())
        }

        # --------------- calculation --------------- multithreaded
        # with ThreadPoolExecutor(max_worker=self.ncpus) as executor:           
        #     futures = {
        #         executor.submit(self.calculation_worker, idx=idx): case_name for idx, case_name in enumerate(self.file_dict.keys())
        #     }

        #     fs = {
        #         futures[future]: future.result() for future in concurrent.futures.as_completed(futures)
        #     }

        # --------------- calculation --------------- multiprocessing - Copy On Write
        # with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
        #     futures = {
        #         executor.submit(calculation_CopyOnWrite, 
        #                         idx=idx,
        #                         true_genes=self.true_gene,
        #                         pseudo_genes=self.pseudo_gene,
        #                         alignment_table=self.alignment_table,
        #                         true_gene_pos=self.true_gene_pos,
        #                         pseudo_gene_pos=self.pseudo_gene_pos,
        #                         scale_factor=self.scale_factor,
        #                         num_gene=self.num_gene,
        #                         true_exception=self.true_exception,
        #                         pseudo_exception=self.pseudo_exception,
        #                         tmp_dir=self.tmp_dir): case_name 
        #         for idx, case_name in enumerate(self.file_dict.keys())
        #     }
        #     fs = {
        #         futures[future]: future.result() for future in concurrent.futures.as_completed(futures)
        #     }

        # --------------- calculation --------------- multiprocessing - Raw Array
        # with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
        #     futures = {
        #         self.executor.submit(self.calculation_worker, idx=idx): case_name for idx, case_name in enumerate(self.file_dict.keys())
        #     }

        #     fs = {
        #         futures[future]: future.result() for future in concurrent.futures.as_completed(futures)
        #     }
        #     return fs, self.scale_factor

        # control_gene_cov = pd.DataFrame(self.control_gene_cov, index=self.file_dict.values(), columns=C.POSITIONS[self.ref]["GENES"].keys())
        # control_gene_cov.to_csv(f'{self.gene_group}_control_gene_cov.csv')
        # gene_region_cov = pd.DataFrame(self.gene_region_cov, index=self.file_dict.values(), columns=self.true_gene + self.pseudo_gene)
        # gene_region_cov.to_csv(f'{self.gene_group}_gene_region_cov.csv')
        # np.savetxt(f'{self.gene_group}_control_gene_cov.csv', self.control_gene_cov,  fmt='%3f', delimiter=',')
        # np.savetxt(f'{self.gene_group}_gene_region_cov.csv', self.gene_region_cov, fmt='%3f', delimiter=',')

        return fs, self.scale_factor