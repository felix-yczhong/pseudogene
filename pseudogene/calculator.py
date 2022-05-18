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
from weakref import ref

import numpy as np
import pandas as pd
import pysam

import smaca.constants as C
import pseudogene
from pseudogene.bam import PseudoBam
from pseudogene.utils import *

def read_input_RawArray(file_path, ref, true_gene_region, pseudo_gene_region, true_gene_pos, pseudo_gene_pos, true_genes, pseudo_genes, idx):
    bam_file = PseudoBam(file_path)

    global control_gene_cov, gene_region_cov, true_table, pseudo_table
    control_gene_cov[idx, ...] = bam_file.get_cov_ranges(C.POSITIONS[ref]["GENES"])
    gene_region_cov[idx, ...] = bam_file.get_cov_ranges_(true_gene_region, pseudo_gene_region)

    offset = 0
    for true_gene in true_genes:
        for j, range in enumerate(true_gene_pos[true_gene]):
            true_table[idx, offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(range[1:3])).squeeze()
        offset += len(true_gene_pos[true_gene])

    offset = 0
    for pseudo_gene in pseudo_genes:
        for j, range in enumerate(pseudo_gene_pos[pseudo_gene]):
            pseudo_table[idx, offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(range[:2])).squeeze()
        offset += len(pseudo_gene_pos[pseudo_gene])

def read_input_CopyOnWrite(file_path, ref, true_gene_region, pseudo_gene_region, true_gene_pos, pseudo_gene_pos, true_genes, pseudo_genes):
    """
    1. calculate H_ik
    2. calclulate c_ix (true/pseudo gene cov)

    3. true/pseudo gene cov (per case)
    4. ture/pseudo alignment points
    """
    bam_file = PseudoBam(file_path)

    num_gene = len(true_genes) + len(pseudo_genes)

    control_gene_cov = np.zeros((len(C.POSITIONS[ref]["GENES"])))
    gene_region_cov = np.zeros((num_gene))

    l1 = sum([len(true_gene_pos[true_gene]) for true_gene in true_genes])
    l2 = sum([len(pseudo_gene_pos[pseudo_gene]) for pseudo_gene in pseudo_genes])
    true_table = np.zeros((l1, 4), dtype=np.int64)
    pseudo_table = np.zeros((l2, 4), dtype=np.int64)

    control_gene_cov = bam_file.get_cov_ranges(C.POSITIONS[ref]["GENES"])
    gene_region_cov = bam_file.get_cov_ranges_(true_gene_region, pseudo_gene_region)

    offset = 0
    for true_gene in true_genes:
        for j, range in enumerate(true_gene_pos[true_gene]):
            true_table[offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(range[1:3])).squeeze()
        offset += len(true_gene_pos[true_gene])

    offset = 0
    for pseudo_gene in pseudo_genes:
        for j, range in enumerate(pseudo_gene_pos[pseudo_gene]):
            pseudo_table[offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(range[:2])).squeeze()
        offset += len(pseudo_gene_pos[pseudo_gene])
    return control_gene_cov, gene_region_cov, true_table, pseudo_table

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
            with open(self.config['data']['fasta']['true_seq'].format(gene=true_gene, genome_ver=ref)) as fasta:
                header = fasta.readline()
            ma = prog.match(header)
            if ma is None:
                raise ValueError(f"Cannot recognize {true_gene}_{self.ref}.fa header")
            chr, start, end = ma.group('chr', 'start', 'end')
            start, end = convert_coordiante_fasta_to_sam(int(start), int(end))
            return chr, start, end
        
        pattern = ">(?P<chr>chr[XY\d]+):(?P<start>[\d]+)-(?P<end>[\d]+)"
        prog = re.compile(pattern)

        true_gene_region = [read_fasta(true_gene, self.ref) for true_gene in self.true_gene]
        pseudo_gene_region = [read_fasta(pseudo_gene, self.ref) for pseudo_gene in self.pseudo_gene]
        return true_gene_region, pseudo_gene_region

    def get_aligned_pos(self, file_path, true_exception, pseudo_exception):
        samfile = pysam.AlignmentFile(file_path, 'rb')
        read_genomic_ranges = list()
        ref_genomic_ranges = list()

        ref_chr, ref_start, _ = analyze_ref_seq_name(samfile.get_reference_name(0))
        for read in samfile.fetch():
            for read_offset, ref_offset, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                read_seq = read.get_forward_sequence()

                exon_num, read_chr, read_start, _ = analyze_read_seq_name(read.query_name) # with chr as prefix
                read_base = BASE_MATCH[read_seq[::-1][read_offset]] if read.is_reverse else read_seq[read_offset]

                read_pos = read_start + read_offset
                ref_pos = ref_start + ref_offset

                if read_base.upper() != ref_base.upper() or \
                   true_exception is not None and (read_chr, read_pos) in true_exception.keys() or \
                   pseudo_exception is not None and (ref_chr, ref_pos) in pseudo_exception.keys():
                    read_genomic_pos = (exon_num, read_chr, read_pos, read_base)
                    read_genomic_ranges.append(read_genomic_pos)

                    if true_exception is not None and (read_chr, read_pos) in true_exception.keys():
                        alt = true_exception[(read_chr, read_pos)]
                    elif pseudo_exception is not None and (ref_chr, ref_pos) in pseudo_exception.keys():
                        alt = pseudo_exception[(ref_chr, ref_pos)]
                    else:
                        alt = ref_base.upper()
                    if alt == '':
                        alt = ref_base.upper()

                    ref_genomic_pos = (ref_chr, ref_pos, ref_base.upper(), alt)
                    ref_genomic_ranges.append(ref_genomic_pos)

        return read_genomic_ranges, ref_genomic_ranges

    def construct_alignment_table(self, alignments):
        # use DFS for connected component for many-many alignment relations
        pos_col = ['exon', 'chr', 'pos', 'ref']

        for true_gene in self.true_gene:
            true_to_pseudos_alignment = []
            true_gene_col_names = itertools.product(['true'], [true_gene], pos_col)

            for pseudo_gene in self.pseudo_gene:
                pseudo_gene_col_names = itertools.product(['pseudo'], [pseudo_gene], pos_col[1:])
                col_names = pd.MultiIndex.from_tuples(itertools.chain(true_gene_col_names, pseudo_gene_col_names))
                true_to_pseudo_alignment = pd.DataFrame([true_pos + pseudo_pos[:3] for (true_pos, pseudo_pos) in alignments[true_gene][pseudo_gene]], columns=col_names)
                true_to_pseudos_alignment.append(true_to_pseudo_alignment)
            join_col_names = list(pd.MultiIndex.from_product([['true'], [true_gene], pos_col]))
            true_to_pseudos_alignment = functools.reduce(lambda left, right: left.merge(right, how='outer', left_on=join_col_names, right_on=join_col_names), true_to_pseudos_alignment)
            return true_to_pseudos_alignment
        
    def get_alignment_info(self):
        """
        read bam files created by query_alignment() to get true gene(s) and pseudo gene(s) alignment information

        :param self:
        :return true_gene_pos: a dictionary with each true gene name as key and a list of their alignment points as value
        :return pseudo_gene_pos: a dictionary with each pseudo gene name as key and a list of their alignment points as value
        :return alignment_table: a pandas dataframe where each row represents how true gene(s) alignment point(s) align to pseudo gene(s) alignment point(s)
        """
        genomic_range_alignments = defaultdict(lambda: defaultdict(lambda: list()))
        true_genomic_ranges = defaultdict(lambda: set())
        pseudo_genomic_ranges = defaultdict(lambda: set())

        for true_gene in self.true_gene:
            for pseudo_gene in self.pseudo_gene:
                sam_path = self.config["data"]["sam"]["aligned_seq"].format(true_gene=true_gene, pseudo_gene=pseudo_gene, genome_ver=self.ref)

                read_genomic_ranges, ref_genomic_ranges = self.get_aligned_pos(sam_path, self.true_exception[true_gene], self.pseudo_exception[pseudo_gene])

                # collect genomic ranges
                true_genomic_ranges[true_gene] |= set(read_genomic_ranges)
                pseudo_genomic_ranges[pseudo_gene] |= set(ref_genomic_ranges)

                # build alignment info
                genomic_range_alignments[true_gene][pseudo_gene] = list(zip(read_genomic_ranges, ref_genomic_ranges))

        true_genomic_ranges = {key: sorted(list(value), reverse=True) for key, value in true_genomic_ranges.items()}
        pseudo_genomic_ranges = {key: sorted(list(value), reverse=True) for key, value in pseudo_genomic_ranges.items()}

        return genomic_range_alignments, true_genomic_ranges, pseudo_genomic_ranges

    def handle_nirvana(self, alignment_info):
        vcf = create_vcf(alignment_info, self.gene_group, self.true_gene, self.pseudo_gene)
        with tempfile.NamedTemporaryFile(mode='r+', dir=self.tmp_dir) as vcf_file, \
             tempfile.NamedTemporaryFile(mode='r+', dir=self.tmp_dir) as nirvana_output:
            vcf_file.write(VCF_HEADER)
            vcf.to_csv(vcf_file, sep='\t', index=False, mode='a')
            vcf_file.seek(0)

            tool_path = self.config['tools']['nirvana_path']

            run_nirvana(tool_path, vcf_file.name, nirvana_output.name, self.ref)

            # will have problem with multiple true genes
            NM_num = self.config[self.gene_group][self.ref][self.true_gene[0]]["NM_number"]
            o = pathlib.Path(os.path.abspath(nirvana_output.name))
            aa_change = get_aa_change_from_nirvana(NM_num, vcf, o)
        return aa_change

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

        self.executor = ThreadPoolExecutor(max_workers=ncpus)

        # self.get_gene_region_output_future = self.executor.submit(self.get_gene_region)
        self.true_gene_region, self.pseudo_gene_region = self.get_gene_region()
        alignment_info, self.true_gene_pos, self.pseudo_gene_pos = self.get_alignment_info()
        self.aa_change_future = self.executor.submit(self.handle_nirvana, alignment_info=alignment_info)

        self.alignment_table = self.construct_alignment_table(alignment_info)

        # aa_change = self.handle_nirvana(alignment_info)

        # self.alignment_table = pd.concat([alignment_table, aa_change], axis=1)

        self.control_gene_cov = np.zeros((self.n_bam, len(C.POSITIONS[self.ref]["GENES"])))
        self.gene_region_cov = np.zeros((self.n_bam, self.num_gene))

        l1 = sum([len(self.true_gene_pos[true_gene]) for true_gene in self.true_gene])
        l2 = sum([len(self.pseudo_gene_pos[pseudo_gene]) for pseudo_gene in self.pseudo_gene])
        self.true_table = np.zeros((self.n_bam, l1, 4), dtype=np.int64)
        self.pseudo_table = np.zeros((self.n_bam, l2, 4), dtype=np.int64)

        # ----------------------- RawArray
        # import ctypes
        # global control_gene_cov
        # control_gene_cov = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * len(C.POSITIONS[self.ref]["GENES"]))).reshape((self.n_bam, len(C.POSITIONS[self.ref]["GENES"])))
        # global gene_region_cov
        # gene_region_cov = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * self.num_gene)).reshape((self.n_bam, self.num_gene))
        # global true_table
        # true_table = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * l1 * 4)).reshape((self.n_bam, l1, 4))
        # global pseudo_table
        # pseudo_table = np.asarray(multiprocessing.RawArray(ctypes.c_double, self.n_bam * l2 * 4)).reshape((self.n_bam, l2, 4))

    def build_raw_table(self, idx):
        """
        combine alignment table and alignment position info
        """
        # need to handle multi-index
        i = ['exon', 'chr', 'pos', 'ref']
        alignment_table = self.alignment_table
        offset = 0
        for true_gene in self.true_gene:
            join_index = list(pd.MultiIndex.from_product([['true'], [true_gene], i]))
            pos_info = pd.DataFrame(self.true_gene_pos[true_gene], columns=i)
            cov_info = pd.DataFrame(self.true_table[idx, offset: offset + len(self.true_gene_pos[true_gene]), :], columns=BASES)
            to_merged = pd.concat([pos_info, cov_info], axis=1, keys=[('true', true_gene), ('true', true_gene)])
            alignment_table = alignment_table.merge(to_merged, how='left', left_on=join_index, right_on=join_index)
            offset += len(self.true_gene_pos[true_gene])

        offset = 0
        for pseudo_gene in self.pseudo_gene:
            join_index = list(pd.MultiIndex.from_product([['pseudo'], [pseudo_gene], i[1:]]))
            pos_info = pd.DataFrame([ele[0:3] for ele in self.pseudo_gene_pos[pseudo_gene]], columns=i[1:])
            cov_info = pd.DataFrame(self.pseudo_table[idx, offset: offset + len(self.pseudo_gene_pos[pseudo_gene]), :], columns=BASES)
            to_merged = pd.concat([pos_info, cov_info], axis=1, keys=[('pseudo', pseudo_gene), ('pseudo', pseudo_gene)])
            alignment_table = alignment_table.merge(to_merged, how='left', left_on=join_index, right_on=join_index)
            offset += len(self.pseudo_gene_pos[pseudo_gene])

        # re-organize table index
        idx = list(pd.MultiIndex.from_product([['true'], self.true_gene, i + BASES])) + list(pd.MultiIndex.from_product([['pseudo'], self.pseudo_gene, i[1:] + BASES])) + [('ratio', 'ratio', 'aa_change')]
        alignment_table = alignment_table[idx]

        return alignment_table

    def calculate_sum(self, raw_table):
        true_idx = [pd.MultiIndex.from_product([['true'], [true_gene], BASES]) for true_gene in self.true_gene]
        pseudo_idx = [pd.MultiIndex.from_product([['pseudo'], [pseudo_gene], BASES]) for pseudo_gene in self.pseudo_gene]

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
                    chr, pos = row[('true', true_gene, 'chr')], row[('true', true_gene, 'pos')]
                    if (e := true_exception[true_gene]) is not None and (base := e.get((chr, pos))) is not None:
                        refs[base] # make sure base in the key of refs, but not adding new elements to its value
                    ref = row[('true', true_gene, 'ref')]
                    refs[ref].append(true_gene)
                for pseudo_gene in self.pseudo_gene:
                    chr, pos = row[('pseudo', pseudo_gene, 'chr')], row[('pseudo', pseudo_gene, 'pos')]
                    if (e := pseudo_exception[pseudo_gene]) is not None and (base := e.get((chr, pos))) is not None:
                        refs[base] # make sure base in the key of refs, but not adding new elements to its value
                    ref = row[('pseudo', pseudo_gene, 'ref')]
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
        """
        1. calculate H_ik
        2. calclulate c_ix (true/pseudo gene cov)

        3. true/pseudo gene cov (per case)
        4. ture/pseudo alignment points
        """
        bam_file = PseudoBam(file_path)
        l1 = sum([len(self.true_gene_pos[true_gene]) for true_gene in self.true_gene])
        l2 = sum([len(self.pseudo_gene_pos[pseudo_gene]) for pseudo_gene in self.pseudo_gene])
        self.true_table[idx, ...] = np.zeros((l1, 4), dtype=np.int64)
        self.pseudo_table[idx, ...] = np.zeros((l2, 4), dtype=np.int64)

        self.control_gene_cov[idx, ...] = bam_file.get_cov_ranges(C.POSITIONS[self.ref]["GENES"])
        self.gene_region_cov[idx, ...] = bam_file.get_cov_ranges_(self.true_gene_region, self.pseudo_gene_region)

        offset = 0
        for true_gene in self.true_gene:
            for j, range in enumerate(self.true_gene_pos[true_gene]):
                self.true_table[idx, offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(range[1:3])).squeeze()
            offset += len(self.true_gene_pos[true_gene])

        offset = 0
        for pseudo_gene in self.pseudo_gene:
            for j, range in enumerate(self.pseudo_gene_pos[pseudo_gene]):
                self.pseudo_table[idx, offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(range[:2])).squeeze()
            offset += len(self.pseudo_gene_pos[pseudo_gene])

    def calculation_worker(self, idx):
        """ 
        format to pandas table (raw_summary, i.e. without sum & ratio) 
        
        """
        raw_table = self.build_raw_table(idx)
        sum_df = self.calculate_sum(raw_table)
        ratio_df = self.calculate_ratio(raw_table, sum_df, self.scale_factor[idx])
        summary = pd.concat([ratio_df, sum_df, raw_table], axis=1)

        # re-order table column
        aa_change_col = summary.pop(('ratio', 'ratio', 'aa_change'))
        summary.insert(0, aa_change_col.name, aa_change_col)

        # re-order table row
        reorder_idx = list(pd.MultiIndex.from_product([['true'], [self.true_gene[0]], ['exon', 'chr', 'pos']]))
        # summary_exception = summary[]
        summary = summary.sort_values(by=reorder_idx)

        f = tempfile.NamedTemporaryFile(dir=self.tmp_dir, delete=False)
        summary.to_csv(f, sep='\t', index=False)
        return f.name

    def calculate(self):
        """
        1. read bam, get true/pseudo gene region (mean) coverage -> calculate scale factor
        2. read bam, get control gene (mean) coverage -> calculate scale factor
        3. read bam, get true/pseudo specific region coverage -> merge with alignment table -> make sum -> make ratio -> output
        """

        # ---------------------- multiprocessing - copy-on-wirte
        # with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
        #     futures = dict()
        #     for idx, file_path in enumerate(self.file_dict.values()):
        #         future = executor.submit(read_input_CopyOnWrite, file_path=file_path,
        #                                                             ref=self.ref, 
        #                                                             true_gene_region=self.true_gene_region,
        #                                                             pseudo_gene_region=self.pseudo_gene_region, 
        #                                                             true_gene_pos=self.true_gene_pos, 
        #                                                             pseudo_gene_pos=self.pseudo_gene_pos, 
        #                                                             true_genes=self.true_gene, 
        #                                                             pseudo_genes=self.pseudo_gene)
        #         futures[future] = idx
        #     for future in concurrent.futures.as_completed(futures):
        #         idx = futures[future]
        #         self.control_gene_cov[idx, ...], self.gene_region_cov[idx, ...], self.true_table[idx, ...], self.pseudo_table[idx, ...] = future.result()

        # ---------------------- multiprocessing - RawArray
        # with ProcessPoolExecutor(max_workers=self.ncpus) as executor:
        #     futures = dict()
        #     for idx, file_path in enumerate(self.file_dict.values()):
        #         future = executor.submit(read_input_RawArray, file_path=file_path,
        #                                                             ref=self.ref, 
        #                                                             true_gene_region=self.true_gene_region,
        #                                                             pseudo_gene_region=self.pseudo_gene_region, 
        #                                                             true_gene_pos=self.true_gene_pos, 
        #                                                             pseudo_gene_pos=self.pseudo_gene_pos, 
        #                                                             true_genes=self.true_gene, 
        #                                                             pseudo_genes=self.pseudo_gene,
        #                                                             idx=idx)
        #         futures[future] = idx
        #     for future in concurrent.futures.as_completed(futures):
        #         future.result()

        # ---------------------- serialized 
        # for idx, file_path in enumerate(self.file_dict.values()):
        #     self.read_input(idx, file_path)

        # ---------------------- multithreaded
        # with ThreadPoolExecutor(max_workers=self.ncpus) as executor:
        fs = [self.executor.submit(self.read_input, idx=idx, file_path=file_path) for idx, file_path in enumerate(self.file_dict.values())]
        for future in concurrent.futures.as_completed(fs):
            future.result()

        # cov_ratio = average target gene(s) coverage over average control genes coverage
        self.cov_ratio = self.gene_region_cov.sum(axis=1).reshape((self.n_bam, 1)) / self.control_gene_cov
        # self.std_k = np.std(self.H_ik, axis=0)
        # self.std_i = np.std(self.H_ik, axis=1)
        self.mean_cov_ratio = self.cov_ratio.sum(axis=0) / self.n_bam
        self.scale_factor = (self.cov_ratio / self.mean_cov_ratio).sum(axis=1) / len(C.POSITIONS[self.ref]["GENES"])

        # ---------------------- multiprocessing - RawArray
        # self.cov_ratio = gene_region_cov.sum(axis=1).reshape((self.n_bam, 1)) / control_gene_cov
        # self.mean_cov_ratio = self.cov_ratio.sum(axis=0) / self.n_bam
        # self.scale_factor = (self.cov_ratio / self.mean_cov_ratio).sum(axis=1) / len(C.POSITIONS[self.ref]["GENES"])
        # self.true_table = true_table
        # self.pseudo_table = pseudo_table

        for future in concurrent.futures.as_completed([self.aa_change_future]):
            aa_change = self.aa_change_future.result()
        self.alignment_table = pd.concat([self.alignment_table, aa_change], axis=1)
        # prioritize exception points in alignment table

        control_gene_cov = pd.DataFrame(self.control_gene_cov, index=self.file_dict.values(), columns=C.POSITIONS[self.ref]["GENES"].keys())
        control_gene_cov.to_csv(f'{self.gene_group}_control_gene_cov.csv')
        gene_region_cov = pd.DataFrame(self.gene_region_cov, index=self.file_dict.values(), columns=self.true_gene + self.pseudo_gene)
        gene_region_cov.to_csv(f'{self.gene_group}_gene_region_cov.csv')
        # np.savetxt(f'{self.gene_group}_control_gene_cov.csv', self.control_gene_cov,  fmt='%3f', delimiter=',')
        # np.savetxt(f'{self.gene_group}_gene_region_cov.csv', self.gene_region_cov, fmt='%3f', delimiter=',')

        # with ThreadPoolExecutor(max_workers=self.ncpus) as executor:
        fs = {
            self.executor.submit(self.calculation_worker, idx=idx): case_name for idx, case_name in enumerate(self.file_dict.keys())
        }
        caseName_to_outputFilePath_mapping = {
            fs[future]: future.result() for future in concurrent.futures.as_completed(fs)
        }
        return caseName_to_outputFilePath_mapping, self.scale_factor