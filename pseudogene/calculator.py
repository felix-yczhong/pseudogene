import subprocess
import concurrent
import tempfile
import pathlib
from concurrent.futures import ThreadPoolExecutor
from itertools import product

import numpy as np
import pandas as pd

from pseudogene.sma_extended import SmaCalculatorExtended
from pseudogene.bam import PseudoBam
from pseudogene.utils import *

class PseudoGeneCalculator(SmaCalculatorExtended):
    def __init__(self, bam_list, case_names, ref, config, gene, n_jobs=1, tmp_dir=None, debug=False):
        super().__init__(bam_list, ref, config, gene, n_jobs)

        self.case_names = case_names
        self.refs = self.get_refs()
        self.debug = debug
        self.tmp_dir = pathlib.Path(tempfile.mkdtemp()) if tmp_dir is None else tmp_dir

    def get_ref(self, row):
        c, start, stop = row[['Chrom', 'Start', 'End']]
        ref = subprocess.run(['samtools', 'faidx', self.fasta_loc, f'chr{c}:{start + 1}-{stop}'], stdout=subprocess.PIPE).stdout
        ref = pd.Series(list(chr(n) for n in ref if chr(n) in BASES), name="Ref")
        return ref

    def get_refs(self):
        true_ref = self.seq_align_map['True Gene'][self.true_gene].apply(self.get_ref, axis=1, result_type='expand').T
        pseudo_ref = {name: df[name].apply(self.get_ref, axis=1, result_type='expand').T for name, df in self.seq_align_map['Pseudo Gene'].groupby(level=0, axis=1)}
        
        refs = list()
        for i in range(self.num_seq):
            true = pd.concat({"Ref": true_ref[i].dropna(axis=0, how='any')}, axis=1)
            true_df = pd.concat({self.true_gene: true}, axis=1)
            pseudo = {name: pd.concat({"Ref": pseudo_ref[name][i].dropna(axis=0, how='any')}, axis=1) for name in self.pseudo_genes}
            pseudo_df = pd.concat(pseudo, axis=1)

            refs.append(pd.concat({"True Gene": true_df, "Pseudo Gene": pseudo_df}, axis=1))

        return refs

    def get_cov_detail(self, row, bam_file, key):
        c, start, stop = bam_file.get_genomic_range(row[['Chrom', 'Start', 'End']])
        pos = pd.Series(range(start + 1, stop + 1), name="Pos")
        chrom = pd.Series(np.full(stop - start, c), name="Chrom")
        cov_detail = pd.DataFrame(bam_file.get_cov([c, start, stop]), columns=BASES + ["Sum"])
        to_concat = [chrom, pos, cov_detail]
        df = pd.concat(to_concat, axis=1, keys=[key] * len(to_concat))
        return df

    def calculate_df(self, row, bam_file: PseudoBam):      
        cov_detail_true = self.get_cov_detail(row['True Gene'][self.true_gene], bam_file, key=self.true_gene)
        cov_detail_pseudoes = pd.concat([self.get_cov_detail(row['Pseudo Gene'][gene], bam_file, key=gene) for gene in self.pseudo_genes], axis=1)
        df = pd.concat([cov_detail_true, cov_detail_pseudoes], axis=1, keys=["True Gene", "Pseudo Gene"])
        return df

    def filter(self, row, df):
        c, start, stop = convert_coordinate_zero_to_one(row["True Gene"][self.true_gene][['Chrom', 'Start', 'End']])
        return df[df[('True Gene', self.true_gene, 'Pos')].between(start, stop, inclusive='left')]

    def call_summary(self, df):
        filter = (df[('True Gene', self.true_gene, 'Ref')] != df[('Pseudo Gene', self.pseudo_genes[0], 'Ref')])
        if self.force.empty:
            return df[filter]
        else:
            force_output = pd.concat(self.force.apply(self.filter, args=[df], axis=1).tolist(), axis=0)
            return pd.concat([df[filter], force_output], axis=0).drop_duplicates()

    def calculate_refs(self, intermediate_summary):
        def _calculate_ref(row):
            ref_true = row[('True Gene', self.true_gene, 'Ref')]
            ref_pseudo = row[('Pseudo Gene', self.pseudo_genes[0], 'Ref')]
            if ref_true != ref_pseudo:
                return pd.Series([ref_true, ref_pseudo])
            sorted = row.loc[list(("Sum", "Sum", base) for base in BASES)].argsort().tolist()[::-1]
            base_i = INDEX_TO_BASE_MAPPING[sorted[0]]
            base_j = INDEX_TO_BASE_MAPPING[sorted[1]]
            if base_i == ref_true:
                ref_pseudo = base_j
            else:
                ref_pseudo = base_i
            return pd.Series([ref_true, ref_pseudo])

        refs = intermediate_summary.apply(_calculate_ref, axis=1)
        refs.columns = ["ref_true", "ref_pseudo"]

        ref = refs["ref_true"].str.cat(refs["ref_pseudo"], sep=':', na_rep='-')
        ref.rename("ref_change", inplace=True)
        return ref, refs

    def calculate_bases(self, sum_df, refs):
        return pd.concat([
            pd.Series([sum_df.loc[index, ("Sum", "Sum", col_name)] for index, col_name in enumerate(refs["ref_true"])], name="base_num_true"),
            pd.Series([sum_df.loc[index, ("Sum", "Sum", col_name)] for index, col_name in enumerate(refs["ref_pseudo"])], name="base_num_pseudo")], 
            axis=1)

    def calculate_simple_ratio(self, bases):
        def _calculate_simple_ratio(row):
            if row["base_num_true"] != 0 and row["base_num_pseudo"] != 0:
                return find_closest_reduced_fraction(row["base_num_true"] / row["base_num_pseudo"], n=self.CNV_ratio[0], m=self.CNV_ratio[1])
            return 0, 0
        def merge_simple_ratio(row):
            if row["simple_true"] != 0 and row["simple_pseudo"] != 0:
                return f'{row["simple_true"]}:{row["simple_pseudo"]}'
            return "NA"

        simple_ratioes = bases.apply(_calculate_simple_ratio, axis=1, result_type='expand')
        simple_ratioes.columns = ["simple_true", "simple_pseudo"]
        simple_ratio = simple_ratioes.apply(merge_simple_ratio, axis=1)
        simple_ratio.rename("reduced_ratio", inplace=True)
        return simple_ratio, simple_ratioes

    def calculate_scaled_ratio(self, bases, scale_factor, num_gene):
        def _calculate_scaled_ratio(row):
            if row["base_num_true"] != 0 and row["base_num_pseudo"] != 0:
                return str(round(row["base_num_true"] / (row["base_num_true"] + row["base_num_pseudo"]) * num_gene * scale_factor * PLOIDY, 2)) + \
                    ':' + str(round(row["base_num_pseudo"] / (row["base_num_true"] + row["base_num_pseudo"]) * num_gene * scale_factor * PLOIDY, 2))
            return "NA"

        scaled_ratio = bases.apply(_calculate_scaled_ratio, axis=1)
        scaled_ratio.rename("scaled_copy_numbers", inplace=True)
        return scaled_ratio
    
    def calculate_delta(self, bases, simple_ratioes):
        def _calculate_delta(row):
            if row["base_num_pseudo"] == 0 or row["simple_pseudo"] == 0:
                return np.nan
            return round(abs(row["base_num_true"] / row["base_num_pseudo"] - row["simple_true"] / row["simple_pseudo"]), 4)
        delta = pd.concat([bases, simple_ratioes], axis=1).apply(_calculate_delta, axis=1)
        delta.rename("delta", inplace=True)
        return delta

    def calculate_sum(self, raw_summary):
        sum_each_base = raw_summary['True Gene'][self.true_gene][BASES] + sum(raw_summary['Pseudo Gene'][pseudo_gene][BASES] for pseudo_gene in self.pseudo_genes)
        sum_df = pd.DataFrame(sum_each_base, columns=BASES)
        sum_df["Sum"] = sum_df.sum(axis=1)
        return pd.concat([sum_df], axis=1, keys=[("Sum", "Sum")])

    def make_summary(self, details, scale_factor, case_name):     
        detail = pd.concat(details, axis=0, ignore_index=True)
        detail = detail.sort_values(by=("True Gene", self.true_gene, "Pos"), axis=0, ascending=True, ignore_index=True)
        raw_summary = self.call_summary(detail)
        sorted_raw_summary = raw_summary.sort_values(by=("True Gene", self.true_gene, "Pos"), axis=0, ascending=True, ignore_index=True)

        # make Sum sub-dataframe
        _sum = self.calculate_sum(sorted_raw_summary)

        intermediate_summary = pd.concat([_sum, sorted_raw_summary], axis=1)

        # make Ratio sub-dataframe
        alt_ref, alt_refs = self.calculate_refs(intermediate_summary)
        bases = self.calculate_bases(_sum, alt_refs)
        simple_ratio, simple_ratioes = self.calculate_simple_ratio(bases)
        scaled_ratio = self.calculate_scaled_ratio(bases, scale_factor, self.num_gene)
        delta = self.calculate_delta(bases, simple_ratioes)

        input_path = self.tmp_dir / self.vcf_path_template.format(case_name=case_name, gene=self.gene)
        output_path = self.tmp_dir / self.nirvana_output_path_template.format(case_name=case_name, gene=self.gene)

        vcf = create_vcf(intermediate_summary, alt_refs, case_name, self.true_gene, self.pseudo_genes[0], input_path)
        run_nirvana(self.nirvana_path, input_path, output_path, self.ref)
        aa_change = get_aa_change_from_nirvana(self.NM_number, vcf, output_path)

        ratio = pd.concat([alt_ref, simple_ratio, scaled_ratio, delta, aa_change], axis=1)
        ratio.columns = pd.MultiIndex.from_product([["Ratio"], ["Ratio"], ratio.columns])
        df = pd.concat([ratio, intermediate_summary], axis=1)
        return df

    def calculate_worker(self, scale_factor, case_name, refs, bam):
        details = self.seq_align_map.apply(self.calculate_df, axis="columns", args=[PseudoBam(bam)])
        details = [pd.concat([ref, detail], axis=1) for ref, detail in zip(refs, details)]

        first_level_col_name = ["Ref", "Chrom", "Pos"] + BASES + ["Sum"]
        top_level_col_name = list(product(["True Gene"], [self.true_gene], first_level_col_name)) + list(product(["Pseudo Gene"], self.pseudo_genes, first_level_col_name))
        col_name = pd.MultiIndex.from_tuples(top_level_col_name)
        details = [detail[col_name] for detail in details]
        summary = self.make_summary(details, scale_factor, case_name)
        summary.set_index(('Ratio', 'Ratio', 'aa_change'), inplace=True)

        summary_path = self.tmp_dir / self.summary_template.format(gene=self.gene, case_name=case_name)
        detail_path = self.tmp_dir / self.details_template.format(gene=self.gene, case_name=case_name)

        with pd.ExcelWriter(summary_path) as writer:
            summary.to_excel(writer, sheet_name=self.gene, freeze_panes=(4, 0))
        if self.debug:
            with pd.ExcelWriter(detail_path) as writer:
                for num, detail in enumerate(details):
                    detail.to_excel(writer, sheet_name=f"{self.gene}_{num}", freeze_panes=(3, 0))
        return summary_path, detail_path

    def calculate(self, n_jobs=1):
        """
        interface to calculate results, delegate actual work to calculate_worker
        """
        results = dict()
        summaries = dict()
        details = dict()

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            for case_name, scale_factor, bam in zip(self.case_names, self.theta_i, self.bam_list):
                result = executor.submit(self.calculate_worker, scale_factor, case_name, self.refs, bam)
                results[result] = case_name

            for result in concurrent.futures.as_completed(results):
                case_name = results[result]
                summaries[case_name], detail = result.result()
                if self.debug:
                    details[case_name] = detail
        return summaries, details, self.theta_i