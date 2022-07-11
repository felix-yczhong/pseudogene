
from collections import defaultdict
import tempfile

import numpy as np
import pandas as pd

import smaca.constants as C
from pseudogene.bam import PseudoBam
from pseudogene.utils import *

def return_list():
    return list()

def read_input_CopyOnWrite(file_path, ref, \
                            true_gene_region, pseudo_gene_region, \
                            true_gene_pos, pseudo_gene_pos, \
                            true_genes, pseudo_genes,
                            true_gene_pos_len, pseudo_gene_pos_len):
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
    true_table = np.zeros((true_gene_pos_len, 4), dtype=np.int64)
    pseudo_table = np.zeros((pseudo_gene_pos_len, 4), dtype=np.int64)

    control_gene_cov = bam_file.get_cov_ranges(C.POSITIONS[ref]["GENES"])
    gene_region_cov = bam_file.get_cov_ranges_(true_gene_region, pseudo_gene_region)

    offset = 0
    for true_gene, pos in true_gene_pos.items():
        for j, pos_info in enumerate(pos):
            true_table[offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(pos_info[:2])).squeeze()
        offset += len(pos)

    offset = 0
    for pseudo_gene, pos in pseudo_gene_pos.items():
        for j, pos_info in enumerate(pos):
            pseudo_table[offset + j, :] = np.asarray(bam_file.get_alignment_cov_util(pos_info[:2])).squeeze()
        offset += len(pos)
    
    return control_gene_cov, gene_region_cov, true_table, pseudo_table

def build_raw_table(idx, alignment_table, true_gene_pos, pseudo_gene_pos, true_table, pseudo_table, true_genes, pseudo_genes):
    pos_col = ['exon', 'chr', 'pos', 'ref', 'alt', 'aa_change']
    ref_sub_pos_col = ['chr', 'pos', 'ref']
    alt_sub_pos_col = ['chr', 'pos', 'alt']

    offset = 0
    for true_gene, pos in true_gene_pos.items():
        join_index = list(pd.MultiIndex.from_product([['True'], [true_gene], ref_sub_pos_col]))
        pos_info = pd.DataFrame(pos, columns=ref_sub_pos_col)
        cov_info = pd.DataFrame(true_table[idx, offset: offset + len(pos), :], columns=BASES)
        to_merged = pd.concat([pos_info, cov_info], axis=1, keys=[('True', true_gene), ('True', true_gene)])
        alignment_table = alignment_table.merge(to_merged, how='left', left_on=join_index, right_on=join_index)
        offset += len(pos)

    offset = 0
    for pseudo_gene, pos in pseudo_gene_pos.items():
        join_index = list(pd.MultiIndex.from_product([['Pseudo'], [pseudo_gene], alt_sub_pos_col]))
        pos_info = pd.DataFrame(pos, columns=alt_sub_pos_col)
        cov_info = pd.DataFrame(pseudo_table[idx, offset: offset + len(pos), :], columns=BASES)
        to_merged = pd.concat([pos_info, cov_info], axis=1, keys=[('Pseudo', pseudo_gene), ('Pseudo', pseudo_gene)])
        alignment_table = alignment_table.merge(to_merged, how='left', left_on=join_index, right_on=join_index)
        offset += len(pos)

    # re-order table columns
    ridx = [('ratio', 'exception', 'is_exception')] + \
            list(pd.MultiIndex.from_product([['True'], true_genes, pos_col[:-2] + BASES])) + \
            list(pd.MultiIndex.from_product([['Pseudo'], pseudo_genes, pos_col[1:] + BASES]))
    alignment_table = alignment_table[ridx]

    alignment_table.sort_values(by=[('True', true_genes[0], 'pos')], inplace=True, ignore_index=True)

    return alignment_table

def calculate_sum(raw_table, true_genes, pseudo_genes):
    true_idx = [pd.MultiIndex.from_product([['True'], [true_gene], BASES]) for true_gene in true_genes]
    pseudo_idx = [pd.MultiIndex.from_product([['Pseudo'], [pseudo_gene], BASES]) for pseudo_gene in pseudo_genes]

    sum_each_base = sum(raw_table[col].to_numpy() for col in true_idx + pseudo_idx)
    sum_df = pd.DataFrame(sum_each_base, columns=pd.MultiIndex.from_product([["Sum"], ["bases"], BASES]))
    return sum_df

def get_refs(row, true_genes, pseudo_genes):
    refs = defaultdict(lambda: list())
    for true_gene in true_genes:
        ref = row[('True', true_gene, 'ref')]
        refs[ref].append(true_gene)

    for pseudo_gene in pseudo_genes:
        ref = row[('Pseudo', pseudo_gene, 'ref')]
        refs[ref].append(pseudo_gene)

    if row[('ratio', 'exception', 'is_exception')]:
        for pseudo_gene in pseudo_genes:
            alt = row[('Pseudo', pseudo_gene, 'alt')]
            refs[alt] # make sure base in the key of refs, but not adding new elements to its value

    for key, value in refs.items():
        if len(value) == 0:
            refs[key].append('other')

    return pd.Series({' & '.join(gene_names): ref for ref, gene_names in refs.items()})

def calculate_refs(raw_table, true_genes, pseudo_genes):
    refs = raw_table.apply(get_refs, axis=1, args=(true_genes, pseudo_genes))
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

def get_ratio(row, scale_factor, num_genes):
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

def calculate_sub_ratio(bases, scale_factor, num_genes):

    ratio = bases.apply(get_ratio, axis=1, args=(scale_factor, num_genes))
    ratio.columns = pd.MultiIndex.from_product([['ratio'], ['Copy Number'], [base[-1] for base in bases.columns]])
    return ratio

def calculate_ratio(raw_table, sum_df, scale_factor, num_gene, true_genes, pseudo_genes):
    """
    1. get the refs per alignment point (from raw_table)
    2. get the bases sum per alignment point (from sum_df)
    3. calculate ratio (infer ratio)
    """
    refs = calculate_refs(raw_table, true_genes, pseudo_genes)
    bases = calculate_bases(refs, sum_df)
    ratios = calculate_sub_ratio(bases, scale_factor, num_gene)

    ratio_table = pd.concat([refs, bases, ratios], axis=1)
    return ratio_table

def calculation_worker_CopyOnWrite(idx, alignment_table, true_gene_pos, pseudo_gene_pos, true_table, pseudo_table, true_genes, pseudo_genes, scale_factor, num_gene, tmp_dir):
    """ 
    format to pandas table (raw_summary, i.e. without sum & ratio) 
    
    """
    raw_table = build_raw_table(idx, alignment_table, true_gene_pos, pseudo_gene_pos, true_table, pseudo_table, true_genes, pseudo_genes)
    sum_df = calculate_sum(raw_table, true_genes, pseudo_genes)
    ratio_df = calculate_ratio(raw_table, sum_df, scale_factor[idx], num_gene, true_genes, pseudo_genes)

    raw_table_without_exception = raw_table.drop(columns=('ratio', 'exception', 'is_exception'))
    is_exception = raw_table[('ratio', 'exception', 'is_exception')]
    summary = pd.concat([is_exception, ratio_df, sum_df, raw_table_without_exception], axis=1)

    f = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)
    summary.to_csv(f, sep='\t', index=False)
    return f.name
