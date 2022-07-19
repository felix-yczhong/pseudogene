from ctypes import alignment
import os
import pathlib
import json
from fractions import Fraction
import functools

import numpy as np
import pandas as pd

BASES = ['A', 'C', 'G', 'T']
BASE_MATCH = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}

INDEX_TO_BASE_MAPPING = dict(enumerate(BASES))
BASE_TO_INDEX_MAPPING = {value: key for key, value in INDEX_TO_BASE_MAPPING.items()}
PLOIDY = 2 # assume control genes all diploid

VCF_HEADER = "##fileformat=VCFv4.3\n"

def extract_case_name(input):
    input = input.rsplit(sep='/', maxsplit=1)[1]
    return input.split(sep='.', maxsplit=1)[0]

def analyze_seq_name(name):
    gene_name, exon_num, chr, start, end = name.split('-')
    start, end = int(start), int(end)
    if exon_num != '':
        exon_num = int(exon_num)
    return gene_name, exon_num, chr, start, end

def get_project_root() -> pathlib.Path:
    return pathlib.Path(__file__).parent.parent

def find_closest_reduced_fraction_F(R, n=5, m=5):
    """
    find with Farey Sequence the closest reduced fraction i/j closest to a real number R <= 1
    :param n: candidate of numerator from 0 to n, where 1 <= i <=n
    :param m: candidate of denominator from 1 to m, where 1 <= j <= m
    """
    a_num = 0
    a_denom = b_num = b_denom = 1
    c_num = a_num + b_num
    c_denom = a_denom + b_denom
    
    while c_num <= n and c_denom <= m:
        if c_num / c_denom < R:
            a_num = c_num
            a_denom = c_denom
        else:
            b_num = c_num
            b_denom = c_denom
        c_num = a_num + b_num
        c_denom = a_denom + b_denom
        
    return (a_num, a_denom) if R - a_num / a_denom <= b_num / b_denom - R else (b_num, b_denom)

def find_closest_reduced_fraction_B(R, n, m):
    """
    find by binary search the closest reduced fraction i/j closest to a real number R <= 1
    :param n: candidate of numerator from 0 to n, where 1 <= i <=n
    :param m: candidate of denominator from 1 to m, where 1 <= j <= m
    """
    candidates = list(set(Fraction(i, j) for i in n for j in m))
    candidates.sort()

    lower = start = 0
    upper = end = len(candidates)

    # lower inclusive, upper exclusive

    while lower < upper:
        mid = (lower + upper) // 2
        value = float(candidates[mid])
        if value < R:
            lower = mid + 1
        elif value > R:
            upper = mid
        else:
            return candidates[mid].as_integer_ratio()
        
    # now lower == upper
    if lower == start:
        return candidates[start].as_integer_ratio()
    elif lower == end:
        return candidates[end - 1].as_integer_ratio()
    else:
        return candidates[lower - 1].as_integer_ratio() if R - candidates[lower - 1] <= candidates[lower] - R else candidates[lower].as_integer_ratio()

def find_closest_reduced_fraction(R, n, m):
    if isinstance(n, int) and n < 1:
        raise ValueError('fraction numerator limit (true gene copy number limit) must be higher than 0')
    if isinstance(m, int) and m < 1:
        raise ValueError('fraction denominator limit (pseudog gene copy number limit) must be higher than 0')
    if isinstance(n, int) and isinstance(m, int) and R < 1:
        return find_closest_reduced_fraction_F(R, n, m)
    
    if isinstance(n, int):
        n = range(n + 1)
    if isinstance(m, int):
        m = range(1, m + 1)
    return find_closest_reduced_fraction_B(R, n, m)

def convert_coordiante_fasta_to_sam(start, end):
    """
    fasta coordinate: 1-based, both inclusive
    sam coordinate: 0-based, start-inclusive, end-exclusive
    """
    return start - 1, end

def run_nirvana(tool_path, input_path, output_path, ref):
    ref_map = {"hg19": "GRCh37", "hg38": "GRCh38"}
    genomver = ref_map[ref]
    print(output_path)
    command = f"dotnet {tool_path}/bin/Release/netcoreapp2.1/Nirvana.dll \
                -c {tool_path}/Data/Cache/{genomver}/Both \
                --sd {tool_path}/Data/SupplementaryAnnotation/{genomver} \
                -r {tool_path}/Data/References/Homo_sapiens.{genomver}.Nirvana.dat \
                -i {input_path} \
                -o {output_path}"
    os.system(command)
    os.system(f"gunzip -f {output_path}.json.gz")
    #TODO return nirvana target file name/path

def run_vep(case_name):
    vcf_path = f'/home/user/SMAca/smaca/{case_name}_called.vcf'
    output_path = f'/home/user/SMAca/smaca/{case_name}_called_vep.tsv'

    tool_path = pathlib.Path("/raid/ensembl-vep-release-101.0")
    genomver = "GRCh38"
    fa = f"{tool_path}/.vep/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
    # fa = f"{base_path}/.vep/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz" for hg19

    cpus = 8
    command = f"perl {tool_path}/vep -i {vcf_path} -o {output_path} \
                --sift b --polyphen b --ccds --uniprot --symbol --numbers --domains --regulatory \
                --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype \
                --af --af_1kg --af_esp --af_gnomad --max_af --pubmed --var_synonyms --variant_class \
                --assembly {genomver} --force_overwrite --fork {cpus} --dir {tool_path}/.vep \
                --offline --fasta {fa} --merged --buffer_size 131072 --plugin MaxEntScan,{tool_path}/fordownload --tab"
    os.system(command)
    #TODO return vep target file name/path

def run_annovar(case_name, panel):
    vcf_path = f'/home/user/SMAca/smaca/{case_name}_called.vcf'
    output_path = f'/home/user/SMAca/smaca/{case_name}_called_annovar.tsv'

    tool_path = pathlib.Path("/raid/annovar-2020Jun08")
    genomver = "hg38"

    cpus = 8
    all_socket = "numactl --interleave=all" if panel == "WGS" else ""
    command = f"{all_socket} perl {tool_path}/table_annovar.pl {vcf_path} {tool_path}/humandb/ --outfile {output_path} --buildver {genomver} \
    --protocol refGene,1000g2015aug_all,esp6500siv2_all,exac03,nci60,avsnp147,cosmic70,clinvar_20210123,dbnsfp35a,gnomad_exome,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene \
    --operation g,f,f,f,f,f,f,f,f,f,f,f,r,g,g \
    --vcfinput --otherinfo --thread {cpus} --maxgenethread {cpus}"
    os.system(command)
    #TODO return annovar target file name/path

def create_vcf(alignment_relation, gene_group):
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', f'{gene_group}']
    sub_columns = columns[5:]

    rows = list()
    for (ref_chr, ref_pos, read_alt), values in alignment_relation.items():
        for (read_chr, read_pos), (ref_gene_name, ref_exon_num, ref_base, read_gene_name, read_base, is_exception) in values.items():
            if ref_base != read_alt:
                rows.append([ref_chr, ref_pos, '.', ref_base, read_alt] + ['.'] * len(sub_columns))
    df = pd.DataFrame(rows, columns=columns)
    df.drop_duplicates(inplace=True)
    df.sort_values(by=columns[:5], inplace=True, ignore_index=True)
    return df

def get_aa_change_from_nirvana(NM_num, vcf, output_path):
    output_path = output_path.parent / (output_path.name + '.json')
    file = json.load(open(output_path))

    def find_func(row, NM_num):
        for entry in file["positions"]:
            if entry["chromosome"] == row["#CHROM"] and \
               entry["position"] == row["POS"] and \
               entry["refAllele"] == row["REF"] and \
               row["ALT"] in entry["altAlleles"]:
                search = entry["variants"][0]["transcripts"]
                if NM_num != None:
                    for t in search:
                        try:
                            if t["source"] == "RefSeq" and t["transcript"] == NM_num:
                                return t["hgvsp"]
                        except:
                            pass
                    for t in search:
                        try:
                            if t["source"] == "RefSeq" and t["transcript"].startswith(NM_num.split('.')[0]):
                                return t["hgvsp"]
                        except:
                            pass
                    for t in search:
                        try:
                            if t["source"] == "RefSeq" and t["isCanonical"]:
                                return t["hgvsp"]
                        except:
                            pass
                    for t in search:
                        try:
                            if t["source"] == "Ensembl" and t["transcript"] == NM_num:
                                return t["hgvsp"]
                        except:
                            pass
                    for t in search:
                        try:
                            if t["source"] == "Ensembl" and t["transcript"].startswith(NM_num.split('.')[0]):
                                return t["hgvsp"]
                        except:
                            pass                    
                    for t in search:
                        try:
                            if t["source"] == "Ensembl" and t["isCanonical"]:
                                return t["hgvsp"]
                        except:
                            pass
                    if len(search) == 1:
                        return search[0]["hgvsp"]
        return "Not Found"

    aa_change = vcf.apply(find_func, axis=1, args=[NM_num])
    aa_change = pd.DataFrame(aa_change, columns=['aa_change'])
    aa_change.reset_index(drop=True, inplace=True)
    return aa_change
