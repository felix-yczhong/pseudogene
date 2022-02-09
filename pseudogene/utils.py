import os
import pathlib
import json

import numpy as np
import pandas as pd

BASES = ['A', 'C', 'G', 'T', 'N']
BASE_MATCH = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}

INDEX_TO_BASE_MAPPING = dict(enumerate(BASES))
BASE_TO_INDEX_MAPPING = {value: key for key, value in INDEX_TO_BASE_MAPPING.items()}
PLOIDY = 2 # assume control genes all diploid

def extract_case_name(input):
    input = input.rsplit(sep='/', maxsplit=1)[1]
    return input.split(sep='.', maxsplit=1)[0]

def get_project_root() -> pathlib.Path:
    return pathlib.Path(__file__).parent.parent

def find_closest_simple_fraction(R, n=5, m=5):
    """
    find with Farey Sequence simple fraction i/j closest to a real number R <= 1
    :param n: candidate of numerator from 1 to n, where 1 <= i <=n
    :param m: candidate of denominator from 1 to m, where 1 <= j <= m
    """
    a_num = 0
    a_denom = b_num = b_denom = 1

    if isinstance(n, int):
        it1 = range(1, n + 1)
    else:
        it1 = n
    if isinstance(m, int):
        it2 = range(1, m + 1)
    else:
        it2 = m

    for c_num, c_denom in zip(it1, it2):
        c_num = a_num + b_num
        c_denom = a_denom + b_denom

        if c_num > n or c_denom > m:
            if R - a_num/a_denom < b_num/b_denom - R:
                return a_num, a_denom
            else:
                return b_num, b_denom

        if c_num/c_denom < R:
            a_num = c_num
            a_denom = c_denom
        else:
            b_num = c_num
            b_denom = c_denom
    return 0, 0

def convert_coordinate_one_to_zero(genomic_range):
    """
    convert 1-based coordinate to 0-based coordinate
    :param start: 1-based start inclusive
    :param stop: 1-based stop exclusive
    :return start, stop: start and stop position in 0-based coordinate
    """
    c, start, stop = genomic_range
    return [c, start - 1, stop - 1]

def convert_coordinate_zero_to_one(genomic_range):
    """
    convert 0-based coordinate to 1-based coordinate
    :param start: 0-based start inclusive
    :param stop: 0-based stop exclusive
    :return start, stop: start and stop position in 1-based coordinate
    """
    c, start, stop = genomic_range
    return [c, start + 1, stop + 1]

def convert_coordinate_one_to_zero_df(row):
    genomic_range = [row["Chrome"], row["Start"], row["End"]]
    new_row = row.copy()
    c, new_row["Start"], new_row["End"] = convert_coordinate_one_to_zero(genomic_range)
    return new_row

def reconstruct(df, convert_func=convert_coordinate_one_to_zero_df):
    new_sub_dfs = []
    for sub_df in [df["True Gene"], df["Pseudo Gene"]]:
        names = []
        subsub_dfs = []
        for name, subsub_df in sub_df.groupby(axis=1, level=0):
            new_subsub_df = subsub_df.apply(convert_func, axis=1)
            subsub_dfs.append(new_subsub_df)
            names.append(name)
        new_sub_df = pd.concat(subsub_dfs, axis=1, names=names)
        new_sub_dfs.append(new_sub_df)
    new_df = pd.concat(new_sub_dfs, axis=1, names=["True Gene", "Pseudo Gene"])
    return new_df

def run_nirvana(tool_path, input_path, output_path, ref):
    ref_map = {"hg19": "GRCh37", "hg38": "GRCh38"}
    genomver = ref_map[ref]
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

def create_vcf(summary, alt_refs, case_name, true_gene, pseudo_gene, vcf_path):
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', f'{case_name}']
    df = pd.DataFrame()
    df[['#CHROM', 'POS']] = summary[[('True Gene', true_gene, 'Chrom'), ('True Gene', true_gene, 'Pos')]]
    df[['REF', 'ALT']] = alt_refs  
    for col_name in ['ID', 'QUAL', 'FILTER', 'INFO' , 'FORMAT', f'{case_name}']:
        df[col_name] = pd.Series(['.'] * len(summary))

    df = df[columns]

    header = "##fileformat=VCFv4.3\n"
    with open(vcf_path, 'w+') as f:
        f.write(header)
    df.to_csv(vcf_path, sep='\t', index=False, mode='a')
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
    aa_change.rename("aa_change", inplace=True)
    return aa_change
