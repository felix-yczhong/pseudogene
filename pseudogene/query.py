import logging
import subprocess
import os
import tempfile
import re

import click
import gffutils

def find_transcript(db, gene_id, transcript_id):
    for transcript in db.children(gene_id, level=1):
        if transcript.id == transcript_id:
            return (transcript.seqid, transcript.start, transcript.end)
    raise KeyError(f'transcript_id {transcript_id} not found')

def find_exons(db, gene_id, transcript_id):
    genomic_ranges = dict()
    try:
        db[gene_id]
    except gffutils.exceptions.FeatureNotFoundError:
        raise gffutils.exceptions.FeatureNotFoundError(f'gene_id {gene_id} not found! Make sure you get the gene id version right')

    for transcript in db.children(gene_id, level=1):
        if transcript.id == transcript_id:
            seqs = db.children(transcript.id, level=1)
            for seq in seqs:
                if seq.id.startswith('exon'):
                    exon_num = int(seq.attributes['exon_number'][0])
                    genomic_ranges[exon_num] = (seq.seqid, seq.start, seq.end)
            if len(genomic_ranges) == 0:
                raise KeyError('no exons found!')
            return genomic_ranges
    raise KeyError(f'transcript_id {transcript_id} not found! Make sure you get the transcript id version right')
    
def construct_search_info(config, gene_group, reference, gene_type_key):
    search_info = dict()
    for gene_name in config[gene_group][gene_type_key]:
        gene_id = config[gene_group][reference][gene_name]['gene_id']
        transcript_id = config[gene_group][reference][gene_name]['transcript_id']
        search_info[gene_name] = (gene_id, transcript_id)
    return search_info

def call_BLAST(true_gene_exon, pseudo_gene_transcript, config, reference, debug, gene_group):
    pattern = re.compile('Subject_(?P<exon_num>[\d]+)')
    with tempfile.TemporaryDirectory() as dir:
        # for true_gene, info in true_gene_exon.items():
        #     true_exon_file_name = config["data"]["fasta"]["true_exon"].format(gene=true_gene, genome_ver=reference)
        #     for pseudo_gene, (pseudo_chr, pseudo_start, pseudo_end) in pseudo_gene_transcript.items():
        #         pseudo_seq_file_name = config["data"]["fasta"]["pseudo_seq"].format(gene=pseudo_gene, genome_ver=reference)
        #         aligned_seq_file_name = config["data"]["sam"]["aligned_seq"].format(true_gene=true_gene, pseudo_gene=pseudo_gene, genome_ver=reference)
        true_exon_file_name = config["data"]["fasta"]["all_true_exon"].format(gene_group=gene_group, genome_ver=reference)
        pseudo_seq_file_name = config["data"]["fasta"]["all_pseudo_seq"].format(gene_group=gene_group, genome_ver=reference)
        aligned_seq_file_name = config["data"]["sam"]["aligned_seq"].format(gene_group=gene_group, genome_ver=reference)
        with tempfile.NamedTemporaryFile(mode='w+', dir=dir) as blast_output:
        # with open('test_.sam', 'w+') as blast_output:
            #  tempfile.NamedTemporaryFile(mode='w+', dir=dir) as blast_output_format_fixed:
            blast_command = config["tools"]["BLAST+_command"].format(true_exon_file_name=true_exon_file_name, pseudo_seq_file_name=pseudo_seq_file_name, blast_output=blast_output.name, )
            os.system(blast_command)

            # process blast output by replacing subject & query to real genomic positions
            # new_lines = list()
            # for line in blast_output:
            #     if line.startswith('@SQ'):
            #         line = line.replace('Query_1', f'{pseudo_chr}:{pseudo_start}-{pseudo_end}')
            #     elif line.startswith('Subject_'):
            #         m = pattern.match(line)
            #         if m:
            #             exon_num = int(m.group('exon_num'))
            #             true_chr, true_start, true_end = info[exon_num]
            #             line = line.replace(f'Subject_{exon_num}', f'{exon_num}:{true_chr}:{true_start}-{true_end}')
            #             line = line.replace('Query_1', f'{pseudo_chr}:{pseudo_start}-{pseudo_end}')
            #     new_lines.append(line)
            # blast_output_format_fixed.writelines(new_lines)
            # blast_output_format_fixed.seek(os.SEEK_SET)

            # calmd
            calmd_command = f"samtools calmd -b {blast_output.name} {true_exon_file_name} > {aligned_seq_file_name}"
            os.system(calmd_command)
            # print(pseudo_seq_file_name)
            if debug:
                debug_calmd_command = f"samtools calmd {blast_output.name} {true_exon_file_name} > {aligned_seq_file_name[:-4] + '.sam'}"
                os.system(debug_calmd_command)
            # index
            index_command = f"samtools index {aligned_seq_file_name}"
            os.system(index_command)

@click.command(name="query_alignment")
@click.option('--debug',
              is_flag=True,
              default=False,
              type=bool,
              help="output aligned sequence in sam format")
@click.option('-r', '--reference', 
              nargs=1,
              type=click.Choice(["hg19", "hg38"]),
              default="hg38",
              show_default=True)
# @click.option('--pseudo_gene', 
#               nargs=2,
#               type=click.STRING,
#               multiple=True,
#               help="takes 2 arguments, (pseudo gene name, pseudo gene id). accept multiple times for multiple pseudo genes")
# @click.option('--true_gene', 
#               nargs=3,
#               type=click.STRING,
#             #   multiple=True,
#               help="takes 3 arguments, (true gene name, true gene id, true gene transcript id)")
@click.option('--gene_groups', 
              nargs=1,
              type=click.STRING,
              help="query only for specific gene group(s), gene group name(s) must be in config.json")
@click.pass_context
def query_alignment(context, reference, gene_groups, debug):
    """
    Query sequence alignment.
    """
    config = context.obj['config']
    db = gffutils.FeatureDB(config["data"]["database_path"][reference].format(genome_ver=reference))
    logger = logging.getLogger()

    # print(reference, gene_group, debug)
    if gene_groups is None:
        gene_groups = config['gene_group']
    else:
        gene_groups = [gene_groups]

    for gene_group in gene_groups:
        # construct search info
        true_gene_search_info = construct_search_info(config, gene_group, reference, 'true_gene')
        if len(true_gene_search_info) > 1:
            raise NotImplementedError('gene group with more than 1 true genes is not implemented.')
        pseudo_gene_search_info = construct_search_info(config, gene_group, reference, 'pseudo_gene')

        # find true gene & pseudo gene info from db
        true_gene_exon = {true_gene_name: find_exons(db, true_gene_id, true_gene_transcript_id) for true_gene_name, (true_gene_id, true_gene_transcript_id) in true_gene_search_info.items()}
        true_gene_transcript = {true_gene_name: find_transcript(db, true_gene_id, true_transcript_id) for true_gene_name, (true_gene_id, true_transcript_id) in true_gene_search_info.items()}  
        pseudo_gene_transcript = {pseudo_gene_name: find_transcript(db, pseudo_gene_id, pseudo_transcript_id) for pseudo_gene_name, (pseudo_gene_id, pseudo_transcript_id) in pseudo_gene_search_info.items()}

        # find seq and save to disk
        # save true gene exons
        true_file_path = config["data"]["fasta"]["all_true_exon"].format(gene_group=gene_group, genome_ver=reference)
        with open(true_file_path, 'w+') as true_file:
            for true_gene, info in true_gene_exon.items():
                output_path = config["data"]["fasta"]["true_exon"].format(gene=true_gene, genome_ver=reference)
                with open(output_path, 'w+') as file:
                    for exon_num in sorted(info.keys()):
                        chr, start, end = info[exon_num]
                        true_ref = subprocess.run(['samtools', 'faidx', config["tools"]["fasta_loc"][reference], f'{chr}:{start}-{end}'], stdout=subprocess.PIPE).stdout.decode()
                        header = f'>{true_gene}-{exon_num}-{chr}-{start}-{end}\n'
                        true_ref = true_ref.split('\n', maxsplit=1)[1]
                        file.write(header + true_ref)
                        true_file.write(header + true_ref)
                subprocess.run(['samtools', 'faidx', output_path])
            subprocess.run(['samtools', 'faidx', true_file_path])

        # save entire true gene
        for true_gene, (chr, start, end) in true_gene_transcript.items():
            true_ref = subprocess.run(['samtools', 'faidx', config["tools"]["fasta_loc"][reference], f'{chr}:{start}-{end}'], stdout=subprocess.PIPE).stdout.decode()
            output_path = config["data"]["fasta"]["true_seq"].format(gene=true_gene, genome_ver=reference)
            with open(output_path, 'w+') as f:
                f.write(true_ref)
        subprocess.run(['samtools', 'faidx', output_path])

        # save entire pseudo gene
        pseudo_file_path = config["data"]["fasta"]["all_pseudo_seq"].format(gene_group=gene_group, genome_ver=reference)
        with open(pseudo_file_path, 'w+') as pseudo_file:
            for pseudo_gene, (chr, start, end) in pseudo_gene_transcript.items():
                pseudo_ref = subprocess.run(['samtools', 'faidx', config["tools"]["fasta_loc"][reference], f'{chr}:{start}-{end}'], stdout=subprocess.PIPE).stdout.decode()
                output_path = config["data"]["fasta"]["pseudo_seq"].format(gene=pseudo_gene, genome_ver=reference)
                with open(output_path, 'w+') as f:
                    f.write(pseudo_ref)
                pseudo_ref = pseudo_ref.split('\n', maxsplit=1)[1]
                header = f'>{pseudo_gene}--{chr}-{start}-{end}\n'
                pseudo_file.write(header + pseudo_ref)
                subprocess.run(['samtools', 'faidx', output_path])
            subprocess.run(['samtools', 'faidx', pseudo_file_path])

        call_BLAST(true_gene_exon, pseudo_gene_transcript, config, reference, debug, gene_group)