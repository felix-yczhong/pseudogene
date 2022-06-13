import tempfile
import pathlib
import os
import subprocess
import concurrent

from collections import defaultdict

import pandas as pd
import click

from pseudogene.calculator import PseudoGeneCalculator
from pseudogene.utils import *

def preprocess_input_files(file_list):
    """
    handle duplicate file paths and duplicate file names
    * for duplicate file paths: remove duplicates
    * for duplicate file names: append _dup{num} to file name
    """
    file_list = {os.path.normpath(file_path) for file_path in file_list}

    caseName_to_filePath_mapping = dict()
    counts = defaultdict(lambda: 0)
    for file_path in file_list:
        case_name = extract_case_name(file_path)
        if (count := counts[case_name]) != 0:
            case_name = case_name + f'_dup{count}'
            counts[case_name] += 1
            # log warning
        caseName_to_filePath_mapping[case_name] = file_path
    return caseName_to_filePath_mapping

def style_func(df, true_gene, pseudo_gene):
    import itertools
    p1 = itertools.product(['Sum'], ['bases'], BASES)
    p2 = itertools.product(['true'], true_gene, ['ref'] + BASES)
    p3 = itertools.product(['pseudo'], pseudo_gene, ['chr', 'pos', 'ref'] + BASES)
    df.style.highlight_max(axis=0).to_excel('styled.xlsx', engine='openpyxl')
    return df

def worker(file_dict, reference, config, gene_group, tmp_dir, ncpus, debug):
    # p = subprocess.Popen(
    #     [path_of_exe, task_id, task_delay],
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE,        
    # )
    return PseudoGeneCalculator(file_dict, reference, config, gene_group=gene_group, tmp_dir=tmp_dir, ncpus=ncpus, debug=debug).calculate()


@click.command()
@click.option('--ncpus', 
              default=os.cpu_count(), 
              type=int, 
              help='number of cores to use.',
              show_default=True)
@click.option('-r', '--reference',
              default="hg38",
              type=click.Choice(["hg19", "hg38"]),
              help='reference genome that was used for alignment',
              show_default=True)
@click.option('--debug',
              is_flag=True,
              default=False,
              type=bool,
              help='write out scale factors and all sequence alignment details instead of just summary. Do NOT run debug mode with large number of genes.')
@click.option('--profile',
              is_flag=True,
              default=False,
              type=bool,
              help='enable profiling')
@click.option('-o', '--output',
              type=click.Path(exists=True, file_okay=False, writable=True, path_type=pathlib.Path),
              default=pathlib.Path(os.getcwd()),
              nargs=1,
              show_default=True,
              help="specify output directory.")
@click.argument("bam_list",
                type=click.Path(exists=True),
                nargs=-1,
                required=True)
@click.pass_context
def run(context, bam_list, output, reference, ncpus, debug, profile):
    """
    Run Pseudogene calculation tool.
    """
    print(bam_list, output, debug, reference, ncpus)

    config = context.obj['config']
    gene_groups = config["gene_group"]
    output_path_template = config["output"]
    summary_tab_template = config["summary_tab"]

    file_dict = preprocess_input_files(bam_list)

    if profile:
        import cProfile, pstats, atexit, io
        prf = cProfile.Profile()
        prf.enable()
        def exit():
            prf.disable()
            print("Profiling completed")
            # ios = io.StringIO()
            with open('profile.txt', 'w+') as f:
                pstats.Stats(prf, stream=f).sort_stats("cumulative").print_stats()
            # print(ios.getvalue())

        atexit.register(exit)

    # process pool to combine & output files
    with tempfile.TemporaryDirectory() as tmp_dir, concurrent.futures.ProcessPoolExecutor(max_workers=min(len(gene_groups), ncpus)) as executor:
        meta_file_dict = dict()
        scale_factors = dict()

        # --------------- multiprocessing ------------------------
        results = {
            executor.submit(worker, file_dict=file_dict, reference=reference, config=config, gene_group=gene_group, tmp_dir=tmp_dir, ncpus=ncpus, debug=debug)
            : gene_group for gene_group in gene_groups
        }

        for result in concurrent.futures.as_completed(results):
            gene_group = results[result]
            # meta_file_dict[gene_group], scale_factors[gene_group] = result.result()

        # ---------------- single process to evaluate performance -------------------
        # gene_groups = ['GBA']
        # for gene_group in gene_groups:
        #     result = PseudoGeneCalculator(file_dict, reference, config, gene_group=gene_group, tmp_dir=tmp_dir, ncpus=ncpus, debug=debug).calculate()
        #     meta_file_dict[gene_group], scale_factors[gene_group] = result

        # write out scale factors
        # with pd.ExcelWriter(output / output_path_template.format(case_name='scale_factors')) as writer:
        #     df = pd.concat([pd.Series(scale_factors[gene_group], index=file_dict.keys()) for gene_group in gene_groups], axis=1, keys=gene_groups)
        #     df.to_excel(writer, sheet_name=summary_tab_template.format(gene='scale_factors'), freeze_panes=(1, 0))        

        # for case_name, file_path in file_dict.items():
        #     # write info (original file path, time stamp) to excel sheet as header
        #     # output and style
        #     with pd.ExcelWriter(output / output_path_template.format(case_name=case_name), engine='openpyxl') as writer:
        #         for gene in gene_groups:
        #             file = pd.read_csv(meta_file_dict[gene][case_name], sep='\t', header=[0, 1, 2])
        #             file.to_excel(writer, sheet_name=summary_tab_template.format(gene=gene), freeze_panes=(3, 0), index=True, merge_cells=True)