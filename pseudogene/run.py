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

def worker(test_file_dict, control_file_dict, reference, config, gene_group, output_path, tmp_dir, ncpus, debug):
    # p = subprocess.Popen(
    #     [path_of_exe, task_id, task_delay],
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE,        
    # )
    return PseudoGeneCalculator(test_file_dict, control_file_dict, reference, config, gene_group, output_path, tmp_dir=tmp_dir, ncpus=ncpus, debug=debug).calculate()

def output_summary(output_path, sheet_name, gene_groups, meta_file_dict):
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        for gene_group in gene_groups:
            file = pd.read_csv(meta_file_dict[gene_group], sep='\t', header=[0, 1, 2])
            # TODO: style df before output
            file.to_excel(writer, sheet_name=sheet_name.format(gene=gene_group), freeze_panes=(3, 0), index=True, merge_cells=True)

@click.command()
@click.option('--ncpus', 
              default=os.cpu_count(), 
              type=int, 
              help='number of cores to use.',
              show_default=True)
@click.option('-r', '--reference',
              default="hg38",
              type=click.Choice(["hg19", "hg38"]),
              help='reference genome used for alignment',
              show_default=True)
@click.option('--debug',
              is_flag=True,
              default=False,
              type=bool,
              help='write out details such as scale factors, control genes average coverage and true gene & pseudogene average coverage.')
@click.option('--profile_path',
              type=click.Path(exists=True, dir_okay=False, writable=True, path_type=pathlib.Path),
              default=pathlib.Path(os.getcwd()) / 'profile.txt',
              show_default=True,
              nargs=1,
              help='specify output path of profiling result, not used if profiling is off')
@click.option('--profile',
              is_flag=True,
              default=False,
              type=bool,
              help='enable profiling.')
@click.option('-o', '--output',
              type=click.Path(exists=True, file_okay=False, writable=True, path_type=pathlib.Path),
              default=pathlib.Path(os.getcwd()),
              nargs=1,
              show_default=True,
              help="specify output directory.")
@click.option('-c', '--control',
              type=click.Path(exists=True, dir_okay=False, readable=True, path_type=pathlib.Path),
              required=False,
              help="takes a file path which contains all control samples' file paths")
@click.argument("bam_list",
                type=click.Path(exists=True),
                nargs=-1,
                required=True)
@click.pass_context
def run(context, bam_list, control, output, profile, profile_path, reference, ncpus, debug):
    """
    Run Pseudogene calculation tool.
    """
    print(bam_list, output, debug, reference, ncpus, profile)

    if profile:
        import cProfile, pstats, atexit, io
        prf = cProfile.Profile()
        prf.enable()
        def exit():
            prf.disable()
            print("Profiling completed")
            # ios = io.StringIO()
            with open(profile_path, 'w+') as f:
                pstats.Stats(prf, stream=f).sort_stats("cumulative").print_stats()
            # print(ios.getvalue())

        atexit.register(exit)

    config = context.obj['config']
    gene_groups = config["gene_group"]

    test_file_dict = preprocess_input_files(bam_list)

    if control is not None:
        with open(control) as c:
            control_list = [pathlib.Path(line.strip()) for line in c]
        for file_path in control_list:
            if not file_path.exists():
                raise FileNotFoundError(f'{file_path} not found!')
            if not file_path.is_file():
                raise Exception(f'{file_path} is not a file!')
        control_file_dict = preprocess_input_files(control_list)
    else:
        control_file_dict = dict()
    
    # process pool to combine & output files
    with tempfile.TemporaryDirectory() as tmp_dir:
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(len(gene_groups), ncpus)) as executor:
            meta_file_dict = dict()

            # --------------- serialized for debugging
            # for gene_group in gene_groups:
            #     cal = PseudoGeneCalculator(test_file_dict, control_file_dict, reference, config, gene_group, output, ncpus=ncpus, tmp_dir=tmp_dir, debug=debug)
            #     meta_file_dict[gene_group] = cal.calculate()

            # --------------- multiprocessing
            results = {
                executor.submit(worker, test_file_dict=test_file_dict, control_file_dict=control_file_dict, reference=reference, config=config, gene_group=gene_group, output_path=output, tmp_dir=tmp_dir, ncpus=ncpus, debug=debug)
                : gene_group for gene_group in gene_groups
            }

            for result in concurrent.futures.as_completed(results):
                gene_group = results[result]
                meta_file_dict[gene_group] = result.result()

        with concurrent.futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
            summaries = {
                executor.submit(output_summary, output_path=output / config["output"].format(case_name=case_name), \
                                sheet_name=config["summary_tab"], \
                                gene_groups=gene_groups, meta_file_dict={gene_group: meta_file_dict[gene_group][case_name] for gene_group in gene_groups}) \
                : case_name for case_name in test_file_dict.keys()
            }

            for summary in concurrent.futures.as_completed(summaries):
                case_name = summaries[summary]
                summary.result()
                #TODO log completion or incompletion messages