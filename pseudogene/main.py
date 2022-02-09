import json
import tempfile
import pathlib
import os
import concurrent
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
import click

from pseudogene.calculator import PseudoGeneCalculator
from pseudogene.utils import *

def worker_func(gene, bam_list, case_names, reference, config, ncpus, tmp_dir, debug):
    calculator = PseudoGeneCalculator(bam_list=bam_list, case_names=case_names, ref=reference, 
                                        config=config, gene=gene, n_jobs=ncpus, tmp_dir=tmp_dir, debug=debug)
    return calculator.calculate(ncpus)

@click.command()
@click.option('--ncpus', 
              default=1, 
              type=int, 
              help='number of cores to use',
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
@click.option('-o', '--output',
              type=click.Path(exists=True, file_okay=False, writable=True, path_type=pathlib.Path),
              default=pathlib.Path(os.getcwd()),
              nargs=1,
              help="specify output directory.")
@click.argument("bam_list",
                type=click.Path(exists=True),
                nargs=-1,
                required=True)
def main(bam_list, ncpus, reference, output, debug):
    """
    Pseudogene calculation tool.
    """
    if not bam_list:
        ctx = click.get_current_context()
        ctx.get_help()
        ctx.exit()

    config = json.load(open(get_project_root() / 'config.json'))
    genes = config["genes"]
    output_path_template = config["output"]
    summary_tab_template = config["summary_tab"]
    detail_tab_template = config["detail_tab"]

    bam_list = sorted(bam_list, key=lambda bam: extract_case_name(bam))
    case_names = [extract_case_name(bam) for bam in bam_list]
    
    results = dict()
    summaries = dict()
    details = dict()
    scale_factors = dict()
    with tempfile.TemporaryDirectory() as tmp_dir, ProcessPoolExecutor(max_workers=min(len(config["genes"]), ncpus)) as executor:
        tmp_dir_path = pathlib.Path(tmp_dir)
        for gene in genes:
            result = executor.submit(worker_func, bam_list=bam_list, case_names=case_names, reference=reference, 
                                    config=config, gene=gene, ncpus=ncpus, tmp_dir=tmp_dir_path, debug=debug)
            results[result] = gene

        for result in concurrent.futures.as_completed(results):
            gene = results[result]
            summaries[gene], details[gene], scale_factors[gene] = result.result()
    
        # combine summaries & details
        for case_name in case_names:
            with pd.ExcelWriter(output / output_path_template.format(case_name=case_name)) as writer: 
                for gene in genes:
                    dfs = pd.read_excel(summaries[gene][case_name], sheet_name=None, header=[0, 1, 2], index_col=0, engine="openpyxl")
                    dfs[gene].to_excel(writer, sheet_name=summary_tab_template.format(gene=gene), freeze_panes=(4, 0))

                    if debug:
                        dfs = pd.read_excel(details[gene][case_name], sheet_name=None, header=[0, 1, 2], index_col=0, engine="openpyxl")
                        for num, df in dfs.items():
                            df.to_excel(writer, sheet_name=detail_tab_template.format(gene=gene, num=num), freeze_panes=(3, 0))
    
    if debug:                    
        with pd.ExcelWriter(output / "scale_factor.xlsx") as writer:
            df = pd.concat([pd.Series(scale_factors[gene], index=case_names) for gene in genes], axis=1, keys=genes)
            df.to_excel(writer, sheet_name="scale_factors", freeze_panes=(1, 0))

if __name__ == "__main__":
    main()