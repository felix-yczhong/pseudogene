import json
import os
import pathlib

import click

from pseudogene.build import build_db
from pseudogene.query import query_alignment
from pseudogene.run import run

from pseudogene.utils import *

@click.group()
@click.version_option()
@click.pass_context
def main(context):

    # preprocessing config.json
    config_loc = os.path.join(get_project_root(), 'config.json')
    config = json.load(open(config_loc))

    for dict_ in [config["data"]["fasta"], config["data"]["sam"]]:
        for key, value in dict_.items():
            if not os.path.isabs(value):
                value = os.path.join(get_project_root(), value)
            dict_[key] = value

    value = config["data"]["database_path"]
    if not os.path.isabs(value):
        value = os.path.join(get_project_root(), value)
    config["data"]["database_path"] = value

    for gene_group in config['gene_group']:
        for true_gene in config[gene_group]['true_gene']:
            for genome_ver in ['hg38']:
                exception = config[gene_group][genome_ver][true_gene]['exception']
                if len(exception) == 0:
                    exception = None
                elif any(len(e) != 3 for e in exception):
                    raise ValueError("exception point format error")
                else:
                    exception = {tuple(e[:2]): e[2] for e in exception}
                config[gene_group][genome_ver][true_gene]['exception'] = exception
        for pseudo_gene in config[gene_group]['pseudo_gene']:
            for genome_ver in ['hg38']:
                exception = config[gene_group][genome_ver][pseudo_gene]['exception']
                if len(exception) == 0:
                    exception = None
                elif any(len(e) != 3 for e in exception):
                    raise ValueError("exception point format error")
                else:
                    exception = {tuple(e[:2]): e[2] for e in exception}
                config[gene_group][genome_ver][pseudo_gene]['exception'] = exception
    # load config.json and initialize context
    context.ensure_object(dict)
    context.obj['config'] = config

    # initialize logging

main.add_command(build_db)
main.add_command(query_alignment)
main.add_command(run)

if __name__ == "__main__":
    main()