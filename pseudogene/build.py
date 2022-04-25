import pathlib
import logging
import json

import click
import gffutils

@click.command(name="build_db")
@click.option('-r', '--reference', 
              nargs=1,
              type=click.Choice(["hg19", "hg38"]), 
              default="hg38",
              show_default=True)
@click.option('--database_path',
              type=click.Path(writable=True, path_type=pathlib.Path),
              nargs=1,
              help="path where the database will be built. Will use the one in config.json if empty.")
@click.argument('database_source',
                type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path),
                nargs=1)
@click.pass_context
def build_db(context, reference, database_source, database_path):
    """
    Build database for sequence alignment
    """
    config = context.obj['config']
    logger = logging.getLogger()

    write_back_flag = True
    if database_source is not None:
        config["database_source"][reference] = database_source
        write_back_flag = True
    if database_path is not None:
        config["database_path"][reference] = database_path
        write_back_flag = True

    gffutils.create_db(config["database_source"][reference], dbfn=config["database_path"][reference], merge_strategy="create_unique")
    logger.info('Building Gencode database success!')
    if write_back_flag:
        logger.info('Overriding database_source and database_path')
        json.dump(config, open('123.json', 'w+'), separators=(',', ': '), indent='\t')