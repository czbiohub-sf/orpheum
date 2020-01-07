from collections import defaultdict
import glob
import logging
import time


import click
import click_log
import pandas as pd
from pyarrow.lib import ArrowIOError
from tqdm import tqdm
import vaex

logger = logging.getLogger(__file__)

@click.command()
@click.argument('parquet')
@click.option('--output-folder', default='.',
              help='Where to output the converted files, by default this '
                   'folder')
@click.option('--output-format', default='hdf5',
              help="Output format for ")
@click_log.simple_verbosity_option(logger)
def cli(parquet, output_folder=None):
    """Make a peptide bloom filter for your peptides

    \b
    Parameters
    ----------
    reads : str
        Sequence file of reads to filter
    peptides : str
        Sequence file of peptides
    peptide_ksize : int
        Number of characters in amino acid words
    long_reads
    verbose

    \b
    Returns
    -------

    """
    # \b above prevents rewrapping of paragraph
    times = {}
    prefix = parquet.split('.parquet')[0]
    # logger.info(prefix)
    try:
        t0 = time.time()
        df = pd.read_parquet(parquet)
        t1 = time.time()
    except ArrowIOError:
        logging.warning(f'Could not read: {parquet}')
    times['pandas__read_parquet'] = t1 - t0

    t6 = time.time()
    df_vx = vaex.from_pandas(df)
    t7 = time.time()
    times['vaex__from_pandas'] = t7 - t6

    arrow = f'{output_folder}/{prefix}.arrow'

    t8 = time.time()
    df_vx.export_arrow(arrow, progress=True)
    t9 = time.time()
    times['vaex__export_arrow'] = t9 - t8

    times = pd.Series(times, name=parquet)
    logger.debug(times.to_csv(index=True, header=True))


if __name__ == "__main__":
    cli()
