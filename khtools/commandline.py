# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from .os_utils import os  # noqa: F401

# Python standard library imports
from functools import partial

# 3rd party libraries
import click

# Within-module imports
from khtools.extract_coding import cli as extract_coding
from khtools.bloom_filter import cli as bloom_filter
from khtools.pandas2vaex import cli as pandas2vaex

click.option = partial(click.option, show_default=True)

settings = dict(help_option_names=['-h', '--help'])


@click.group(options_metavar='',
             subcommand_metavar='<command>',
             context_settings=settings)
def cli():
    """
    Kmer hashing tools contains data cleaning and visualization code for
    analyzing sequencing datasets at the k-mer level
    Kmer hashing tools contains data cleaning and visualization code for
    analyzing kmer-hashing similarity matrices
    """
    pass


cli.add_command(extract_coding, name='extract-coding')
cli.add_command(bloom_filter, name='bloom-filter')
cli.add_command(pandas2vaex, name='pandas2vaex')

if __name__ == "__main__":
    cli()
