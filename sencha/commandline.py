# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain

# Python standard library imports
from functools import partial

# 3rd party libraries
import click

# Within-module imports
from sencha.compare_kmer_content import cli as compare_kmers
from sencha.translate import cli as translate
from sencha.index import cli as index

click.option = partial(click.option, show_default=True)

settings = dict(help_option_names=["-h", "--help"])


@click.group(
    options_metavar="", subcommand_metavar="<command>", context_settings=settings
)
def cli():
    """
    Kmer hashing tools contains data cleaning and visualization code for
    analyzing sequencing datasets at the k-mer level
    Kmer hashing tools contains data cleaning and visualization code for
    analyzing kmer-hashing similarity matrices
    """
    pass


cli.add_command(compare_kmers, name="compare-kmers")
cli.add_command(index, name="index")
cli.add_command(translate, name="translate")

if __name__ == "__main__":
    cli()
