"""
test_homology.py

Tests comparing homologous genes on k-mer content
"""
from click.testing import CliRunner
import os
from io import StringIO

import pandas as pd
import pandas.util.testing as pdt
import pytest


@pytest.fixture
def homologues_tsv(data_folder):
    return os.path.join(data_folder, 'homology',
                        'human-chicken-orthologs-subset.tsv')

@pytest.fixture
def homologues_tsv_gz(data_folder):
    return os.path.join(data_folder, 'homology',
                        'human-chicken-orthologs-subset.tsv.gz')

#
# class TestHomologyTable:
#
#     pass


def test_cli_tsv_all_pairs(homologues_tsv):
    from khtools.homology import cli

    runner = CliRunner()
    result = runner.invoke(cli, ['--n-background', "2",
                                 "human", "chicken", homologues_tsv])
    assert result.exit_code == 0


def test_cli_tsv_n_subset5(homologues_tsv):
    from khtools.homology import cli

    runner = CliRunner()
    result = runner.invoke(cli, ['--n-subset', '20', "--n-background", "2",
                                 "human", "chicken", homologues_tsv])
    assert result.exit_code == 0


def test_cli_tsv_gz(homologues_tsv_gz):
    from khtools.homology import cli

    runner = CliRunner()
    result = runner.invoke(cli, ['--n-subset', '10',
                                 "human", "chicken", homologues_tsv_gz])
    assert result.exit_code == 0

