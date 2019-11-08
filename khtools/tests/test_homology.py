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

@pytest.fixture
def homology_scores_allpairs_csv(data_folder):
    return os.path.join(data_folder, 'homology',
                        'human-chicken-orthologs-subset-kmer-homology-scores'
                        '__allpairs_nbackground2.csv')


@pytest.fixture
def homology_scores_allpairs(homology_scores_allpairs_csv):
    return pd.read_csv(homology_scores_allpairs_csv)

#
# class TestHomologyTable:
#
#     pass


def test_cli_tsv_all_pairs(homologues_tsv, homology_scores_allpairs):
    from khtools.homology import cli

    runner = CliRunner()
    result = runner.invoke(cli, ['--n-background', "2",
                                 "human", "chicken", homologues_tsv])
    assert result.exit_code == 0
    assert homology_scores_allpairs.to_csv() == result.stdout


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

