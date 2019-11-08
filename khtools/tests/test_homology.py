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
def species1():
    return 'human'


@pytest.fixture
def species2():
    return 'chicken'


@pytest.fixture
def homologues_tsv(data_folder):
    return os.path.join(data_folder, 'homology',
                        'human-chicken-orthologs-subset.tsv')

@pytest.fixture
def homologues_tsv_gz(data_folder):
    return os.path.join(data_folder, 'homology',
                        'human-chicken-orthologs-subset.tsv.gz')

@pytest.fixture
def homologues(homologues_tsv):
    return pd.read_csv(homologues_tsv, sep='\t')

@pytest.fixture
def homologues(homologues_tsv):
    return pd.read_csv(homologues_tsv, sep='\t')


@pytest.fixture
def homology_scores_allpairs_csv(data_folder):
    return os.path.join(data_folder, 'homology',
                        'human-chicken-orthologs-subset-kmer-homology-scores'
                        '__allpairs_nbackground2.csv')


@pytest.fixture
def homology_scores_allpairs(homology_scores_allpairs_csv):
    return pd.read_csv(homology_scores_allpairs_csv)


class TestHomologyTable:

    def test___init__(self, species1, species2, homologues):
        from khtools.homology import HomologyTable

        test = HomologyTable(homologues, species1, species2)
        assert test.species1 == species1
        assert test.species2 == species2
        assert test.species1_id_col == 'Query protein or transcript ID'
        assert test.species2_id_col == \
               'Chicken protein or transcript stable ID'
        assert test.homology_type_col == 'Chicken homology type'
        assert test.quantitative_features == \
               ['%id. target Chicken gene identical to query gene',
                '%id. query gene identical to target Chicken gene',
                'Chicken Gene-order conservation score',
                'Chicken Whole-genome alignment coverage',
                'dN with Chicken', 'dS with Chicken']
        assert test.protein_coding.shape == (11, 21)
        assert test.non_coding.shape == (9, 21)


def test_cli_tsv_all_pairs(species1, species2, homologues_tsv, homology_scores_allpairs):
    from khtools.homology import cli

    runner = CliRunner()
    result = runner.invoke(cli, ['--n-background', "2",
                                 species1, species2, homologues_tsv])
    assert result.exit_code == 0
    assert homology_scores_allpairs.to_csv() == result.stdout


def test_cli_tsv_n_subset5(species1, species2, homologues_tsv):
    from khtools.homology import cli

    runner = CliRunner()
    result = runner.invoke(cli, ['--n-subset', '20', "--n-background", "2",
                                 species1, species2, homologues_tsv])
    assert result.exit_code == 0


def test_cli_tsv_gz(species1, species2, homologues_tsv_gz):
    from khtools.homology import cli

    runner = CliRunner()
    result = runner.invoke(cli, ['--n-subset', '10',
                                 species1, species2, homologues_tsv_gz])
    assert result.exit_code == 0

