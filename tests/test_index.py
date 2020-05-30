from click.testing import CliRunner
import khmer
import numpy as np
import pytest

from sencha.constants_index import (
    DEFAULT_PROTEIN_KSIZE,
    DEFAULT_DAYHOFF_KSIZE,
    DEFAULT_HP_KSIZE,
)


def test_per_translation_false_positive_rate():
    from sencha.index import per_translation_false_positive_rate

    n_kmers_in_translation = 14
    n_total_kmers = 4e7
    test = per_translation_false_positive_rate(n_kmers_in_translation, n_total_kmers)
    assert test == 3.275785921922898e-06


def test_per_read_false_positive_coding_rate():
    from sencha.index import per_read_false_positive_coding_rate

    read_length = 60
    peptide_ksize = 7
    test = per_read_false_positive_coding_rate(read_length, peptide_ksize)
    assert test == 3.884682410293139e-05


# Tie the alphabet name to its default ksize to make sure we keep getting the
# right sequences
@pytest.fixture(
    params=[
        ("protein", DEFAULT_PROTEIN_KSIZE),
        ("dayhoff", DEFAULT_DAYHOFF_KSIZE),
        ("hydrophobic-polar", DEFAULT_HP_KSIZE),
    ],
    ids=[
        f"protein_ksize{DEFAULT_PROTEIN_KSIZE}",
        f"dayhoff_ksize{DEFAULT_DAYHOFF_KSIZE}",
        f"hp_ksize{DEFAULT_HP_KSIZE}",
    ],
)
def alphabet_ksize_index(request):
    return request.param


@pytest.fixture
def peptide_ksize_index(alphabet_ksize_index):
    return alphabet_ksize_index[1]


@pytest.fixture
def alphabet_index(alphabet_ksize_index):
    return alphabet_ksize_index[0]


def test_make_peptide_index(
    variable_peptide_fasta, alphabet_index, peptide_ksize_index
):
    from sencha.index import make_peptide_index

    test = make_peptide_index(
        variable_peptide_fasta,
        peptide_ksize_index,
        alphabet_index,
        n_tables=4,
        tablesize=1e6,
        max_observed_fraction=1e-1,
    )
    if "first1000lines" in variable_peptide_fasta:
        # This is the adversarial fasta with short sequences
        TRUE_N_UNIQUE_KMERS = {
            ("protein", DEFAULT_PROTEIN_KSIZE): 14448,
            ("dayhoff", DEFAULT_DAYHOFF_KSIZE): 14094,
            ("hydrophobic-polar", DEFAULT_HP_KSIZE): 12244,
        }
    else:
        # This is the "normal" fasta
        TRUE_N_UNIQUE_KMERS = {
            ("protein", DEFAULT_PROTEIN_KSIZE): 517240,
            ("dayhoff", DEFAULT_DAYHOFF_KSIZE): 506018,
            ("hydrophobic-polar", DEFAULT_HP_KSIZE): 522461,
        }
    true_n_unique_kmers = TRUE_N_UNIQUE_KMERS[(alphabet_index, peptide_ksize_index)]

    # For now, assert that the number of kmers is within 0.001% of the true value
    np.testing.assert_allclose(test.n_unique_kmers(), true_n_unique_kmers, rtol=1e-4)


def test_error_if_index_tables_too_small(adversarial_peptide_fasta,):
    from sencha.index import make_peptide_index

    with pytest.raises(ValueError) as pytest_wrapped_error:
        make_peptide_index(
            adversarial_peptide_fasta,
            peptide_ksize=9,
            molecule="protein",
            n_tables=2,
            tablesize=1e2,
        )
        assert pytest_wrapped_error.type == ValueError
        assert pytest_wrapped_error.value.code == 42


def test_error_if_too_many_observed_kmers(peptide_fasta):
    from sencha.index import make_peptide_index

    with pytest.raises(ValueError):
        make_peptide_index(
            peptide_fasta,
            peptide_ksize=7,
            molecule="dayhoff",
            n_tables=4,
            tablesize=1e6,
        )


def test_maybe_make_peptide_index(peptide_bloom_filter_path, alphabet, peptide_ksize):
    from sencha.index import maybe_make_peptide_index

    test = maybe_make_peptide_index(
        peptide_bloom_filter_path, peptide_ksize, alphabet, peptides_are_index=True,
    )
    # No assertion, just check that it ran
    assert isinstance(test, khmer._oxli.graphs.Nodegraph)


def test_cli_minimum(peptide_fasta):
    from sencha.index import cli

    runner = CliRunner()
    result = runner.invoke(cli, [peptide_fasta,])
    assert result.exit_code == 0


def test_cli_options(peptide_fasta, alphabet, peptide_ksize):
    from sencha.index import cli

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--peptide-ksize",
            peptide_ksize,
            "--alphabet",
            alphabet,
            "--tablesize",
            "1e7",
            "--max-observed-fraction",
            "1e-1",
            peptide_fasta,
        ],
    )
    assert result.exit_code == 0


def test_get_peptide_ksize_default(alphabet):
    from sencha.index import get_peptide_ksize
    from sencha.constants_index import (
        DEFAULT_PROTEIN_KSIZE,
        DEFAULT_HP_KSIZE,
        DEFAULT_DAYHOFF_KSIZE,
    )

    test = get_peptide_ksize(alphabet, peptide_ksize=None)
    if alphabet == "protein":
        assert test == DEFAULT_PROTEIN_KSIZE
    elif alphabet == "dayhoff":
        assert test == DEFAULT_DAYHOFF_KSIZE
    elif alphabet == "hydrophobic-polar":
        assert test == DEFAULT_HP_KSIZE


def test_get_peptide_ksize_with_ksize(alphabet):
    from sencha.index import get_peptide_ksize

    peptide_ksize = 123
    test = get_peptide_ksize(alphabet, peptide_ksize)
    assert test == peptide_ksize


def test_get_peptide_ksize_with_bad_alphabet():
    from sencha.index import get_peptide_ksize

    with pytest.raises(ValueError):
        get_peptide_ksize("not a real alphabet type", None)
