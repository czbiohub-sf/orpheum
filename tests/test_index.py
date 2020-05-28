from click.testing import CliRunner
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
        ("dayhoff", DEFAULT_PROTEIN_KSIZE),
        ("hydrophobic-polar", DEFAULT_HP_KSIZE),
        ("hydrophobic-polar", DEFAULT_PROTEIN_KSIZE)
    ],
    ids=[
        "protein_default_ksize",
        "dayhoff_default_ksize",
        "dayhoff_protein_ksize_xfail",
        "hp_default_ksize",
        "hp_protein_ksize_xfail",
    ],
)
def alphabet_ksize_index(request):
    return request.param


@pytest.fixture
def peptide_ksize_index(alphabet_ksize):
    return alphabet_ksize[1]


@pytest.fixture
def alphabet_index(alphabet_ksize):
    return alphabet_ksize[0]


def test_make_peptide_index(variable_peptide_fasta,
                            alphabet_index,
                            peptide_ksize_index):
    from sencha.index import make_peptide_index

    test = make_peptide_index(
        variable_peptide_fasta, peptide_ksize_index, alphabet_index,
        n_tables=4, tablesize=1e6
    )
    if "first1000lines" in variable_peptide_fasta:
        TRUE_N_UNIQUE_KMERS = {
            ("protein", 7): 13966,
            ("dayhoff", 7): 10090,
            ("dayhoff", 11): 13605,
            ("dayhoff", 12): 13816,
            ("hydrophobic-polar", 31): 12888,
            ("hydrophobic-polar", 21): 13076,
            ("hydrophobic-polar", 7): 136,
        }
    else:
        TRUE_N_UNIQUE_KMERS = {
            ("protein", 7): 506352,
            ("dayhoff", 7): 99863,
            ("dayhoff", 11): 472197,
            ("dayhoff", 12): 488469,
            ("hydrophobic-polar", 31): 515863,
            ("hydrophobic-polar", 21): 434810,
            ("hydrophobic-polar", 7): 170,
        }
    true_n_unique_kmers = TRUE_N_UNIQUE_KMERS[(alphabet_index, peptide_ksize_index)]

    # For now, assert that the number of kmers is within 0.1% of the true value
    np.testing.assert_allclose(test.n_unique_kmers(), true_n_unique_kmers, rtol=0.001)


def test_error_if_index_tables_too_small(adversarial_peptide_fasta,):
    from sencha.index import make_peptide_index
    with pytest.raises(SystemExit) as pytest_wrapped_error:
        make_peptide_index(
            adversarial_peptide_fasta, peptide_ksize=9, molecule='protein',
            n_tables=2, tablesize=1e2
        )
        assert pytest_wrapped_error.type == SystemExit
        assert pytest_wrapped_error.value.code == 42


def test_error_if_too_many_observed_kmers(peptide_fasta,):
    from sencha.index import make_peptide_index
    with pytest.raises(ValueError):
        make_peptide_index(
            peptide_fasta, peptide_ksize=7, molecule='dayhoff',
            n_tables=4, tablesize=1e6,
        )


def test_maybe_make_peptide_index(
    peptide_bloom_filter_path, alphabet, peptide_ksize
):
    from sencha.index import maybe_make_peptide_index

    maybe_make_peptide_index(
        peptide_bloom_filter_path,
        peptide_ksize,
        alphabet,
        peptides_are_bloom_filter=True,
    )
    # No assertion, just check that it ran
    # assert isinstance(test, khmer.Nodegraph)


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
            "1e4",
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
