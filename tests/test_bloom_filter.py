from click.testing import CliRunner
import pytest


def test_make_peptide_bloom_filter(variable_peptide_fasta,
                                   molecule, peptide_ksize):
    from khtools.bloom_filter import make_peptide_bloom_filter

    test = make_peptide_bloom_filter(variable_peptide_fasta,
                                     peptide_ksize,
                                     molecule,
                                     n_tables=4,
                                     tablesize=1e6)
    if 'first1000lines' in variable_peptide_fasta:
        TRUE_N_UNIQUE_KMERS = {
            ("protein", 7): 13966,
            ("dayhoff", 7): 10090,
            ("dayhoff", 11): 13605,
            ("hydrophobic-polar", 21): 13076,
            ("hydrophobic-polar", 7): 136
        }
    else:
        TRUE_N_UNIQUE_KMERS = {
            ("protein", 7): 506352,
            ("dayhoff", 7): 99863,
            ("dayhoff", 11): 472197,
            ("dayhoff", 12): 488469,
            ("hydrophobic-polar", 21): 434810,
            ("hydrophobic-polar", 7): 170
        }
    true_n_unique_kmers = TRUE_N_UNIQUE_KMERS[(molecule, peptide_ksize)]

    # For now, assert that the number of kmers is within 5% of the true value
    assert test.n_unique_kmers() > true_n_unique_kmers * 0.95
    assert test.n_unique_kmers() < true_n_unique_kmers * 1.05


def test_maybe_make_peptide_bloom_filter(peptide_bloom_filter_path,
                                         molecule, peptide_ksize):
    from khtools.bloom_filter import maybe_make_peptide_bloom_filter

    maybe_make_peptide_bloom_filter(peptide_bloom_filter_path,
                                    peptide_ksize,
                                    molecule,
                                    peptides_are_bloom_filter=True)
    # No assertion, just check that it ran
    # assert isinstance(test, khmer.Nodegraph)


def test_cli_minimum(peptide_fasta):
    from khtools.bloom_filter import cli

    runner = CliRunner()
    result = runner.invoke(cli, [
        peptide_fasta,
    ])
    assert result.exit_code == 0


def test_cli_options(peptide_fasta, molecule, peptide_ksize):
    from khtools.bloom_filter import cli

    runner = CliRunner()
    result = runner.invoke(cli, [
        '--peptide-ksize', peptide_ksize, '--molecule', molecule,
        "--tablesize", "1e4",
        peptide_fasta,
    ])
    assert result.exit_code == 0


def test_get_peptide_ksize_default(molecule):
    from khtools.bloom_filter import get_peptide_ksize, \
        DEFAULT_PROTEIN_KSIZE, DEFAULT_HP_KSIZE, DEFAULT_DAYHOFF_KSIZE

    test = get_peptide_ksize(molecule, peptide_ksize=None)
    if molecule == 'protein':
        assert test == DEFAULT_PROTEIN_KSIZE
    elif molecule == 'dayhoff':
        assert test == DEFAULT_DAYHOFF_KSIZE
    elif molecule == 'hydrophobic-polar':
        assert test == DEFAULT_HP_KSIZE


def test_get_peptide_ksize_with_ksize(molecule):
    from khtools.bloom_filter import get_peptide_ksize

    peptide_ksize = 123
    test = get_peptide_ksize(molecule, peptide_ksize)
    assert test == peptide_ksize


def test_get_peptide_ksize_with_bad_molecule():
    from khtools.bloom_filter import get_peptide_ksize

    peptide_ksize = 123
    with pytest.raises(ValueError):
        get_peptide_ksize("not a real molecule type", peptide_ksize)
