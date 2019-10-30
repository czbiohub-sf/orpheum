import os

from khmer import Nodegraph
import pytest
"""
conftest.py contains fixtures or functions-turned-variables that can be
used in any test
"""
from khtools.bloom_filter import DEFAULT_PROTEIN_KSIZE, \
    DEFAULT_DAYHOFF_KSIZE, DEFAULT_HP_KSIZE


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), './data')


@pytest.fixture
def peptide_fasta(data_folder):
    filename = os.path.join(data_folder, 'bloom_filter',
                            'Homo_sapiens.GRCh38.pep.subset.fa.gz')
    return filename


@pytest.fixture
def adversarial_peptide_fasta(data_folder):
    filename = os.path.join(data_folder, 'bloom_filter',
                            'Homo_sapiens.GRCh38.pep.first1000lines.fa')
    return filename


@pytest.fixture(params=['normal', 'adversarial'])
def variable_peptide_fasta(request, peptide_fasta, adversarial_peptide_fasta):
    if request.param == 'normal':
        return peptide_fasta
    else:
        return adversarial_peptide_fasta


# Tie the molecule name to its default ksize to make sure we keep getting the
# right sequences
@pytest.fixture(params=[('protein', DEFAULT_PROTEIN_KSIZE),
                        ('dayhoff', DEFAULT_DAYHOFF_KSIZE),
                        pytest.param(('dayhoff', DEFAULT_PROTEIN_KSIZE),
                                     marks=pytest.mark.xfail),
                        ('hydrophobic-polar', DEFAULT_HP_KSIZE),
                        pytest.param(
                            ('hydrophobic-polar', DEFAULT_PROTEIN_KSIZE),
                            marks=pytest.mark.xfail)],
                ids=[
                    'protein_default_ksize', 'dayhoff_default_ksize',
                    'dayhoff_protein_ksize_xfail', 'hp_default_ksize',
                    'hp_protein_ksize_xfail'
])
def molecule_ksize(request):
    return request.param


@pytest.fixture
def peptide_ksize(molecule_ksize):
    return molecule_ksize[1]


@pytest.fixture
def molecule(molecule_ksize):
    return molecule_ksize[0]


@pytest.fixture
def peptide_bloom_filter_path(data_folder, molecule, peptide_ksize):
    filename = os.path.join(
        data_folder, 'bloom_filter',
        f'Homo_sapiens.GRCh38.pep.subset.molecule-{molecule}_'
        f'ksize-{peptide_ksize}.bloomfilter.nodegraph'
    )
    return filename


@pytest.fixture
def peptide_bloom_filter(peptide_bloom_filter_path, peptide_fasta, molecule,
                         peptide_ksize):
    """Load bloom filter from path if exists, otherwise, make it"""
    try:
        return Nodegraph.load(peptide_bloom_filter_path)
    except (FileNotFoundError, OSError):
        from khtools.bloom_filter import make_peptide_bloom_filter

        bloom_filter = make_peptide_bloom_filter(peptide_fasta,
                                                 peptide_ksize,
                                                 molecule,
                                                 tablesize=1e6)
        bloom_filter.save(peptide_bloom_filter_path)
        return bloom_filter
