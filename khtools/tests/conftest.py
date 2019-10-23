import os

from khmer import Nodegraph
import pytest


"""
conftest.py contains fixtures or functions-turned-variables that can be
used in any test
"""
from khtools.bloom_filter import DEFAULT_PROTEIN_KSIZE, DEFAULT_DAYHOFF_KSIZE, \
    DEFAULT_HP_KSIZE


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        './data')


@pytest.fixture
def peptide_fasta(data_folder):
    filename = os.path.join(data_folder, 'bloom_filter',
                            'Homo_sapiens.GRCh38.pep.subset.fa.gz')
    return filename


# Tie the molecule name to its default ksize to make sure we keep getting the
# right sequences
@pytest.fixture(params=[('protein', DEFAULT_PROTEIN_KSIZE),
                        ('dayhoff', DEFAULT_DAYHOFF_KSIZE),
                        pytest.param(('dayhoff', DEFAULT_PROTEIN_KSIZE),
                                     marks=pytest.mark.xfail),
                        ('hydrophobic-polar', DEFAULT_HP_KSIZE),
                        pytest.param(('hydrophobic-polar', DEFAULT_PROTEIN_KSIZE),
                                     marks=pytest.mark.xfail)
                        ],
                ids=['protein_default_ksize', 'dayhoff_default_ksize',
                     'dayhoff_protein_ksize_xfail',
                     'hp_default_ksize', 'hp_protein_ksize_xfail'])
def molecule_ksize(request):
    return request.param


@pytest.fixture
def peptide_ksize(molecule_ksize):
    return molecule_ksize[1]


@pytest.fixture
def molecule(molecule_ksize):
    return molecule_ksize[0]


@pytest.fixture
def peptide_bloom_filter(data_folder, molecule, peptide_ksize):
    filename = os.path.join(data_folder, 'bloom_filter',
                            f'Homo_sapiens.GRCh38.pep.subset.molecule-{molecule}_ksize-{peptide_ksize}.bloomfilter.nodegraph')
    return Nodegraph.load(filename)



@pytest.fixture
def peptide_bloom_filter(data_folder, molecule, peptide_ksize):
    filename = os.path.join(data_folder, 'bloom_filter',
                            f'Homo_sapiens.GRCh38.pep.subset.molecule-{molecule}_ksize-{peptide_ksize}.bloomfilter.nodegraph')
    return Nodegraph.load(filename)

