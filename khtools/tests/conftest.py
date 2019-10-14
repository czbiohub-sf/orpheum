import os

import pytest


"""
conftest.py contains fixtures or functions-turned-variables that can be
used in any test
"""


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


@pytest.fixture(params=['protein', 'dayhoff', 'hydrophobic-polar'])
def molecule(request):
    return request.param


@pytest.fixture(params=[7, 8])
def peptide_ksize(request):
    return request.param
