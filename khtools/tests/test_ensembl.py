import json
import os
import pickle
import re

from httmock import all_requests, HTTMock
import pytest


@pytest.fixture
def ensembl_protein_id():
    return "ENSP00000354687"


@all_requests
def ensembl_mock(url, request):
    match = re.match(r'.*/(.*?)$', url.path)

    if match is None:
        raise ValueError(f"URL {url} doesn't end in an Ensembl ID. Is this a "
                         "valid Ensembl REST URL?")

    ensembl_id = match.group(1)

    try:
        filename = os.path.join(pytest.config.rootdir,
                                f'khtools/tests/data/{ensembl_id}.pkl')
        with open(filename, 'rb') as testing_file:
            testing_data = pickle.load(testing_file)
    except FileNotFoundError:
        raise ValueError(f"Testing file for transcript {ensembl_id} "
                         "not found")

    print(f"Loaded test data for {ensembl_id}")
    print(f"Test data generated on {testing_data['time_generated']}")

    if 'sequence' in url.path:
        return testing_data['sequence_response']
    elif 'lookup' in url.path:
        return testing_data['lookup_response']
    else:
        raise NotImplementedError(f"Endpoint {url.path} not implemented for "
                                  "mocking")


def test_lookup(ensembl_protein_id):
    from khtools.ensembl import lookup

    with HTTMock(ensembl_mock):
        test = lookup(ensembl_protein_id)
        s = '{"start":3307,"Parent":"ENST00000361390","end":4262,"species":"homo_sapiens","id":"ENSP00000354687","object_type":"Translation","db_type":"core","length":318}'
        true = json.loads(s)
        assert test == true
