import json
import os
import pickle
import re

from httmock import all_requests, HTTMock
import pytest


def make_pickled_test_data(ensembl_id, responses):
    """

    Parameters
    ----------
    ensembl_id : str
        ENSEMBL identifier to write down
    responses : dict
        Mapping of response names to the data, e.g. "sequence_response" or
        "lookup_response"

    Returns
    -------

    """
    import time
    from datetime import datetime
    import pickle

    ts = time.time()
    time_generated = datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S")

    d = {}
    for key, value in responses.items():
        d[key] = value
    d["time_generated"] = time_generated
    d["status_code"] = 200

    folder = os.path.join(os.path.dirname(__file__), "data")
    filename = os.path.join(folder, f"{ensembl_id}.pkl")
    with open(filename, "wb") as f:
        pickle.dump(d, f)


@pytest.fixture
def ensembl_protein_id():
    return "ENSP00000354687"


@pytest.fixture
def ensembl_transcript_id():
    return "ENST00000361390"


@pytest.fixture
def ensembl_gene_id():
    return "ENSG00000198888"


@all_requests
def ensembl_mock(url, request):
    match = re.match(r".*/(.*?)$", url.path)

    if match is None:
        raise ValueError(
            f"URL {url} doesn't end in an Ensembl ID. Is this a "
            "valid Ensembl REST URL?"
        )

    ensembl_id = match.group(1)

    try:
        filename = os.path.join(
            pytest.config.rootdir, f"sencha/tests/data/{ensembl_id}.pkl"
        )
        with open(filename, "rb") as testing_file:
            testing_data = pickle.load(testing_file)
    except FileNotFoundError:
        raise ValueError(f"Testing file for {ensembl_id} not found")

    print(f"Loaded test data for {ensembl_id}")
    print(f"Test data generated on {testing_data['time_generated']}")

    if "sequence" in url.path:
        return testing_data["sequence_response"]
    elif "lookup" in url.path:
        return testing_data["lookup_response"]
    elif "homology" in url.path:
        return testing_data["homology_response"]
    else:
        raise NotImplementedError(f"Endpoint {url.path} not implemented for " "mocking")


@pytest.mark.skip
def test_lookup(ensembl_protein_id):
    from sencha.ensembl import lookup

    with HTTMock(ensembl_mock):
        test = lookup(ensembl_protein_id)
        s = '{"start":3307,"Parent":"ENST00000361390","end":4262,"species":"homo_sapiens","id":"ENSP00000354687","object_type":"Translation","db_type":"core","length":318}'
        true = json.loads(s)
        assert test == true


@pytest.mark.skip
def test_lookup_expand_true(ensembl_transcript_id):
    from sencha.ensembl import lookup

    with HTTMock(ensembl_mock):
        test = lookup(ensembl_transcript_id, expand=True)
        s = """{"Exon": [{"assembly_name": "GRCh38",
           "db_type": "core",
           "end": 4262,
           "id": "ENSE00001435714",
           "object_type": "Exon",
           "seq_region_name": "MT",
           "species": "homo_sapiens",
           "start": 3307,
           "strand": 1,
           "version": 2}],
 "Parent": "ENSG00000198888",
 "Translation": {"Parent": "ENST00000361390",
                 "db_type": "core",
                 "end": 4262,
                 "id": "ENSP00000354687",
                 "length": 318,
                 "object_type": "Translation",
                 "species": "homo_sapiens",
                 "start": 3307},
 "assembly_name": "GRCh38",
 "biotype": "protein_coding",
 "db_type": "core",
 "display_name": "MT-ND1-201",
 "end": 4262,
 "id": "ENST00000361390",
 "is_canonical": 1,
 "logic_name": "mt_genbank_import_homo_sapiens",
 "object_type": "Transcript",
 "seq_region_name": "MT",
 "source": "ensembl",
 "species": "homo_sapiens",
 "start": 3307,
 "strand": 1,
 "version": 2}"""
        true = json.loads(s)
        assert test == true


@pytest.mark.skip
def test_get_orthologues(ensembl_gene_id):
    from sencha.ensembl import get_orthologues

    with HTTMock(ensembl_mock):
        test = get_orthologues(ensembl_gene_id, "mouse")
        s = """{"data": [{"id": "ENSG00000198888",
   "homologies": [{"method_link_type": "ENSEMBL_ORTHOLOGUES",
     "type": "ortholog_one2one",
     "source": {"id": "ENSG00000198888",
      "perc_id": 77.044,
      "align_seq": "MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYGLLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGLLFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYEVTLAIILLSTLLMSGSFNLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTAYPRFRYDQLMHLLWKNFLPLTLALLMWYVSMPITISSIPPQT",
      "taxon_id": 9606,
      "species": "homo_sapiens",
      "perc_pos": 88.6792,
      "cigar_line": "318M",
      "protein_id": "ENSP00000354687"},
     "target": {"species": "mus_musculus",
      "perc_pos": 88.6792,
      "cigar_line": "318M",
      "protein_id": "ENSMUSP00000080991",
      "perc_id": 77.044,
      "id": "ENSMUSG00000064341",
      "align_seq": "MFFINILTLLVPILIAMAFLTLVERKILGYMQLRKGPNIVGPYGILQPFADAMKLFMKEPMRPLTTSMSLFIIAPTLSLTLALSLWVPLPMPHPLINLNLGILFILATSSLSVYSILWSGWASNSKYSLFGALRAVAQTISYEVTMAIILLSVLLMNGSYSLQTLITTQEHMWLLLPAWPMAMMWFISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFALFFMAEYTNIILMNALTTIIFLGPLYYINLPELYSTNFMMEALLLSSTFLWIRASYPRFRYDQLMHLLWKNFLPLTLALCMWHISLPIFTAGVPPYM",
      "taxon_id": 10090},
     "taxonomy_level": "Euarchontoglires",
     "dn_ds": null}]}]}"""
        true = json.loads(s)
        assert test == true
