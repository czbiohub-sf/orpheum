import requests
import sys


def maybe_get_cds(transcript_id):
    try:
        return get_sequence(transcript_id, type='cds', verbose=False)
    except requests.exceptions.HTTPError:
        return None

def get_rna_sequence_from_protein_id(protein_id, verbose=False, type='cdna'):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{protein_id}?content-type=application/json"

    r = requests.get(server + ext,
                     headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    # print(repr(decoded))
    if verbose:
        pprint(decoded)

    transcript_id = decoded['Parent']
    sequence = get_sequence(transcript_id, type)
    return sequence


def get_sequence(ensembl_id, verbose=False):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{ensembl_id}"

    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.text
    if verbose:
        print("first 1000 residues:", decoded[:1000])
    return decoded


def get_orthologues(ensembl_id, target_species, verbose=False):
    server = "https://rest.ensembl.org"
    ext = f"/homology/id/{ensembl_id}?target_species={target_species};type=orthologues"

    r = requests.get(server + ext,
                     headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    # print(repr(decoded))
    if verbose:
        pprint(decoded)
    return decoded


def lookup(ensembl_id, verbose=False):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}"

    r = requests.get(server + ext,
                     headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    # print(repr(decoded))
    if verbose:
        pprint(decoded)
    return decoded
