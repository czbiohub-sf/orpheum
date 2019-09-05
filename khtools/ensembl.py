import requests
import sys


def get_sequence(ensembl_id, verbose=False, type='cds'):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{ensembl_id}?type={type}"

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
