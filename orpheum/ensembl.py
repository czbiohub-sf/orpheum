import logging
from pprint import pprint
import sys

import requests


# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def maybe_get_cds(transcript_id):
    try:
        return get_sequence(transcript_id, type="cds", verbose=False)
    except requests.exceptions.HTTPError:
        return None


def get_rna_sequence_from_protein_id(
    protein_id, ignore_errors=False, verbose=False, type="cdna"
):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{protein_id}?content-type=application/json"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        if ignore_errors:
            logger.warning(
                f"{protein_id} was not found, likely deprecated. " "Skipping ..."
            )
            return None
        else:
            r.raise_for_status()
            sys.exit()

    decoded = r.json()
    # print(repr(decoded))
    if verbose:
        pprint(decoded)

    transcript_id = decoded["Parent"]
    sequence = get_sequence(transcript_id, type)
    return sequence


def get_sequence(ensembl_id, type=None, ignore_errors=True, verbose=False):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{ensembl_id}"
    if type is not None:
        ext += f"?type={type}"

    r = requests.get(server + ext, headers={"Content-Type": "text/plain"})

    if not r.ok:
        if ignore_errors:
            logger.warning(
                f"{ensembl_id} was not found, likely deprecated. " "Skipping ..."
            )
            return None
        else:
            r.raise_for_status()
            sys.exit()

    decoded = r.text
    if verbose:
        print("first 1000 residues:", decoded[:1000])
    return decoded


def get_orthologues(ensembl_id, target_species, verbose=False):
    server = "https://rest.ensembl.org"
    ext = (
        f"/homology/id/{ensembl_id}?target_species={target_species}" ";type=orthologues"
    )

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    # print(repr(decoded))
    if verbose:
        pprint(decoded)
    return decoded


def lookup(ensembl_id, expand=False, verbose=False):
    # Intify the boolean to 0 (False) and 1 (True)
    expand = int(expand)
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?expand={expand}"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    if verbose:
        pprint(decoded)
    return decoded
