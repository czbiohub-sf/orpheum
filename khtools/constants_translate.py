import itertools
from collections import namedtuple

from khtools.sequence_encodings import ALPHABET_ALIASES

DEFAULT_JACCARD_THRESHOLD = 0.5
DEFAULT_HP_JACCARD_THRESHOLD = 0.8
SEQTYPE_TO_ANNOUNCEMENT = {
    "noncoding_nucleotide":
        "nucleotide sequence from reads WITHOUT matches to "
        "protein-coding peptides",
    "coding_nucleotide":
        "nucleotide sequence from reads WITH protein-coding translation"
        " frame nucleotides",
    "low_complexity_nucleotide":
        "nucleotide sequence from low complexity (low entropy) reads",
    "low_complexity_peptide":
        "peptide sequence from low "
        "complexity (low entropy) translated"
        " reads"
}
SCORING_DF_COLUMNS = [
    'read_id', 'jaccard_in_peptide_db', 'n_kmers', 'classification'
]

LOW_COMPLEXITY_PER_ALIAS = [
    list((
        alias,
        "Low complexity peptide in {} alphabet".format(alphabet))
        for alias in aliases) for alphabet,
    aliases in ALPHABET_ALIASES.items()]
LOW_COMPLEXITY_CATEGORIES = dict(
    list(itertools.chain(*LOW_COMPLEXITY_PER_ALIAS)))


PROTEIN_CODING_CATEGORIES = {
    "too_short_peptide":
        "All translations shorter than peptide k-mer size + 1",
    "stop_codons": "All translation frames have stop codons",
    "coding": "Coding",
    'non_coding': 'Non-coding',
    'low_complexity_nucleotide': "Low complexity nucleotide",
    'too_short_nucleotide':
        'Read length was shorter than 3 * peptide k-mer size'
}

SingleReadScore = namedtuple(
    "SingleReadScore",
    ['max_fraction_kmers_in_peptide_db',
     'max_n_kmers',
     'special_case'])

COMPLEXITY_THRESHOLD = 0.3
