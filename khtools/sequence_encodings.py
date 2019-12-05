DNA_ALPHABET = "A", "C", "G", "T"
AMINO_ACID_SINGLE_LETTERS = "R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", \
                            "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"
DAYHOFF_MAPPING = {
    "C": "a",

    # Small
    "A": "b",
    "G": "b",
    "P": "b",
    "S": "b",
    "T": "b",

    # Acid and amide
    "D": "c",
    "E": "c",
    "N": "c",
    "Q": "c",

    # Basic
    "H": "d",
    "K": "d",
    "R": "d",

    # Hydrophobic
    "I": "e",
    "L": "e",
    "M": "e",
    "V": "e",

    # Aromatic
    "F": "f",
    "W": "f",
    "Y": "f"
}
DAYHOFF_v2_MAPPING = {
    "C": "a",

    # Small
    "A": "b",
    "G": "b",
    "P": "b",

    # Phosphorylateable
    "S": "B",
    "T": "B",

    # Acid and amide
    "D": "c",
    "E": "c",
    "N": "c",
    "Q": "c",

    # Basic
    "H": "d",
    "K": "d",
    "R": "d",

    # Hydrophobic
    "I": "e",
    "L": "e",
    "M": "e",
    "V": "e",

    # Aromatic
    "F": "f",
    "W": "f",
    "Y": "f"
}
HP_MAPPING = {
    # Hydrophobic
    "A": "h",
    "F": "h",
    "G": "h",
    "I": "h",
    "L": "h",
    "M": 'h',
    "P": "h",
    "V": "h",
    "W": "h",
    "Y": "h",

    # Hydrophilic - polar
    "C": 'p',
    "D": 'p',
    "E": "p",
    "H": "p",
    "K": "p",
    "N": "p",
    "Q": "p",
    "R": "p",
    "S": "p",
    "T": "p"
}

# GBMR4, SDM12, HSDM17 from the following paper:
# Peterson, E. L., Kondev, J., Theriot, J. A., & Phillips, R. (2009).
# Reduced amino acid alphabets exhibit an improved sensitivity and
# selectivity in fold assignment. Bioinformatics, 25(11), 1356–1362.
# http://doi.org/10.1093/bioinformatics/btp164
GBMR4 = {
    # Small/polar?
    'A': 'a',
    'D': 'a',
    'K': 'a',
    'E': 'a',
    'R': 'a',
    'N': 'a',
    'T': 'a',
    'S': 'a',
    'Q': 'a',

    # Hydrophobic/large?
    'Y': 'b',
    'F': 'b',
    'L': 'b',
    'I': 'b',
    'V': 'b',
    'M': 'b',
    'C': 'b',
    'W': 'b',
    'H': 'b',

    "G": "c",

    "P": "d"
}


SDM12_MAPPING = {
    "A": "a",
    "D": "b",

    'K': 'c', 'E': 'c', 'R': 'c',

    "N": "d",

    'T': 'e', 'S': 'e', 'Q': 'e',

    'Y': 'f', 'F': 'f',

    'L': 'g', 'I': 'g', 'V': 'g', 'M': 'g',

    "C": "h",

    "W": "i",

    "H": "j",

    "G": "k",

    "P": "l"
}

HSMD17 = {
    "A": "a",
    "D": "b",
    'K': 'c', 'E': 'c',

    "R": "d",
    "N": "e",
    "T": "f",
    "S": "g",
    "Q": "h",
    "Y": 'i',
    'F': 'j',

    'L': 'k', 'I': 'k', 'V': 'k',

    "M": "l",
    "C": "m",
    "G": "n",
    "P": "o"
}

BOTVINNIK_MAPPING = {
    # Small and hydrophobic
    "A": "a",
    "G": "a",

    # Hydrophobic
    "L": "b",
    "I": "b",
    "V": "b",

    # Aromatic, not W
    "F": "c",
    "Y": "c",

    # Polar or charged
    # Phosphorylate-able
    "S": "d",
    "T": "d",

    # Polar, uncharged
    "N": "e",
    "Q": "e",

    # Polar, negatively charged
    "D": "f",
    "E": "f",

    # Polar, positively charged
    # Not histidine
    "R": "g",
    "K": "g",

    # Special
    "C": "h",
    "M": "i",
    "W": "j",
    "H": "k",
    "P": "m"
}

PURINE_PYRIMIDINE_MAPPING = {"A": "R", "C": "Y", "G": "R", "T": "Y"}
AMINO_KETO_MAPPING = {"A": "M", "C": "M", "G": "K", "T": "K"}
WEAK_STRONG_MAPPING = {"A": "W", "C": "S", "G": "S", "T": "W"}
AMINO_KETO_TRANSLATION = str.maketrans(AMINO_KETO_MAPPING)
WEAK_STRONG_TRANSLATION = str.maketrans(WEAK_STRONG_MAPPING)
PURINE_PYRIMIDINE_TRANSLATION = str.maketrans(PURINE_PYRIMIDINE_MAPPING)
DAYHOFF_TRANSLATION = str.maketrans(DAYHOFF_MAPPING)
DAYHOFF_V2_TRANSLATION = str.maketrans(DAYHOFF_v2_MAPPING)
HP_TRANSLATION = str.maketrans(HP_MAPPING)
BOTVINNIK_TRANSLATION = str.maketrans(BOTVINNIK_MAPPING)


VALID_PEPTIDE_MOLECULES = 'protein', 'peptide', 'dayhoff', \
                          'hydrophobic-polar', 'hp'
# Nucleic acid mappings


def amino_keto_ize(seq):
    return seq.translate(AMINO_KETO_TRANSLATION)


def weak_strong_ize(seq):
    return seq.translate(WEAK_STRONG_TRANSLATION)


def purine_pyrimidize(seq):
    return seq.translate(PURINE_PYRIMIDINE_TRANSLATION)


# Amino acid mappings
def dayhoffize(seq):
    return seq.translate(DAYHOFF_TRANSLATION)


def dayhoff_v2_ize(seq):
    return seq.translate(DAYHOFF_V2_TRANSLATION)


def hpize(seq):
    return seq.translate(HP_TRANSLATION)


def botvinnikize(seq):
    return seq.translate(BOTVINNIK_TRANSLATION)


def encode_peptide(peptide_sequence, molecule):
    if molecule == 'dayhoff':
        return dayhoffize(peptide_sequence)
    elif molecule == 'hydrophobic-polar' or molecule == 'hp':
        return hpize(peptide_sequence)
    elif molecule in VALID_PEPTIDE_MOLECULES:
        return peptide_sequence
    else:
        raise ValueError(f"{molecule} is not a valid amino acid encoding, "
                         "only "
                         "{', '.join(VALID_PEPTIDE_MOLECULES} can be used")
