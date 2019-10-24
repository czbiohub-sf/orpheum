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
    "N": 'p',
    "C": 'p',
    "S": "p",
    "T": "p",
    "D": "p",
    "E": "p",
    "R": "p",
    "H": "p",
    "K": "p",
    "Q": "p"
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
