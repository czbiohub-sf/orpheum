import click

# khmer Nodegraph features
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = int(1e8)

# Cribbed from https://click.palletsprojects.com/en/7.x/parameters/
class BasedIntParamType(click.ParamType):
    name = "integer"

    def convert(self, value, param, ctx):
        try:
            if isinstance(value, int):
                return value
            if "e" in value:
                sigfig, exponent = value.split("e")
                sigfig = float(sigfig)
                exponent = int(exponent)
                return int(sigfig * 10 ** exponent)
            return int(value, 10)
        except TypeError:
            self.fail(
                "expected string for int() conversion, got "
                f"{value!r} of type {type(value).__name__}",
                param,
                ctx,
            )
        except ValueError:
            self.fail(f"{value!r} is not a valid integer", param, ctx)


BASED_INT = BasedIntParamType()

# '*' = stop codon
# 'X' = unknown amino acid
# 'U' = Selenocystine amino acid
RESIDUES_TO_IGNORE = "*", "X", "U"

MAX_FRACTION_OBSERVED_TO_THEORETICAL_KMERS = 1e-4

# Maximum collisions allowed in the bloom filter
MAX_BF_FALSE_POSITIVES = 0.2


# Default k-mer sizes for different alphabets
DEFAULT_PROTEIN_KSIZE = 9
DEFAULT_DAYHOFF_KSIZE = 15
DEFAULT_HP_KSIZE = 39


BEST_KSIZES = {
    "protein": DEFAULT_PROTEIN_KSIZE,
    "peptide": DEFAULT_PROTEIN_KSIZE,
    "protein20": DEFAULT_PROTEIN_KSIZE,
    "peptide20": DEFAULT_PROTEIN_KSIZE,
    "aa20": DEFAULT_PROTEIN_KSIZE,
    "dayhoff": DEFAULT_DAYHOFF_KSIZE,
    "dayhoff6": DEFAULT_DAYHOFF_KSIZE,
    "botvinnik": 13,
    "botvinnik8": 13,
    "hydrophobic-polar": DEFAULT_HP_KSIZE,
    "hydrophobic-polar2": DEFAULT_HP_KSIZE,
    "hp": DEFAULT_HP_KSIZE,
    "hp2": DEFAULT_HP_KSIZE,
    "aa9": 13,
    "gbmr4": 20,
    "sdm12": 11,
    "hsdm17": 10,
}
