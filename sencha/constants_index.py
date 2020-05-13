import click

# khmer Nodegraph features
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = int(1e8)

# Default k-mer sizes for different alphabets
DEFAULT_PROTEIN_KSIZE = 7
DEFAULT_DAYHOFF_KSIZE = 12
DEFAULT_HP_KSIZE = 31


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
