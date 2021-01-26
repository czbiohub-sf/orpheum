import pandas as pd


def combine_cell_ontology_free_annotation(row):
    if pd.notnull(row["free_annotation"]):
        return "{cell_ontology_class} ({free_annotation})".format(**row)
    else:
        return row["cell_ontology_class"]


def extract_cell_metadata(name_column, pattern=r"(?P<column>\w+):(?P<value>[\w-]+)"):
    expanded = name_column.str.extractall(pattern)
    expanded_index = expanded.reset_index()
    annotations = expanded_index.pivot(
        columns="column", values="value", index="level_0"
    )
    annotations["cell_ontology_free_annotation"] = annotations.apply(
        combine_cell_ontology_free_annotation, axis=1
    )
    return annotations


def to_key_value_pair(attribute):
    if len(attribute) > 1:
        try:
            return attribute[0], int(attribute[1])
        except ValueError:
            return attribute[0], attribute[1]
    else:
        return "comparison_sequence", attribute[0]


def extract_experiment_metadata(basename):
    key = basename.split(".csv")[0]
    split = key.split("_")
    attributes = [x.split("=") for x in split]
    attributes = dict(to_key_value_pair(x) for x in attributes)
    return key, attributes
