import logging

import pandas as pd
from tqdm import tqdm

from .compare_kmer_content import compare_all_seqs
from .ensembl import get_sequence, get_rna_sequence_from_protein_id

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

QUANTITATIVE_KEYWORDS = {
    "Gene-order conservation score",
    "alignment coverage",
    "dN with",
    "dS with",
    "%id",
}
ORTHOLOGY_ORDER = [
    "No homology",
    "ortholog_one2one",
    "ortholog_one2many",
    "ortholog_many2many",
]

# Colors from seaborn's default color_palette()
ORTHOLOGY_PALETTE = dict(
    zip(
        ORTHOLOGY_ORDER,
        [
            "grey",
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
        ],
    )
)


class HomologyTable:
    def __init__(self, data, species1, species2):
        """

        Parameters
        ----------
        data : pandas.DataFrame
            Homology table from ENSEMBL Biomart (all columns)
        species1 : str
            Common name for species1, e.g. "mouse"
        species2 : str
            Common name for species2, e.g. "fly"
        """
        self.data = data
        self.species1 = species1
        self.species2 = species2

        # Extract column for homology type (e.g. one2one, one2many, many2many)
        self.homology_type_col = [
            x for x in self.data.columns if x.endswith("homology type")
        ][0]

        self.data["is_homologue"] = self.data[self.homology_type_col].notnull()

        self.species1_id_col = "Query protein or transcript ID"
        self.species2_id_col = [
            x
            for x in self.data.columns
            if x.endswith("protein or transcript stable ID")
        ][0]
        self.quantitative_features = [
            x
            for x in self.data.columns
            if any(keyword in x for keyword in QUANTITATIVE_KEYWORDS)
        ]

        gene_order_col = [
            col
            for col in self.quantitative_features
            if "Gene-order conservation score" in col
        ][0]
        self._protein_coding_rows = self.data[gene_order_col].notnull()
        self.protein_coding = self.data.loc[self._protein_coding_rows]
        self.non_coding = self.data.loc[~self._protein_coding_rows]

    @staticmethod
    def get_sequences_from_ids(df, id_column, moltype, seqtype):

        # ignore_errors=True skips deprecated IDs
        if moltype == "protein" and seqtype != "protein":
            seqs = [
                get_rna_sequence_from_protein_id(x, type=seqtype, ignore_errors=True)
                for x in tqdm(df[id_column])
            ]
        else:
            seqs = [get_sequence(x, ignore_errors=True) for x in tqdm(df[id_column])]
        # Sanitize output based on deprecated ENSEMBL IDs that don't have
        # sequences
        id_seqs = [(ID, seq) for ID, seq in zip(df[id_column], seqs) if seq is not None]
        return id_seqs

    def _get_cross_species(self, random_subset, kmer_comparisons):
        """Add species columns and subset when species are different"""
        id_to_species1 = pd.Series(
            self.species1, index=random_subset[self.species1_id_col]
        )
        id_to_species2 = pd.Series(
            self.species2, index=random_subset[self.species2_id_col]
        )

        id_to_species = pd.concat([id_to_species1, id_to_species2]).to_dict()

        kmer_comparisons["species1"] = kmer_comparisons.id1.map(id_to_species)
        kmer_comparisons["species2"] = kmer_comparisons.id2.map(id_to_species)
        kmer_comparisons["species_species"] = (
            kmer_comparisons.species1 + "_" + kmer_comparisons.species2
        )
        cross_species = kmer_comparisons.query("species1 != species2")
        del kmer_comparisons
        return cross_species

    def _add_orthology_metadata(self, cross_species, random_subset):
        """Join with original metadata to get homology information"""
        left_on = ["id1", "id2"]
        right_on = [self.species1_id_col, self.species2_id_col]
        cross_species_metadata = cross_species.merge(
            random_subset, left_on=left_on, right_on=right_on, how="left"
        )
        return cross_species_metadata

    def _subset_non_orthologous(self, cross_species_metadata, random_state):
        """Take random subsets of the non-homologous data and add back"""
        # Take random subsets of the non-homologous data as there's several
        # orders of magnitude more non-homologous pairs than homologous pairs
        cross_species_metadata_subset_non_homologues = cross_species_metadata.groupby(
            ["id1", "ksize", "alphabet"], as_index=False, group_keys=False
        ).apply(
            lambda x: x.loc[x["is_homologue"].isnull()].sample(
                10, random_state=random_state
            )
        )
        # Add the randomly sampled non homologous data back to the data that is
        # homologous
        cross_species_metadata_subset = pd.concat(
            [
                cross_species_metadata_subset_non_homologues,
                cross_species_metadata.query("is_homologue == True"),
            ],
            ignore_index=True,
        )

        return cross_species_metadata_subset

    def compare_orthology(
        self,
        datatype,
        n_subset=200,
        random_state=0,
        n_jobs=32,
        n_background=100,
        ksizes=list(range(2, 41)),
    ):
        """

        Parameters
        ----------
        datatype : str
            Either 'protein_coding_peptide', 'protein_coding_cdna', or
            'protein_coding_cds', 'non_coding'
        n_subset
        random_state
        n_jobs
        n_background : int
            Number of background comparisons to do, per species1 sequence
        ksizes

        Returns
        -------
        homology_jaccard : pandas.DataFrame
            Table of jaccard similarities of random subsets of proteins or
            transcripts across species1 and species2, joined with the original
            metadata table

        """

        if datatype == "protein_coding_peptide":
            data = self.protein_coding
            moltype = "protein"
            seqtype = "protein"
        elif datatype == "protein_coding_cdna":
            data = self.protein_coding
            moltype = "DNA"
            seqtype = "cdna"
        elif datatype == "protein_coding_cds":
            data = self.protein_coding
            moltype = "DNA"
            seqtype = "cds"
        elif datatype == "non_coding":
            data = self.non_coding
            moltype = "DNA"
            seqtype = "cdna"
        else:
            raise ValueError(
                "Only 'protein_coding_peptide',"
                " and 'protein_coding_cdna', 'protein_coding_"
                "cds', and 'non_coding' datatypes are accepted"
            )

        logger.info(f"datatype: {datatype}, moltype: {moltype}, " f"seqtype: {seqtype}")

        logger.info("Subsetting data")
        random_subset = data.sample(n_subset, random_state=random_state)
        logger.info("Getting sequences from IDs")
        species1_id_seqs = self.get_sequences_from_ids(
            random_subset, self.species1_id_col, moltype, seqtype
        )
        species2_id_seqs = self.get_sequences_from_ids(
            random_subset, self.species2_id_col, moltype, seqtype
        )

        logger.info("K-merizing and calculating jaccard comparisons")
        kmer_comparisons = compare_all_seqs(
            species1_id_seqs,
            species2_id_seqs,
            n_jobs,
            ksizes,
            moltype=moltype,
            n_background=n_background,
        )

        logger.info("Cleaning up k-mer comparisons for cross-species data")
        cross_species = self._get_cross_species(random_subset, kmer_comparisons)
        cross_species_metadata = self._add_orthology_metadata(
            cross_species, random_subset
        )

        cross_species_metadata_fillna = cross_species_metadata.fillna("No homology")
        return cross_species_metadata_fillna
