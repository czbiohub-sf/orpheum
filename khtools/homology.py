import requests
import sys

import pandas as pd

from .compare_peptide import compare_all_seqs
from .ensembl import get_sequence, maybe_get_cds

QUANTITATIVE_KEYWORDS = set(
    ['conservation score', 'alignment coverage', 'dN with', 'dS with',
     '%id'])


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
        homology_type_col = [x for x in self.data.columns
                             if x.endswith("homology type")][0]

        self.data['is_homologue'] = self.data[homology_type_col].notnull()

        self.species1_id_col = 'Query protein or transcript ID'
        self.species2_id_col = [x for x in self.data.columns if
                           x.endswith("protein or transcript stable ID")][0]
        self.quantitative_features = [x for x in self.data.columns if any(
            keyword in x for keyword in QUANTITATIVE_KEYWORDS)]

        dN_col = [col for col in self.quantitative_features if 'dN with' in col][0]
        self._protein_coding_rows = self.data[dN_col].notnull()
        self.protein_coding = self.data.loc[self._protein_coding_rows]
        self.non_coding = self.data.loc[~self._protein_coding_rows]

    def get_sequences_from_ids(self, df, id_column):
        seqs = [get_sequence(x) for x in
                         df[id_column]]
        id_seqs = list(
            zip(df[id_column], seqs))
        return id_seqs

    def _get_cross_species(self, random_subset, kmer_comparisons):
        """Add species columns and subset when species are different"""
        id_to_species1 = pd.Series(self.species1,
                                   index=random_subset[self.species1_id_col])
        id_to_species2 = pd.Series(self.species2,
                                   index=random_subset[self.species2_id_col])

        id_to_species = pd.concat([id_to_species1, id_to_species2]).to_dict()

        kmer_comparisons[
            'species1'] = kmer_comparisons.id1.map(id_to_species)
        kmer_comparisons[
            'species2'] = kmer_comparisons.id2.map(id_to_species)
        kmer_comparisons[
            'species_species'] = kmer_comparisons.species1 + "_" + kmer_comparisons.species2
        cross_species = kmer_comparisons.query('species1 != species2')
        del kmer_comparisons
        return cross_species

    def _add_orthology_metadata(self, cross_species, random_subset):
        """Join with original metadata to get homology information"""
        left_on = ['id1', 'id2']
        right_on = [self.species1_id_col, self.species2_id_col]
        cross_species_metadata = cross_species.merge(random_subset,
                                                     left_on=left_on,
                                                     right_on=right_on,
                                                     how='left')
        return cross_species_metadata

    def _subset_non_orthologous(self, cross_species_metadata, random_state):
        """Take random subsets of the non-homologous data and add back"""
        # Take random subsets of the non-homologous data as there's several
        # orders of magnitude more non-homologous pairs than homologous pairs
        cross_species_metadata_subset_non_homologues = cross_species_metadata.groupby(
            ['id1', 'ksize', 'molecule'], as_index=False,
            group_keys=False).apply(
            lambda x: x.loc[x['is_homologue'].isnull()].sample(
                10, random_state=random_state))
        # Add the randomly sampled non homologous data back to the data that is
        # homologous
        cross_species_metadata_subset = pd.concat(
            [cross_species_metadata_subset_non_homologues,
             cross_species_metadata.query('is_homologue == True')],
            ignore_index=True)

        return cross_species_metadata_subset

    def compare_orthology(self, datatype, n_subset=200, random_state=0,
                          n_jobs=32, ksizes=list(range(2, 41))):
        if datatype == 'protein':
            data = self.protein_coding
            moltype = 'protein'
        elif datatype == 'non_coding':
            data = self.protein_coding
            moltype = 'DNA'
        else:
            raise ValueError("Only 'protein_coding' and 'non_coding' data "
                             "types are accepted")

        random_subset = data.sample(n_subset, random_state=random_state)
        species1_id_seqs = self.get_sequences_from_ids(random_subset,
                                                        self.species1_id_col)
        species2_id_seqs = self.get_sequences_from_ids(random_subset,
                                                        self.species2_id_col)
        seqlist = species1_id_seqs + species2_id_seqs
        kmer_comparisons = compare_all_seqs(seqlist, n_jobs, ksizes,
                                            moltype=moltype)

        cross_species = self._get_cross_species(random_subset,
                                                kmer_comparisons)
        cross_species_metadata = self._add_orthology_metadata(cross_species,
                                                              random_subset)
        cross_species_metadata_subset = self._subset_non_orthologous(
            cross_species_metadata, random_state)

        cross_species_metadata_fillna = cross_species_metadata_subset.fillna(
            "No homology")
        return cross_species_metadata_fillna
