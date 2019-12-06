from anndata import Anndata
import pandas as pd


class SimilarityMatrix(Anndata):

    def __init__(self, data, metadata=None):
        self.data = data
        self.metadata = metadata

    def _validate_metadata(self):
        # Check that all metadata
        assert


    @classmethod
    def from_csv(csv):
        similarities = pd.read_csv(csv)
        similarities.index = similarities.columns
        return SimilarityMatrix(similarities, filename=csv)

