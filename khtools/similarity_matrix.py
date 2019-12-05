import pandas as pd


class SimilarityMatrix:

    def __init__(self):
        pass

    @classmethod
    def from_csv(csv):
        similarities = pd.read_csv(pd)
        similarities.index = similarities.columns
        return SimilarityMatrix(similarities)
