import csv
import itertools
import json
from collections import Counter

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from sencha.constants_translate import (
    LOW_COMPLEXITY_CATEGORIES,
    PROTEIN_CODING_CATEGORIES,
    SCORING_DF_COLUMNS,
)
from sencha.log_utils import get_logger

logger = get_logger(__file__)


class CreateSaveSummary:
    def __init__(
        self,
        filenames,
        csv,
        parquet,
        json_summary,
        peptide_bloom_filter_filename,
        alphabet,
        peptide_ksize,
        jaccard_threshold,
        coding_scores,
    ):
        self.filenames = filenames
        self.csv = csv
        self.parquet = parquet
        self.json_summary = json_summary
        self.peptide_bloom_filter_filename = peptide_bloom_filter_filename
        self.alphabet = alphabet
        self.peptide_ksize = peptide_ksize
        self.jaccard_threshold = jaccard_threshold
        self.coding_scores = coding_scores
        if self.coding_scores != []:
            (
                self.read_ids,
                self.jaccard_in_peptide_dbs,
                self.n_kmers,
                self.categories,
                self.translation_frames,
                self.filenames,
            ) = map(list, zip(*self.coding_scores))

    def maybe_write_csv(self):
        if self.csv:
            logger.info("Writing coding scores of reads to {}".format(self.csv))
            # writing to csv file
            with open(self.csv, "w") as csvfile:
                # creating a csv writer object
                csvwriter = csv.writer(csvfile, lineterminator="\n")

                # writing the fields
                csvwriter.writerow(SCORING_DF_COLUMNS)

                # writing the data rows
                csvwriter.writerows(self.coding_scores)

    def maybe_write_parquet(self):
        if self.parquet:
            logger.info("Writing coding scores of reads to {}".format(self.parquet))
            batch = pa.RecordBatch.from_arrays(
                [
                    self.read_ids,
                    self.jaccard_in_peptide_dbs,
                    self.n_kmers,
                    self.categories,
                    self.translation_frames,
                    self.filenames,
                ],
                names=SCORING_DF_COLUMNS,
            )
            pq.write_table(pa.Table.from_batches([batch]), self.parquet)

    def make_empty_coding_categories(self):
        coding_categories = dict.fromkeys(PROTEIN_CODING_CATEGORIES.values(), 0)
        molecule_low_complexity_key = LOW_COMPLEXITY_CATEGORIES[self.alphabet]
        coding_categories[molecule_low_complexity_key] = 0
        return coding_categories

    def maybe_write_json_summary(self):
        if not self.json_summary:
            # Early exit if json_summary is not True
            return
        empty_coding_categories = self.make_empty_coding_categories()

        if self.coding_scores == []:
            summary = {
                "input_files": self.filenames,
                "jaccard_info": {
                    "count": 0,
                    "mean": None,
                    "std": None,
                    "min": None,
                    "25%": None,
                    "50%": None,
                    "75%": None,
                    "max": None,
                },
                "categorization_counts": empty_coding_categories,
                "categorization_percentages": empty_coding_categories,
                "histogram_n_coding_frames_per_read": {
                    str(i): 0 for i in range(len(empty_coding_categories))
                },
                "histogram_n_coding_frames_per_read_percentages": {
                    str(i): 0 for i in range(len(empty_coding_categories))
                },
            }
        else:
            summary = self.generate_coding_summary()
        with open(self.json_summary, "w") as f:
            logger.info("Writing translate summary to {}".format(self.json_summary))
            json.dump(summary, fp=f)
        # Delete these attributes once the summary is set
        # For large csv files containing lots of coding_scores use a lot of RAM
        if self.coding_scores != []:
            del self.read_ids
            del self.jaccard_in_peptide_dbs
            del self.n_kmers
            del self.categories
            del self.translation_frames
            del self.filenames
            del self.coding_scores
        return summary

    def generate_coding_summary(self):
        (
            translation_frame_percentages,
            translation_frame_counts,
        ) = self.get_n_translated_frames_per_read()

        files = np.unique(self.filenames).tolist()

        (
            categorization_percentages,
            categorization_counts,
        ) = self.get_n_per_coding_category()
        # Get Jaccard distributions, count, min, max, mean, stddev, median
        jaccard_info = {
            "count": int(np.count_nonzero(~np.isnan(self.jaccard_in_peptide_dbs))),
            "mean": np.nanmean(self.jaccard_in_peptide_dbs),
            "std": np.nanstd(self.jaccard_in_peptide_dbs),
            "min": np.nanmin(self.jaccard_in_peptide_dbs),
            "25%": np.nanpercentile(self.jaccard_in_peptide_dbs, 25),
            "50%": np.nanpercentile(self.jaccard_in_peptide_dbs, 50),
            "75%": np.nanpercentile(self.jaccard_in_peptide_dbs, 75),
            "max": np.nanmax(self.jaccard_in_peptide_dbs),
        }

        summary = {
            "input_files": files,
            "jaccard_info": jaccard_info,
            "categorization_counts": categorization_counts,
            "categorization_percentages": categorization_percentages,
            "histogram_n_coding_frames_per_read": translation_frame_counts,
            "histogram_n_coding_frames_per_read_percentages": translation_frame_percentages,
            "peptide_bloom_filter": self.peptide_bloom_filter_filename,
            "peptide_alphabet": self.alphabet,
            "peptide_ksize": self.peptide_ksize,
            "jaccard_threshold": self.jaccard_threshold,
        }
        return summary

    def get_n_per_coding_category(self):
        # Initialize to all zeros
        counts = self.make_empty_coding_categories()
        read_id_category = [
            (read_id, category)
            for read_id, category in zip(self.read_ids, self.categories)
        ]

        for read_id, categories_for_read_id in itertools.groupby(
            read_id_category, key=lambda x: x[0]
        ):
            categories_for_read_id = [
                read_category[1] for read_category in list(categories_for_read_id)
            ]
            unique_categories = np.unique(categories_for_read_id).tolist()

            # If any of the frames is coding, then that read is called coding,
            # the remaining are non-coding, unless the reads too short
            # in nucleotide or peptide space.
            if "Coding" in unique_categories:
                counts["Coding"] += 1
            elif (
                "Translation is shorter than peptide k-mer size + 1"
                in unique_categories
            ):
                counts["Translation is shorter than peptide k-mer size + 1"] += 1
            elif (
                "Read length was shorter than 3 * peptide k-mer size"
                in unique_categories
            ):
                counts["Read length was shorter than 3 * peptide k-mer size"] += 1
            elif len(unique_categories) == 1:
                counts[unique_categories[0]] += 1
            else:
                counts["Non-coding"] += 1
        total = sum(list(counts.values()))
        percentages = {
            category: 100 * count / total for category, count in counts.items()
        }
        # Replace with observations
        return percentages, counts

    def get_n_translated_frames_per_read(self):
        """Of all coding sequences, get number of possible translations"""
        predicted_coding = [
            self.read_ids[index]
            for index, category in enumerate(self.categories)
            if category == "Coding"
        ]
        n_coding_per_read = Counter(predicted_coding)

        coding_per_read_histogram = Counter(n_coding_per_read.values())
        total = sum(list(coding_per_read_histogram.values()))
        histogram_for_json = {
            "Number of reads with {} putative protein-coding translations".format(
                key
            ): value
            for key, value in coding_per_read_histogram.items()
        }
        percentages_for_json = {
            "Number of reads with {} putative protein-coding translations".format(
                key
            ): 100
            * value
            / total
            for key, value in coding_per_read_histogram.items()
        }

        return percentages_for_json, histogram_for_json
