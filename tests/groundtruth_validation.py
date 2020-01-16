import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import screed
import sklearn.metrics

from khtools.bloom_filter import make_peptide_bloom_filter
from khtools.extract_coding import score_reads

DEFAULT_CMAP = "Oranges"
DEFAULT_FONT_SIZE = 11
DEFAULT_LINE_WIDTH = 0.5
DEFAULT_FIG_SIZE = (9, 9)


def format_element_to_matlab_confusion_matrix(row, col, confusion_matrix):
    """
    Return a string for the element on row, col location for
    confusion_matrix to either
    number of observation\npercentage of observations or
    percentage_correct_classifications\npercentage_incorrect_classifications
    per class

    Args:
        row: int row location inside array
        col: int col locat
        confusion_matrix: numpy array of symmetrical (x, x) shape with target
            class/groundtruth in columns and output/predicted class in rows and
            total for each class including unmatched groundtruth
            (not detected by the model) & unmatched predicted class(due to low
            confidence)

    Returns:
        text: str string depending on row, col location for
            element in the confusion_matrix to either
            number of observation\npercentage of observations or
            percentage_correct_classifications\npercentage_incorrect
            classifications per class.
            If row, col are equal to the last row, col of the array then the
            element represents totals per class,
            percentage_correct_classifications\npercentage_incorrect
            classifications per class is returned else
            it is the number of objects per label in either ground truth or
            predicted and number of observation\npercentage of observations is
            returned

    """
    current_element = confusion_matrix[row][col]
    total_predicted = confusion_matrix[-1][-1]
    percentage_total = (float(current_element) / total_predicted) * 100
    cm_length = confusion_matrix.shape[0]

    # for totals calculate percentage accuracy and percentage error
    if(col == (cm_length - 1)) or (row == (cm_length - 1)):
        # totals and percents
        if(current_element != 0):
            if(col == cm_length - 1) and (row == cm_length - 1):
                total_correct = 0
                for i in range(confusion_matrix.shape[0] - 1):
                    total_correct += confusion_matrix[i][i]
                percentage_correct_classifications = (
                    float(total_correct) / current_element) * 100
            elif(col == cm_length - 1):
                true_positives_for_label = confusion_matrix[row][row]
                true_predicted_per_class = current_element
                percentage_correct_classifications = (
                    float(true_positives_for_label) / true_predicted_per_class
                ) * 100
            elif(row == cm_length - 1):
                true_positives_for_label = confusion_matrix[col][col]
                true_groundtruth_per_class = current_element
                percentage_correct_classifications = (
                    float(
                        true_positives_for_label) / true_groundtruth_per_class
                ) * 100
            percentage_incorrect_classifications = \
                100 - percentage_correct_classifications
        else:
            percentage_correct_classifications = \
                percentage_incorrect_classifications = 0

        percentage_correct_classifications_s = [
            '%.2f%%' % (percentage_correct_classifications), '100%'][
            percentage_correct_classifications == 100]
        txt = '%s\n%.2f%%' % (
            percentage_correct_classifications_s,
            percentage_incorrect_classifications)
    else:
        if(percentage_total > 0):
            txt = '%s\n%.2f%%' % (current_element, percentage_total)
        else:
            txt = '0\n0.0%'
    return txt


def plot_cm(confusion_matrix, labels, output_fig):
    """
    Save confusion matrix, precision, recall scores of each of
    the unique labels to a figure

    Args:
        confusion_matrix: numpy array of symmetrical (x, x) shape,
        along the columns is predicted data
        labels: list of unqiue names of the objects present
        output_fig: str output figure file containing confusion matrix,
            precision, recall per class (format can be png, pdf, eps, svg)

    Returns:
        Plots and saves matlab like confusion matrix with
        precision, recall, and total percentages to output_fig
    The rows correspond to the predicted class (Output Class) and
    the columns correspond to the true class (Target Class).
    The diagonal cells correspond to observations that are correctly
    classified. The off-diagonal cells correspond to incorrectly
    classified observations. Both the number of observations and the
    percentage of the total number of observations are shown in each cell.

    The column on the far right of the plot shows the percentages of all
    the examples predicted to belong to each class that are correctly and
    incorrectly classified. These metrics are often called the precision
    (or positive predictive value) and false discovery rate, respectively.
    The row at the bottom of the plot shows the percentages of all the examples
    belonging to each class that are correctly and incorrectly classified.
    These metrics are often called the recall
    (or true positive rate or sensitivity)
    and false negative rate, respectively.
    The cell in the bottom right of the plot shows the overall accuracy.
    """
    # Transpose to set the ground truth to be along columns
    confusion_matrix = insert_totals(confusion_matrix)
    confusion_matrix = confusion_matrix.T
    fig, ax = plt.subplots()
    im = ax.imshow(confusion_matrix, cmap=DEFAULT_CMAP)

    # set ticklabels rotation
    ax.set_xticks(np.arange(confusion_matrix.shape[1]))
    ax.set_yticks(np.arange(confusion_matrix.shape[0]))
    ax.set_xticklabels(labels, rotation=45, fontsize=10)
    ax.set_yticklabels(labels, rotation=0, fontsize=10)

    # Turn off all the ticks
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    # Loop over data dimensions and create text annotations.
    cm_length = confusion_matrix.shape[0]
    for row in range(cm_length):
        for col in range(cm_length):
            text = ax.text(
                row, col,
                format_element_to_matlab_confusion_matrix(
                    row, col, confusion_matrix),
                ha="center", va="center", color="black", fontsize=8)

    # Turn spines off and create black grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.set_xticks(np.arange(confusion_matrix.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(confusion_matrix.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Titles and labels
    ax.set_title('Confusion matrix', fontweight='bold')
    ax.set_xlabel('Output Class', fontweight='bold')
    ax.set_ylabel('Target Class', fontweight='bold')

    # Save figure
    plt.tight_layout()
    plt.savefig(output_fig, dpi=300)
    del im, text


def insert_totals(confusion_matrix):
    """ insert total column and line (the last ones) """
    x = confusion_matrix.shape[0]
    complete_confusion_matrix = np.zeros(
        (x + 1, x + 1), dtype=confusion_matrix.dtype)
    complete_confusion_matrix[0:-1, 0:-1] = confusion_matrix
    complete_confusion_matrix[-1, -1] = np.sum(confusion_matrix)
    for i in range(x):
        complete_confusion_matrix[-1, i] = np.sum(confusion_matrix[:, i])
        complete_confusion_matrix[i, -1] = np.sum(confusion_matrix[i, :])
    return complete_confusion_matrix


home = os.path.expanduser("~")
DATA_DIR = home + "/Downloads/gencode.v32.groundtruth/"
whole_human_genome = DATA_DIR + "gencode.v32.transcripts.fa"
whole_human_genome_protein_translations = \
    DATA_DIR + "gencode.v32.pc_translations.fa"
whole_human_genome_protein_transcripts = \
    DATA_DIR + "gencode.v32.pc_transcripts.fa"
peptide_ksize = 7
molecule = "protein"
coding_peptide_bloom_filter = make_peptide_bloom_filter(
    whole_human_genome_protein_translations,
    peptide_ksize,
    molecule,
    tablesize=1e6)

labels = ["Coding", "Non-coding"]

scoring_lines = []
gt_protein_coding_count = 0
for record in screed.open(whole_human_genome_protein_transcripts):
    gt_protein_coding_count += 1
    line = [record.name, 'Coding']
    scoring_lines.append(line)
gt_df = pd.DataFrame(
    scoring_lines,
    columns=['read_id', 'classification']
)
gt_df.to_csv(DATA_DIR + "gencode.v32.gt_coding_classification.csv")

predicted_df = score_reads(
    whole_human_genome_protein_transcripts,
    coding_peptide_bloom_filter,
    long_reads=True)

predicted_df.to_csv(
    DATA_DIR + "gencode.v32.protein_coding_prediction_classification.csv")

confusion_matrix = sklearn.metrics.confusion_matrix(
    gt_df["classification"].tolist(), predicted_df["classification"].tolist(),
    labels)

np.save(
    DATA_DIR + "gencode.v32.protein_coding_prediction_classification.npy",
    confusion_matrix)

plot_cm(
    confusion_matrix,
    labels,
    DATA_DIR + "confusion_matrix_protein_coding_long_reads.eps")

predicted_df = score_reads(
    whole_human_genome,
    coding_peptide_bloom_filter,
    long_reads=True)

predicted_df.to_csv(
    DATA_DIR + "gencode.v32.whole_genome_prediction_classification.csv")

protein_coding_biotypes = [
    "protein_coding",
    "nonsense_mediated_decay",
    "non_stop_decay",
    "polymorphic_pseudogene"]

protein_coding_biotypes_special_case = ["IG_gene", "TR_gene"]
scoring_lines = []
predicted = []
for index, row in predicted_df.iterrows():
    read_id = row['read_id']
    annotation = read_id.split("|")[-2]
    if annotation in protein_coding_biotypes:
        gt_classification = "Coding"
    elif ((annotation.startswith("IG") or annotation.startswith("TR")) and annotation.endswith("_gene")): 
        gt_classification = "Coding"
    else:
        gt_classification = "Non-coding"
    if row['classification'] == 'Coding':
        predicted.append('Coding')
    else:
        predicted.append('Non-coding')
    scoring_lines.append([read_id, gt_classification])

gt_df = pd.DataFrame(
    scoring_lines,
    columns=['read_id', 'classification']
)
gt_df.to_csv(DATA_DIR + "gencode.v32.gt_whole_genome_classification.csv")

confusion_matrix = sklearn.metrics.confusion_matrix(
    gt_df["classification"].tolist(), predicted,
    labels)

np.save(
    DATA_DIR + "gencode.v32.whole_genome_prediction_classification.npy",
    confusion_matrix)

plot_cm(
    confusion_matrix, labels, DATA_DIR + "confusion_matrix_whole_genome_long_reads.eps")
