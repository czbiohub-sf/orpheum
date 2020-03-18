from collections import Counter
from tqdm import tqdm

from .sequence_encodings import encode_peptide
from .compare_kmer_content import kmerize

ENCODINGS_TO_COUNT = 'hydrophobic-polar', 'dayhoff', 'protein'
# U = Selenocystine: https://en.wikipedia.org/wiki/Selenocysteine
# X = Unknown
SELENOCYSTEINE = "U"
UNKNOWN = "X"
UNWANTED_AMINO_ACIDS = SELENOCYSTEINE, UNKNOWN
STOP_CODON = "*"

def remove_unwanted_kmers(kmers, unwanted):
    return {kmer for kmer in kmers if all(x not in kmer for x in unwanted)}


def count_kmers_single_alphabet_ksize(filename, ksize, alphabet,
                                      verbose=False):
    all_kmers = Counter()

    n_seqs_with_unknown = 0
    n_seqs_with_selenocysteine = 0

    with screed.open(filename) as records:
        for record in tqdm(records):
            seq = record['sequence']
            # "*" = stop codon, and we don't want no stop codons
            if STOP_CODON in seq:
                continue

            if UNKNOWN in seq:
                n_seqs_with_unknown += 1
            if SELENOCYSTEINE in seq:
                n_seqs_with_selenocysteine += 1

            encoded = encode_peptide(seq, alphabet)

            if verbose:
                print(seq)
                print(encoded)
            try:
                kmers = kmerize(encoded, ksize)
                kmers = remove_unwanted_kmers(kmers, unwanted_amino_acids)

                if verbose:
                    print(kmers)
                all_kmers.update(kmers)
            except ValueError:
                # Typically, sequence is too short to k-merize
                continue

    if verbose:
        print(f'n_seqs_with_x: {n_seqs_with_x}')
        print(f'n_seqs_with_u: {n_seqs_with_u}')
    return all_kmers


def count_kmers_multiple_ksizes_encodings(
        filename,
        encodings=ENCODINGS_TO_COUNT,
        verbose=False,
        unwanted_amino_acids=UNWANTED_AMINO_ACIDS,
        pickle_file_prefix=None):
    encoding_ksize_n_kmers = {}

    for encoding in encodings:
        if verbose:
            print(f'--- encoding: {encoding} ---')
        for ksize in ksizes:
            if verbose:
                print(f'\t--- ksize: {ksize} ---')
            all_kmers = count_kmers_single_encoding_ksize(filename, en)
            n_kmers = len(all_kmers)
            if verbose:
                print(f'\t\tn_kmers: {n_kmers}')
            encoding_ksize_n_kmers[(encoding, ksize)] = n_kmers
            maybe_pickle_kmers(pickle_file_prefix, all_kmers, encoding, ksize)
            del all_kmers

    return encoding_ksize_n_kmers


def maybe_pickle_kmers(pickle_file_prefix, all_kmers, encoding, ksize):
    if pickle_file_prefix is not None:
        pickle_file = f'{pickle_file_prefix}__molecule-{encoding}_ksize-{ksize}.pickle'
        with open(filename, 'wb') as f:
            pickle.dump(all_kmers, f)