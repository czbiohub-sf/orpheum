from tqdm import tqdm

from .sequence_encodings import encode_peptide
from .compare_kmer_content import kmerize

ENCODINGS_TO_COUNT = 'hydrophobic-polar', 'dayhoff', 'protein'
# U = Selenocystine: https://en.wikipedia.org/wiki/Selenocysteine
# X = Unknown
UNWANTED_AMINO_ACIDS = "U", "X"


def remove_unwanted_kmers(kmers, unwanted):
    return {kmer for kmer in kmers if all(x not in kmer for x in unwanted)}


def count_kmers_single_alphabet_ksize(filename, ksize, alphabet,
                                      verbose=False):
    all_kmers = set()

    n_seqs_with_x = 0
    n_seqs_with_u = 0

    with screed.open(filename) as records:
        for record in tqdm(records):
            seq = record['sequence']
            # No stop codons
            if '*' in seq:
                continue

            if "X" in seq:
                n_seqs_with_x += 1
            if "U" in seq:
                n_seqs_with_u += 1

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
        unwanted_amino_acids=UNWANTED_AMINO_ACIDS):
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

            del all_kmers

    return encoding_ksize_n_kmers


def maybe_pickle_kmers(encoding, ksize):
    pickle_file = f'{pickle_folder}/Homo_sapiens.GRCh38.pep.all__kmers__molecule-{encoding}_ksize-{ksize}.pickle'
    with open(filename, 'wb') as f:
        pickle.dump(all_kmers, f)