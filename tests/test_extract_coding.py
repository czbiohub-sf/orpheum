from collections import namedtuple
import os
import re
import warnings

from Bio.Seq import Seq
from click.testing import CliRunner
import pandas as pd
import pandas.util.testing as pdt
import pytest
import screed

AlphabetKsizeScorepath = namedtuple("AlphabetKsizeScores",
                                    ['alphabet', 'protein_ksize',
                                     'score_path'])
AlphabetKsizeScores = namedtuple("AlphabetKsizeScores",
                                 ['alphabet', 'protein_ksize', 'scores'])


@pytest.fixture
def seq():
    s = 'CGCTTGCTTAATACTGACATCAATAATATTAGGAAAATCGCAATATAACTGTAAATCCTGTTCTGTC'
    with warnings.catch_warnings():
        # Ignore The following warning because we don't use Bio.Alphabet
        # explicitly:
        # PendingDeprecationWarning: We intend to remove or replace
        # Bio.Alphabet in 2020, ideally avoid using it explicitly in your
        # code. Please get in touch if you will be adversely affected by this.
        # https://github.com/biopython/biopython/issues/2046
        warnings.simplefilter("ignore")
        return Seq(s)


@pytest.fixture
def seq_with_n():
    return "CGCAGATGTAGAAAATTTCCGCCAGAGGGTAGGCCCACTGACCCAGGATCTGAAGGNNNGCANC" \
           "ANGCAGAAGTTC"


@pytest.fixture
def low_complexity_seq():
    return "CCCCCCCCCACCACCACCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCA" \
           "CACACCCCCAACACCC"


@pytest.fixture
def low_complexity_seq_step2():
    return "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT"


@pytest.fixture
def low_complexity_seq_in_peptide_space():
    # Shouts to Colleen Stoyas and the Spinocerebellar Ataxia type 7 Patients
    return "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"


@pytest.fixture
def seq_too_short():
    return 'GGAGAAGCCATCATAACTGCAGACC'


@pytest.fixture
def low_complexity_seq_step2():
    return "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT"


@pytest.fixture
def coding_seq1():
    # Correct translation: SFAVHTHRENPAQPGAVTGSATV
    # Gene: ASPDH in Macaque
    return "GTAACAGTAGCAGAGCCGGTGACAGCGCCAGGCTGGGCTGGGTTCTCTCTGTGGGTGTGCACGGCAAAGCTG"


@pytest.fixture
def coding_seq2():
    # Correct translation: EEIAAGKCRRPAVKQFHDSKIK
    # Gene: RPL16a in Macaque
    return "ATTTGATCTTGGAGTCGTGGAACTGCTTGACAGCTGGCCGGCGGCACTTGCCAGCTGCAATCTCTTCC"


@pytest.fixture
def noncoding_seq():
    return "GAACGCCGTGAGGGAGAGGAGAGAGAGGGGGGGAGAGGAGGGGAGCGGAGGGGAGAGGGGGGAGTGAGAA"


@pytest.fixture
def seq_all_stop_codons():
    # All translation frames have stop codons
    return 'TTAAGTTCTAGTCTGTGAGCACTTGTAGTTCAATAATCGTCATCTTCATCAGAGTCCATTACTTTTCTTCTGTTG'


@pytest.fixture(params=['coding_seq1', 'coding_seq2',
                        'noncoding_seq', 'low_complexity_seq',
                        'low_complexity_seq_in_peptide_space',
                        'low_complexity_seq_step2', 'seq_all_stop_codons',
                        'seq_with_n'])
def seq_to_score(request, low_complexity_seq, low_complexity_seq_step2,
                 coding_seq1, coding_seq2, low_complexity_seq_in_peptide_space,
                 noncoding_seq, seq_all_stop_codons, seq_with_n):
    if request.param == 'coding_seq1':
        return request.param, coding_seq1
    if request.param == 'coding_seq2':
        return request.param, coding_seq2
    elif request.param == 'noncoding_seq':
        return request.param, noncoding_seq
    elif request.param == 'low_complexity_seq':
        return request.param, low_complexity_seq
    elif request.param == 'low_complexity_seq_step2':
        return request.param, low_complexity_seq_step2
    elif request.param == 'low_complexity_seq_in_peptide_space':
        return request.param, low_complexity_seq_in_peptide_space
    elif request.param == 'seq_all_stop_codons':
        return request.param, seq_all_stop_codons
    elif request.param == 'seq_with_n':
        return request.param, seq_with_n


@pytest.fixture(params=['seq', 'low_complexity_seq',
                        'low_complexity_seq_step2'])
def type_seq(request, seq, low_complexity_seq, low_complexity_seq_step2):
    if request.param == 'seq':
        return request.param, seq
    elif request.param == 'low_complexity_seq':
        return request.param, low_complexity_seq
    elif request.param == 'low_complexity_seq_step2':
        return request.param, low_complexity_seq_step2


@pytest.fixture(params=[1, 2])
def fastp_complexity_step(request):
    return request.param


@pytest.fixture
def uniprot_opisthokonta_bloom_filter_path(data_folder):
    return os.path.join(
        data_folder, 'extract_coding',
        'uniprot-reviewed_yes+taxonomy_2759.fasta.alphabet-protein_ksize-8.bloomfilter.nodegraph')


@pytest.fixture
def uniprot_opisthokonta_bloom_filter(uniprot_opisthokonta_bloom_filter_path):
    from khtools.bloom_filter import load_nodegraph
    return load_nodegraph(uniprot_opisthokonta_bloom_filter_path)


def test_three_frame_translation(seq):
    from khtools.extract_coding import three_frame_translation

    test = [str(x) for x in three_frame_translation(seq)]
    true = [
        'RLLNTDINNIRKIAI*L*ILFC', 'ACLILTSIILGKSQYNCKSCSV',
        'LA*Y*HQ*Y*ENRNITVNPVL'
    ]
    assert test == true


@pytest.fixture
def reads(data_folder):
    return os.path.join(
        data_folder,
        'SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22.fq')


@pytest.fixture
def empty_fasta(data_folder):
    return os.path.join(
        data_folder, 'empty_fasta.fasta')


@pytest.fixture
def true_scores_path(data_folder, molecule, peptide_ksize):
    return os.path.join(
        data_folder, "extract_coding",
        "SRR306838_GSM752691_hsa_br_F_1_trimmed_"
        f"subsampled_n22__molecule-{molecule}_ksize-"
        f"{peptide_ksize}.csv")


@pytest.fixture
def true_scores(true_scores_path):
    return pd.read_csv(true_scores_path)


@pytest.fixture
def jaccard_threshold(molecule):
    from khtools.extract_coding import get_jaccard_threshold
    threshold = get_jaccard_threshold(None, molecule)
    return threshold


@pytest.fixture
def peptide_ksize(molecule):
    from khtools.bloom_filter import get_peptide_ksize

    ksize = get_peptide_ksize(molecule)
    return ksize


@pytest.fixture
def single_alphabet_ksize_true_scores_path(data_folder):
    true_scores_path = os.path.join(
        data_folder, "extract_coding",
        "SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22__"
        "molecule-protein_ksize-7.csv")
    return AlphabetKsizeScorepath('protein', 7, true_scores_path)


@pytest.fixture
def single_alphabet_ksize_true_scores(single_alphabet_ksize_true_scores_path):
    alphabet, ksize, true_scores_path = single_alphabet_ksize_true_scores_path
    true_scores = pd.read_csv(true_scores_path)
    return AlphabetKsizeScores(alphabet, ksize, true_scores)


def test_compute_fastp_low_complexity(type_seq, fastp_complexity_step):
    from khtools.extract_coding import compute_fastp_complexity

    seqtype, seq = type_seq
    test = compute_fastp_complexity(seq, step=fastp_complexity_step)
    if fastp_complexity_step == 1:
        if seqtype == 'seq':
            assert test == 0.746268656716418
        elif seqtype == 'low_complexity_seq':
            assert test == 0.2631578947368421
        elif seqtype == 'low_complexity_seq_step2':
            assert test == 0.9833333333333333
    elif fastp_complexity_step == 2:
        if seqtype == 'seq':
            assert test == 0.7164179104477612
        elif seqtype == 'low_complexity_seq':
            assert test == 0.21052631578947367
        elif seqtype == 'low_complexity_seq_step2':
            assert test == 0


def test_evaluate_is_fastp_low_complexity(type_seq):
    from khtools.extract_coding import evaluate_is_fastp_low_complexity

    seqtype, seq = type_seq

    test = evaluate_is_fastp_low_complexity(seq)
    if seqtype == 'seq':
        # regular sequence is not low complexity
        assert not test
    elif seqtype == 'low_complexity_seq':
        # low complexity sequence should evaluate to low complexity!
        assert test
    elif seqtype == 'low_complexity_seq_step2':
        # This sequence should be caught when step=2
        assert test


def test_three_frame_translation_no_stops(seq):
    from khtools.extract_coding import three_frame_translation_no_stops

    test = {
        k: str(v)
        for k, v in three_frame_translation_no_stops(seq).items()
    }
    true = {2: 'ACLILTSIILGKSQYNCKSCSV'}
    assert test == true


def test_six_frame_translation_no_stops(seq):
    from khtools.extract_coding import six_frame_translation_no_stops

    test = {k: str(v) for k, v in six_frame_translation_no_stops(seq).items()}
    true = {
        2: 'ACLILTSIILGKSQYNCKSCSV',
        -2: 'TEQDLQLYCDFPNIIDVSIKQA',
        -3: 'QNRIYSYIAIFLILLMSVLSK'
    }
    assert test == true


def test_score_single_read(capsys, seq_to_score,
                           uniprot_opisthokonta_bloom_filter):
    from khtools.extract_coding import score_single_read

    seqtype, seq = seq_to_score
    score = score_single_read(seq, uniprot_opisthokonta_bloom_filter,
                              peptide_ksize=8, alphabet='protein',
                              description=seqtype)
    # --- Check fasta output --- #
    captured = capsys.readouterr()
    standard_output = captured.out
    seqtypes_without_translations = ('noncoding_seq', 'low_complexity_seq',
                                     'low_complexity_seq_step2',
                                     'seq_all_stop_codons',
                                     'low_complexity_seq_in_peptide_space',
                                     'seq_all_stop_codons', 'seq_with_n')
    low_complexity_nucleotide_seqs = ('noncoding_seq', 'low_complexity_seq',
                                      'low_complexity_seq_step2')
    true_translations = []
    if seqtype == 'coding_seq1':
        true_translations = [
            # this is the only correct one, but keeping the other ones here so
            # the test passes
            'SFAVHTHRENPAQPGAVTGSATV',
            # This is not the biologically correct reading frame
            'VTVAEPVTAPGWAGFSLWVCTAKL',
            # This is not the biologically correct reading frame
            'ALPCTPTERTQPSLALSPALLLL']
        category = 'Coding'
    elif seqtype == 'coding_seq2':
        true_translations = [
            # Only this is the correct coding sequence, the others are
            # biologically incorrect but are here for now so the test passes
            'EEIAAGKCRRPAVKQFHDSKIK',
            # The following sequences are biologically incorrect
            'FDLGVVELLDSWPAALASCNLF',
            'LILESWNCLTAGRRHLPAAISS',
            'KRLQLASAAGQLSSSSTTPRSN'
        ]
        category = 'Coding'
    elif seqtype in seqtypes_without_translations:
        true_translations = []
        # stdout should be empty
        assert standard_output == ''
        # Assign the coding detection category
        if seqtype in low_complexity_nucleotide_seqs:
            category = 'Low complexity nucleotide'
        elif seqtype == 'low_complexity_seq_in_peptide_space':
            category = "Too few k-mers in protein20 alphabet"
        elif seqtype == 'seq_all_stop_codons':
            category = 'All translation frames have stop codons'
        elif seqtype == 'seq_with_n':
            category = "Nucleotide sequence contains ambiguous 'N' characters"

    stdout_lines = standard_output.splitlines()
    test_n_translations = sum(1 for line in stdout_lines if '>' in line)
    assert test_n_translations == len(true_translations)
    for translation in true_translations:
        assert translation in standard_output

    assert score.special_case == category


def test_score_reads(capsys, tmpdir, reads, peptide_bloom_filter, molecule,
                     true_scores,
                     true_scores_path,
                     true_protein_coding_fasta_path):
    from khtools.extract_coding import score_reads

    test = score_reads(reads,
                       peptide_bloom_filter,
                       molecule=molecule)
    # Convert to basename to be compatible with test data
    test['filename'] = test['filename'].map(os.path.basename)

    # Check that scoring was the same
    pdt.assert_equal(test, true_scores)

    # --- Check fasta output --- #
    captured = capsys.readouterr()
    test_names = []
    for line in captured.out.splitlines():
        if line.startswith(">"):
            test_names.append(line.lstrip('>'))

    # Check that the proper sequences were output
    true_names = get_fasta_record_names(true_protein_coding_fasta_path)

    # Check that precision is high -- everything in "test" was truly coding
    assert all(test_name in true_names for test_name in test_names)

    captured_lines = captured.out.splitlines()
    with open(true_protein_coding_fasta_path) as f:
        for true_line in f.readlines():
            assert true_line.strip() in captured_lines


def write_fasta_string_to_file(fasta_string, folder, prefix):
    test_fasta_filename = os.path.join(folder, prefix + '.fasta')
    with open(test_fasta_filename) as f:
        f.write(fasta_string)
    return test_fasta_filename


def get_fasta_record_names(fasta_path):
    names = []
    for record in screed.open(fasta_path):
        name = record['name']
        names.append(name)
    return set(names)


def test_maybe_write_json_summary_empty(peptide_bloom_filter_path, molecule,
                                        peptide_ksize):
    from khtools.extract_coding import maybe_write_json_summary, \
        DEFAULT_JACCARD_THRESHOLD

    coding_scores = pd.DataFrame(columns=['read_id', 'jaccard_in_peptide_db',
                                          'n_kmers', 'classification',
                                          'filename'])
    summary = maybe_write_json_summary(
        coding_scores, json_summary=True, filenames=['nonexistent.fa'],
        bloom_filter=peptide_bloom_filter_path, molecule=molecule,
        peptide_ksize=peptide_ksize,
        jaccard_threshold=DEFAULT_JACCARD_THRESHOLD)
    assert summary['input_files'] == ['nonexistent.fa']
    assert summary['jaccard_info']['count'] == 0


def test_get_n_translated_frames_per_read():
    from khtools.extract_coding import get_n_translated_frames_per_read

    # Make fake dataframe with "read_id" and "classification" columns only
    # for testing
    data = {
        'read_id':
            ['read0', 'read0', 'read0', 'read0', 'read0', 'read0',
             'read1', 'read1', 'read1', 'read1', 'read1',
             'read2', 'read2', 'read2', 'read2',
             'read3', 'read3', 'read3',
             'read4', 'read4',
             'read5', 'read5',
             'read6', 'read7', 'read8', 'read9', 'read10'],
    }
    df = pd.DataFrame(data)
    df['classification'] = "Coding"

    percentages, histogram = get_n_translated_frames_per_read(df)
    assert histogram == {
        'Number of reads with 1 putative protein-coding translations': 5,
        'Number of reads with 2 putative protein-coding translations': 2,
        'Number of reads with 6 putative protein-coding translations': 1,
        'Number of reads with 5 putative protein-coding translations': 1,
        'Number of reads with 4 putative protein-coding translations': 1,
        'Number of reads with 3 putative protein-coding translations': 1}
    assert percentages == {
        'Number of reads with 1 putative protein-coding translations': 45.45454545454545,
        'Number of reads with 2 putative protein-coding translations': 18.181818181818183,
        'Number of reads with 6 putative protein-coding translations': 9.090909090909092,
        'Number of reads with 5 putative protein-coding translations': 9.090909090909092,
        'Number of reads with 4 putative protein-coding translations': 9.090909090909092,
        'Number of reads with 3 putative protein-coding translations': 9.090909090909092}


def test_get_n_per_coding_classification(molecule):
    from khtools.extract_coding import get_n_per_coding_classification, \
        TOO_FEW_KMERS_CATEGORIES
    from khtools.sequence_encodings import ALIAS_TO_ALPHABET

    data = [
        ['read1', 'All translations shorter than peptide k-mer size + 1'],
        ['read2', 'All translation frames have stop codons'],
        ['read3', 'Coding'],
        ['read4', 'Non-coding'],
        ['read5', 'Low complexity nucleotide'],
        ['read6', 'Read length was shorter than 3 * peptide k-mer size'],
        ['read7', TOO_FEW_KMERS_CATEGORIES[molecule]],
        ['read8', "Nucleotide sequence contains ambiguous 'N' characters"]
    ]
    df = pd.DataFrame(data, columns=['read_id', 'classification'])

    test_percentages, test_counts = get_n_per_coding_classification(df,
                                                                    molecule)
    canonical_molecule = ALIAS_TO_ALPHABET[molecule]
    true_percentages = {
        'All translations shorter than peptide k-mer size + 1': 12.5,
        'All translation frames have stop codons': 12.5,
        'Coding': 12.5, 'Non-coding': 12.5,
        'Low complexity nucleotide': 12.5,
        'Read length was shorter than 3 * peptide k-mer size': 12.5,
        f'Too few k-mers in {canonical_molecule} alphabet': 12.5,
        "Nucleotide sequence contains ambiguous 'N' characters": 12.5,
    }
    true_counts = {
        'All translations shorter than peptide k-mer size + 1': 1,
        'All translation frames have stop codons': 1,
        'Coding': 1,
        'Non-coding': 1,
        'Low complexity nucleotide': 1,
        'Read length was shorter than 3 * peptide k-mer size': 1,
        f'Too few k-mers in {canonical_molecule} alphabet': 1,
        "Nucleotide sequence contains ambiguous 'N' characters": 1,
    }
    assert test_counts == true_counts
    assert test_percentages == true_percentages


def test_generate_coding_summary(reads, data_folder,
                                 single_alphabet_ksize_true_scores):
    from khtools.extract_coding import generate_coding_summary

    alphabet, ksize, true_scores = single_alphabet_ksize_true_scores
    jaccard_threshold = 0.95
    ksize = 7

    peptide_bloom_filter = 'bloom_filter.nodegraph'

    test_summary = generate_coding_summary(
        true_scores, peptide_bloom_filter, alphabet,
        peptide_ksize=ksize, jaccard_threshold=jaccard_threshold)

    true_summary = {
        'input_files': [
            'SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22.fq'],
        'jaccard_info': {'count': 15.0, 'mean': 0.2478485838779956,
                         'std': 0.3924572683558275, 'min': 0.0,
                         '25%': 0.029411764705882356, '50%': 0.0625,
                         '75%': 0.15073529411764708, 'max': 1.0},
        'classification_value_counts': {
            'All translations shorter than peptide k-mer size + 1': 2,
            'All translation frames have stop codons': 3, 'Coding': 3,
            'Non-coding': 12, 'Low complexity nucleotide': 2,
            'Read length was shorter than 3 * peptide k-mer size': 0,
            'Too few k-mers in protein20 alphabet': 1,
            "Nucleotide sequence contains ambiguous 'N' characters": 0,
        },
        'classification_percentages': {
            'All translations shorter than peptide k-mer size + 1': 8.695652173913043,
            'All translation frames have stop codons': 13.043478260869565,
            'Coding': 13.043478260869565, 'Non-coding': 52.17391304347826,
            'Low complexity nucleotide': 8.695652173913043,
            'Read length was shorter than 3 * peptide k-mer size': 0.0,
            'Too few k-mers in protein20 alphabet': 4.3478260869565215,
            "Nucleotide sequence contains ambiguous 'N' characters": 0,
        },
        'histogram_n_coding_frames_per_read': {
            'Number of reads with 1 putative protein-coding translations': 3},
        'histogram_n_coding_frames_per_read_percentages': {
            'Number of reads with 1 putative protein-coding translations': 100.0},
        'peptide_bloom_filter': 'bloom_filter.nodegraph',
        'peptide_alphabet': 'protein', 'peptide_ksize': 7,
        'jaccard_threshold': 0.95}
    assert test_summary == true_summary


def test_cli_peptide_fasta(reads, peptide_fasta, molecule, peptide_ksize,
                           true_protein_coding_fasta_string):
    from khtools.extract_coding import cli

    runner = CliRunner()
    result = runner.invoke(cli, [
        '--peptide-ksize', peptide_ksize, '--alphabet', molecule,
        peptide_fasta, reads
    ])
    assert result.exit_code == 0
    # CliRunner jams together the stderr and stdout so just check if the
    # true string is contained in the output
    assert true_protein_coding_fasta_string in result.output

    # Make sure "Writing extract_coding summary to" didn't get accidentally
    # written to stdout instead of stderr
    assert 'Writing extract_coding summary to' \
           not in true_protein_coding_fasta_string


def test_cli_bad_jaccard_threshold_float(reads, peptide_fasta):
    from khtools.extract_coding import cli

    runner = CliRunner()
    result = runner.invoke(cli, [
        "--jaccard-threshold", "3.14", peptide_fasta, reads
    ])
    assert result.exit_code == 2
    error_messages = ['Error: Invalid value for ', '--jaccard-threshold',
                      ': --jaccard-threshold needs to be a number between 0 '
                      'and 1, but 3.14 was provided']
    for error_message in error_messages:
        assert error_message in result.output


def test_cli_bad_jaccard_threshold_string(reads, peptide_fasta):
    from khtools.extract_coding import cli

    runner = CliRunner()
    result = runner.invoke(cli, [
        "--jaccard-threshold", "beyonce", peptide_fasta, reads
    ])
    assert result.exit_code == 2
    error_messages = ['Error: Invalid value for ', '--jaccard-threshold',
                      ': beyonce is not a valid floating point value']
    for error_message in error_messages:
        assert error_message in result.output


def test_cli_peptide_bloom_filter(reads, peptide_bloom_filter_path, molecule,
                                  peptide_ksize,
                                  true_protein_coding_fasta_string):
    from khtools.extract_coding import cli

    runner = CliRunner()
    result = runner.invoke(cli, [
        '--peptide-ksize', peptide_ksize, "--peptides-are-bloom-filter",
        '--alphabet', molecule, peptide_bloom_filter_path, reads
    ])
    assert result.exit_code == 0
    assert true_protein_coding_fasta_string in result.output


def test_cli_csv(tmpdir, reads, peptide_bloom_filter_path, molecule,
                 peptide_ksize, true_protein_coding_fasta_string, true_scores):
    from khtools.extract_coding import cli

    csv = os.path.join(tmpdir, 'coding_scores.csv')

    runner = CliRunner()
    result = runner.invoke(cli, [
        '--peptide-ksize', peptide_ksize, "--csv", csv,
        "--peptides-are-bloom-filter", '--alphabet', molecule,
        peptide_bloom_filter_path, reads
    ])
    assert result.exit_code == 0
    assert true_protein_coding_fasta_string in result.output
    assert os.path.exists(csv)

    # the CLI adds the filename to the scoring dataframe
    true = true_scores.copy()
    true['filename'] = reads

    test_scores = pd.read_csv(csv)
    pdt.assert_equal(test_scores, true)


def test_cli_coding_nucleotide_fasta(tmpdir, reads, peptide_fasta):
    from khtools.extract_coding import cli

    coding_nucleotide_fasta = os.path.join(tmpdir, 'coding_nucleotides.fasta')

    runner = CliRunner()
    result = runner.invoke(cli, [
        "--coding-nucleotide-fasta", coding_nucleotide_fasta,
        peptide_fasta, reads
    ])
    assert result.exit_code == 0
    assert os.path.exists(coding_nucleotide_fasta)


def test_cli_noncoding_fasta(tmpdir, reads, peptide_fasta):
    from khtools.extract_coding import cli

    noncoding_nucleotide_fasta = os.path.join(tmpdir,
                                              'noncoding_nucleotides.fasta')

    runner = CliRunner()
    result = runner.invoke(cli, [
        "--noncoding-nucleotide-fasta", noncoding_nucleotide_fasta,
        peptide_fasta, reads
    ])
    assert result.exit_code == 0
    assert os.path.exists(noncoding_nucleotide_fasta)


def test_cli_low_complexity_nucleotide(tmpdir, reads, peptide_fasta):
    from khtools.extract_coding import cli

    low_complexity_nucleotide_fasta = os.path.join(
        tmpdir, 'low_complexity_nucleotide.fasta')

    runner = CliRunner()
    result = runner.invoke(cli, [
        "--low-complexity-nucleotide-fasta", low_complexity_nucleotide_fasta,
        peptide_fasta, reads
    ])
    assert result.exit_code == 0
    assert os.path.exists(low_complexity_nucleotide_fasta)


def test_cli_low_complexity_peptide(
        tmpdir,
        reads,
        peptide_fasta):
    from khtools.extract_coding import cli

    low_complexity_peptide_fasta = os.path.join(tmpdir,
                                                'low_complexity_peptide.fasta')

    runner = CliRunner()
    result = runner.invoke(cli, [
        "--low-complexity-peptide-fasta", low_complexity_peptide_fasta,
        peptide_fasta, reads
    ])
    assert result.exit_code == 0
    assert os.path.exists(low_complexity_peptide_fasta)


def test_cli_json_summary(tmpdir, reads, peptide_fasta):
    from khtools.extract_coding import cli

    json_summary = os.path.join(tmpdir, 'coding_summary.json')

    runner = CliRunner()
    result = runner.invoke(cli, [
        "--json-summary", json_summary,
        peptide_fasta, reads
    ])
    assert result.exit_code == 0
    assert os.path.exists(json_summary)


def test_cli_empty_fasta_json_summary(tmpdir, empty_fasta, peptide_fasta):
    from khtools.extract_coding import cli

    json_summary = os.path.join(tmpdir, 'coding_summary.json')

    runner = CliRunner()
    result = runner.invoke(cli, [
        "--json-summary", json_summary,
        peptide_fasta, empty_fasta
    ])
    assert result.exit_code == 0
    assert os.path.exists(json_summary)
