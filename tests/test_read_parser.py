import os

import pytest
from orpheum.read_parser import ReadParser


@pytest.fixture()
def bamfile(data_folder):
    basename = 'possorted_genome_bam__rbfox3_exon13_chr11_118494068-118494243.bam'
    return os.path.join(data_folder, basename)


def test_fasta(reads):
    """Super simple test to make sure fastq/fasta parsing worked"""
    readparser = ReadParser(reads)
    for read_name, read_seq in readparser:
        assert read_name == 'SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1'
        assert read_seq == 'CGCTTGCTTAATACTGACATCAATAATATTAGGAAAATCGCAATATAACTGTAAATCCTGTTCTGTC'
        break

def test_bam(bamfile):
    """Super simple test to make sure bam/sam parsing worked"""
    readparser = ReadParser(bamfile)
    for read_name, read_seq in readparser:
        assert read_name == 'NB501938:262:HG27TBGXK:3:12408:19165:6261'
        assert read_seq == 'GCATTTGAGTGGCTGTTGCTCAGCTGCCATCTGGCACTTTGTCCCAAAGGTTCGTTTGCTCGCTCGGTTTCCTCTCATCCCCATGTACTC'
        break
