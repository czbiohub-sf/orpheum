import pysam
import screed

from orpheum.log_utils import get_logger

logger = get_logger(__file__)


class ReadParser:
    """An abstraction of bam and fasta file format to simply return reads: ids and sequences"""

    def __init__(self, filename, soft_clipping=False):
        """

        :param filename:
        :param soft_clipping: bool
            Applies to bam files only. If true, then return only the aligned segment of the read.
            Doesn't return the bases that were "soft clipped"
        """

        self.filename = filename
        self.soft_clipping = soft_clipping

        try:
            # Peek at the file with screed and see if it's a valid fasta/fastq format
            # Use "with" so it closes automatically
            with screed.open(self.filename) as _:
                pass
            self.filetype = "fastx"
            if self.soft_clipping:
                logger.warn(
                    "soft_clipping set to 'True' on a fastq/fasta file, but "
                    "this has no effect on these files"
                )
        except ValueError:
            try:
                pysam.AlignmentFile(self.filename)
                self.filetype = "bam"
            except ValueError:
                raise ValueError(
                    f"Input file {self.filename} was not fastq/fasta or bam"
                )

    def __iter__(self):
        if self.filetype == "fastx":
            return self._iter_reads_fastx()
        elif self.filetype == "bam":
            return self._iter_reads_bam()

    def _iter_reads_fastx(self):
        with screed.open(self.filename) as records:
            for record in records:
                yield record["name"], record["sequence"]

    def _iter_reads_bam(self):
        bam = pysam.AlignmentFile(self.filename)
        for record in bam:
            yield record.qname, record.query_sequence
