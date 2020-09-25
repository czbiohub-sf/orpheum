import os

import screed

import sencha.constants_translate as constants_translate
from sencha.log_utils import get_logger


logger = get_logger(__file__)


def calculate_chunksize(total_jobs_todo, processes):
    """
    Return integer - chunksize representing the number of jobs
    per process that needs to be run
    total_jobs_todo : int
        total number of jobs
    processes; int
        number of processes to be used for multiprocessing
    Returns
    -------
    Integer reprsenting number of jobs to be run on each process
    """
    chunksize, extra = divmod(total_jobs_todo, processes)
    if extra:
        chunksize += 1
    return chunksize


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def split_fasta_files(reads, chunksize, outdir):
    if not os.path.exists(os.path.abspath(outdir)):
        os.makedirs(outdir)
    basename = os.path.basename(reads).split(".")[0]
    record_iter = iter(screed.open(reads))
    filenames = []
    for i, batch in enumerate(batch_iterator(record_iter, chunksize)):
        filename = os.path.join(outdir, "{}_{}.fasta".format(basename, i + 1))
        with open(filename, "w") as ouput_handle:
            for record in batch:
                description = record["name"]
                sequence = record["sequence"]
                write_fasta(ouput_handle, description, sequence)
        filenames.append(filename)
    return filenames


def write_fasta(file_handle, description, sequence):
    file_handle.write(">{}\n{}\n".format(description, sequence))


def maybe_write_fasta(file_handle, description, sequence):
    """Write fasta to file handle if it is not None"""
    if file_handle is not None:
        write_fasta(file_handle, description, sequence)


def open_and_announce(
    filename, key, announcement_dict=constants_translate.SEQTYPE_TO_ANNOUNCEMENT
):
    """Return an opened file handle to write and announce"""
    announcement = announcement_dict[key]
    logger.info("Writing {} to {}".format(announcement, filename))
    return open(filename, "w")


def maybe_open_fastas(
    fastas, announcement_dict=constants_translate.SEQTYPE_TO_ANNOUNCEMENT
):
    file_handles = {}
    for key, fasta in fastas.items():
        if fasta is not None:
            file_handles[key] = open_and_announce(fasta, key, announcement_dict)
        else:
            file_handles[key] = None
    return file_handles


def maybe_close_fastas(file_handles):
    for file_handle in file_handles.values():
        if file_handle is not None:
            file_handle.close()
