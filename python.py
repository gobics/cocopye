"""
This is a playground for testing new functions. It should not be part of the final application/library.
"""

from src.libcc.preprocessing import kmer
from src.libcc.preprocessing import prodigal
from src.libcc.preprocessing.pfam import uproc, count_pfams
import typing as tp
import time


def test_sequence():
    st = time.time()
    sequences: tp.List[kmer.Sequence] = kmer.read_fasta_file("../concatorfs_prodigal.fasta", kmer.PROTEIN)
    end = time.time()
    print(end-st)
    st = time.time()
    counts = [seq.kmer_count(4) for seq in sequences]
    end = time.time()
    print(end-st)


def test_prodigal():
    prodigal("/home/nemo/exclude_from_backup/prodigal/prodigal.linux",
             "/home/nemo/exclude_from_backup/prodigal/subset.fasta",
             "/home/nemo/exclude_from_backup/prodigal/output.fasta",
             True)


def test_uproc():
    uproc_result = uproc("/home/nemo/exclude_from_backup/uproc/uproc-prot",
          "/home/nemo/exclude_from_backup/uproc/pfam_db",
          "/home/nemo/exclude_from_backup/uproc/model",
          "/home/nemo/exclude_from_backup/prodigal/output.fasta",
        False)

    print(uproc_result.shape)


def test_uproc_orf():
    count_mat, sequences = count_pfams("/home/nemo/.local/share/libcc/uproc/bin/uproc-orf", "/home/nemo/exclude_from_backup/uproc/uproc-prot",
          "/home/nemo/exclude_from_backup/uproc/pfam_db",
          "/home/nemo/exclude_from_backup/uproc/model","/run/media/nemo/74691D6C3C589D1D/set_1/fasta")

    print(count_mat.mat().shape)
    print(sequences)
    count_mat.save_to_file("mat.npy")


test_uproc_orf()
