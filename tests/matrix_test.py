from src.cocopye.preprocessing import prodigal, kmer
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
