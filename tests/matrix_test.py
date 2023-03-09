from src.libcc.preprocessing import sequence
import typing as tp
import time


def test_sequence():
    st = time.time()
    sequences: tp.List[sequence.Sequence] = sequence.read_fasta_file("../concatorfs_prodigal.fasta", sequence.PROTEIN)
    end = time.time()
    print(end-st)
    st = time.time()
    counts = [seq.kmer_count(4) for seq in sequences]
    end = time.time()
    print(end-st)