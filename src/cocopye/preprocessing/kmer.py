"""
This module contains several classes and functions for reading input sequences and transforming them into a count
matrix.

Example usage:
```python
from cocopye.preprocessing.sequence import Sequence, Alphabet, DNA

# Create a new sequence using a given alphabet
seq1 = Sequence("AACCGGTT", alphabet=DNA)
# Create a new sequence with its own alphabet
seq2 = Sequence("ACGT$TTGA", alphabet=Alphabet("ACGT", ignored="$"))

# Count all 4-mers of the created sequences
kmer_count1 = seq1.kmer_count(4)
# Here, due to the alphabet of seq2, k-mers that contain '$' aren't counted.
kmer_count2 = seq2.kmer_count(4)

# Read protein sequences from a file ...
sequences = read_fasta_file("file.fasta", PROTEIN)
# ... and count their 4-mers
count_matrix = sequences_to_count_matrix(sequences, 4)
```
"""

from __future__ import annotations

import os
import sys

import numpy as np
import numpy.typing as npt
import typing as tp
from typing import Tuple, List
from Bio import SeqIO
from numba import njit, types
from numba.typed import Dict
from tqdm import tqdm

from ..matrices import QueryMatrix, DatabaseMatrix


class Sequence:
    """
    This class represents a sequence. Please note that it is mainly intended as a collection of preprocessing functions
    for our tool and not suitable as a general purpose sequence class.
    """
    metadata: str
    """This can be used to store additional information such as FASTA headers."""
    _alphabet: Alphabet
    _seq: npt.NDArray[np.int8]

    def __init__(self, seq: str, alphabet: Alphabet, metadata: str = ""):
        """
        :param seq: Sequence as an ASCII-encoded string. (As long as you only use letters and maybe some common special
        characters such as `*` the encoding shouldn't be a limitation.)
        :param alphabet: An `Alphabet` object. You can use `DNA`/`PROTEIN` or create your own.
        :param metadata: Additional information. Currently this is not in use.
        """
        self.metadata = metadata
        self._alphabet = alphabet
        self._seq = alphabet.translate_sequence(seq)

    def kmer_count(self, k: int) -> npt.NDArray[np.uint8]:
        """
        Count all k-mers. If a k-mer contains symbols that are ignored by the sequence alphabet, it is not counted.
        :param k: k-mer length
        :return: A 1-dimensional numpy array of length n^k, where n is the size of the sequence alphabet. Each element
        represents the number of occurences of one k-mer. The k-mer itself is encoded in the index, which means that it
        is guaranteed that in sequences that are based on the same alphabet the same k-mer is always in the same
        position.
        """
        return _numba_kmer_count(self._alphabet.len(), k, self._seq)


# @njit  # type: ignore
def _numba_kmer_count(alphabet_size: int, k: int, seq: npt.NDArray[np.int8]) -> npt.NDArray[np.uint8]:  # TODO: Stop at max_value
    arr = np.zeros(alphabet_size ** k, dtype=np.uint8)

    for idx in range(len(seq)):
        kmer_idx = _numba_kmer_idx(seq[idx:idx+k], alphabet_size)
        if kmer_idx == -1:
            continue
        if arr[kmer_idx] == 255:
            continue
        arr[kmer_idx] += 1

    return arr


_numba_kmer_count = njit(_numba_kmer_count)


@njit  # type: ignore
def _numba_kmer_idx(kmer: npt.NDArray[np.int8], alphabet_size: int) -> int:
    kmer_idx = 0
    for c in kmer:
        if c == -1:
            return -1  # return None didn't work with numba, so we use this as a workaround

        kmer_idx *= alphabet_size
        kmer_idx += c

    return kmer_idx


class Alphabet:
    symbols: str
    ignored: str
    _map: Dict

    def __init__(self, symbols: str, ignored: str):
        self.symbols = symbols
        self.ignored = ignored

        self._map = Dict.empty(
            key_type=types.byte,
            value_type=types.int8
        )
        for i, c in enumerate(bytes(symbols, "ASCII")):
            self._map[types.byte(c)] = types.int8(i)
        for c in bytes(ignored, "ASCII"):
            self._map[types.byte(c)] = types.int8(-1)

    def translate_sequence(self, seq: str) -> npt.NDArray[types.int8]:
        return _numba_translate_sequence(self._map, bytes(seq, "ASCII"))

    def len(self) -> int:
        return len(self.symbols)


# @njit  # type: ignore
def _numba_translate_sequence(dictx: Dict, seq: bytes) -> npt.NDArray[np.int8]:
    arr = np.empty(len(seq), dtype=types.int8)
    for idx, c in enumerate(seq):
        if c not in dictx:
            print(c)
        arr[idx] = dictx[c]
    return arr


_numba_translate_sequence = njit(_numba_translate_sequence)  # TODO: Check if this is as fast as the decorator


DNA: Alphabet = Alphabet("ACGT", "")
"""`Alphabet("ACGT", "")`"""
PROTEIN: Alphabet = Alphabet("ARNDCEQGHILKMFPSTWYV", "XB$*")
"""`Alphabet("ARNDCEQGHILKMFPSTWYV", "XB$*")`"""


def create_database_matrix(
        fasta_file: str,
        alphabet: Alphabet = PROTEIN,
        sequences: tp.Optional[tp.List[str]] = None,
        k: int = 4
) -> DatabaseMatrix:
    counts = []

    for idx, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        print("\r Reading " + str(idx), end="", flush=True)
        if sequences is not None and record.id not in sequences:
            continue

        seq = Sequence(str(record.seq), alphabet)
        counts.append(seq.kmer_count(k))

    return DatabaseMatrix(np.array(counts, dtype=np.uint8))


def count_kmers(
        bin_folder: str,
        file_extension: str = "fna",
        alphabet: Alphabet = PROTEIN,
        k: int = 4
) -> Tuple[QueryMatrix, List[str]]:
    bins = [file.rpartition(".")[0] for file in os.listdir(bin_folder) if file.rpartition(".")[2] == file_extension]
    counts = []

    for bin_id in tqdm(bins, ncols=100, leave=False, desc="Counting Kmers"):
        seq_str = "$".join([str(record.seq) for record in SeqIO.parse(os.path.join(bin_folder, bin_id + "." + file_extension), "fasta")])
        seq = Sequence(seq_str, alphabet, bin_id)
        counts.append(seq.kmer_count(k))

    return QueryMatrix(np.array(counts, dtype=np.uint8)), bins


def read_fasta_file(filename: str, alphabet: Alphabet) -> tp.List[Sequence]:
    return [Sequence(str(seq_record.seq), alphabet, seq_record.id) for seq_record in SeqIO.parse(filename, "fasta")]


def sequences_to_query_matrix(sequences: tp.List[Sequence], k: int) -> QueryMatrix:
    return QueryMatrix(np.array([seq.kmer_count(k) for seq in sequences]))
