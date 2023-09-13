import _io
import os
import subprocess
from multiprocessing.pool import ThreadPool
from typing import List, Optional, Tuple, Dict

import numpy as np
import numpy.typing as npt
from Bio import SeqIO
from tqdm import tqdm

from .matrices import QueryMatrix, DatabaseMatrix

MAX_PFAM = 17126


# TODO: Differentiate between fasta file and fasta folder instead of DatabaseMatrix and QueryMatrix

def create_database_matrix(
        orf_bin: str,
        prot_bin: str,
        pfam_dir: str,
        model_dir: str,
        fasta_file: str,
        sequences: Optional[List[str]] = None
) -> DatabaseMatrix:
    process_orf = subprocess.Popen(
        orf_bin, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
    )

    process_prot = subprocess.Popen(
        [prot_bin, "-p", "-F", "hf", pfam_dir, model_dir],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=process_orf.stdout
    )

    count_pfams_async = ThreadPool(processes=1).apply_async(_count_pfams, (process_prot.stdout, False))

    for idx, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        print("\r Reading " + str(idx), end="", flush=True)
        if sequences is not None and record.id not in sequences:
            continue

        assert process_orf.stdin is not None  # MyPy

        process_orf.stdin.write(">" + record.id + "\n")
        process_orf.stdin.write(str(record.seq) + "\n")

    assert process_orf.stdin is not None  # MyPy

    process_orf.stdin.close()
    process_orf.wait()

    pfam_counts, sequences = count_pfams_async.get()

    result, errors = process_prot.communicate()

    if process_prot.returncode != 0:
        raise Exception(errors)

    return DatabaseMatrix(pfam_counts)


def count_pfams(
        orf_bin: str,
        prot_bin: str,
        pfam_dir: str,
        model_dir: str,
        bin_folder: str,
        file_extension: str = "fna"
) -> Tuple[QueryMatrix, List[str]]:
    """
    This function takes a directory with bins in FASTA format and creates a Pfam count matrix. Each FASTA file is
    considered to be one bin. Sequence headers inside the files (most likely contig ids) are ignored.

    :param orf_bin: Path to uproc-orf binary
    :param prot_bin: Path to uproc-prot binary
    :param pfam_dir: Path to the UProC Pfam database directory
    :param model_dir: Path to UProC model directory
    :param bin_folder: Bin folder
    :param file_extension: File extension of the bin FASTA files. Probably something like .fna or .fasta. Each file in
    the bin folder that has this extension is considered a bin.
    :return: A 2-tuple where the first element is a QueryMatrix containing the Pfam counts. Each row represents a bin
    and each column a Pfam. The second element of the tuple is a list of bin names (names of the input FASTA files
    without file extension) in the same order as they appear in the QueryMatrix.
    """
    process_orf = subprocess.Popen(
        orf_bin, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
    )

    assert process_orf.stdin is not None  # MyPy

    process_prot = subprocess.Popen(
        [prot_bin, "-p", "-F", "hf", pfam_dir, model_dir],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=process_orf.stdout
    )

    count_pfams_async = ThreadPool(processes=1).apply_async(_count_pfams, (process_prot.stdout,))

    bins = [file.rpartition(".")[0] for file in os.listdir(bin_folder) if file.rpartition(".")[2] == file_extension]

    for bin_id in tqdm(bins, ncols=100, desc="- Counting Pfams"):
        for record in SeqIO.parse(os.path.join(bin_folder, bin_id + "." + file_extension), "fasta"):
            process_orf.stdin.write(">" + bin_id + "$$" + record.id + "\n")
            process_orf.stdin.write(str(record.seq) + "\n")

    process_orf.stdin.close()
    process_orf.wait()

    pfam_counts, sequences = count_pfams_async.get()

    result, errors = process_prot.communicate()

    if process_prot.returncode != 0:
        raise Exception(errors)

    return QueryMatrix(pfam_counts), sequences


def _count_pfams(stdout: _io.BufferedReader, merge: bool = True) -> Tuple[npt.NDArray[np.uint8], List[str]]:
    pfams: Dict[str, List[int]] = {}
    sequences = []

    for line in iter(stdout.readline, ''):
        line = line.decode("utf-8")
        if line == "":
            break

        seq, pfam = line.split(",")[:2]
        if merge:
            seq = seq.rpartition("$$")[0]
        pfam = int(pfam.strip()[2:])

        if seq not in pfams:
            pfams[seq] = []
            sequences.append(seq)

        pfams[seq].append(pfam)

    count_mat = np.zeros((len(sequences), MAX_PFAM + 1), dtype=np.uint8)
    for idx, seq in enumerate(sequences):
        for pfam in pfams[seq]:
            if count_mat[idx, pfam] == 255:
                continue
            count_mat[idx, pfam] += 1

    return count_mat, sequences
