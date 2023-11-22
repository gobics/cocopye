# Copyright 2023 Niklas Birth
#
# This file is part of CoCoPyE.
#
# CoCoPyE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CoCoPyE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CoCoPyE. If not, see <https://www.gnu.org/licenses/>.

import _io
import os
import subprocess
from datetime import datetime
from multiprocessing.pool import ThreadPool
from typing import List, Tuple, Dict

import numpy as np
import numpy.typing as npt
from Bio import SeqIO
from tqdm import tqdm

MAX_PFAM = 17126


def count_pfams(
        orf_bin: str,
        prot_bin: str,
        pfam_dir: str,
        model_dir: str,
        bin_folder: str,
        file_extension: str = "fna",
        num_threads: int = 8,
        print_progress: bool = True
) -> Tuple[npt.NDArray[np.uint8], List[str], List[float]]:
    """
    This function takes a directory with bins in FASTA format and creates a Pfam count matrix. Each FASTA file is
    considered to be one bin. Sequence headers inside the files (most likely contig ids) are ignored.

    :param orf_bin: Path to uproc-orf binary
    :param prot_bin: Path to uproc-prot binary
    :param pfam_dir: Path to the UProC Pfam database directory
    :param model_dir: Path to UProC model directory
    :param bin_folder: Folder containing input bins in FASTA format
    :param file_extension: File extension of the input FASTA files. Probably something like .fna or .fasta. Each file in
    the bin folder that has this extension is considered a bin.
    :param num_threads: Number of threads that UProC should use. It is possible (and likely) that UProC ignores this
    parameter, but you can try.
    :param print_progress: Print a progress bar to stdout
    :return: A 3-tuple. The first element is a QueryMatrix containing the Pfam counts. Each row represents a bin
    and each column a Pfam. The second element is a list of bin names (names of the input FASTA files without file
    extension) in the same order as they appear in the QueryMatrix. The third element is a list of count-ratios of the
    input bins (number of pfams divided by bin size).
    """
    process_orf = subprocess.Popen(
        orf_bin, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
    )

    assert process_orf.stdin is not None  # MyPy

    process_prot = subprocess.Popen(
        [prot_bin, "-p", "-F", "hf", "-t", str(num_threads), pfam_dir, model_dir],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=process_orf.stdout
    )

    count_pfams_async = ThreadPool(processes=1).apply_async(_count_pfams, (process_prot.stdout,))

    bins = [file.rpartition(".")[0] for file in os.listdir(bin_folder) if file.rpartition(".")[2] == file_extension]

    lengths = {}

    for bin_id in tqdm(
            bins,
            ncols=0,
            desc="\033[0;37m[" + str(datetime.now()) + "]\033[0m Counting Pfams",
            disable=not print_progress
    ):
        for record in SeqIO.parse(os.path.join(bin_folder, bin_id + "." + file_extension), "fasta"):
            lengths[bin_id] = len(str(record.seq))
            process_orf.stdin.write(">" + bin_id + "$$" + record.id + "\n")
            process_orf.stdin.write(str(record.seq) + "\n")

    process_orf.stdin.close()
    process_orf.wait()

    pfam_counts, sequences, total_counts = count_pfams_async.get()
    count_ratio = [total_counts[seq] / lengths[seq] for seq in sequences]

    result, errors = process_prot.communicate()

    if process_prot.returncode != 0:
        raise Exception(errors)

    return pfam_counts, sequences, count_ratio


def _count_pfams(
        stdout: _io.BufferedReader,
        merge: bool = True
) -> Tuple[npt.NDArray[np.uint8], List[str], Dict[str, int]]:
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

    total_counts = {key: len(pfams[key]) for key in pfams}

    return count_mat, sequences, total_counts
