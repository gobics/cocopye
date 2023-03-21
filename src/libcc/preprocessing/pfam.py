import subprocess
import numpy as np

from ..matrices import QueryMatrix


def uproc(uproc_bin: str, pfam_dir: str, model_dir: str, infile: str, merge: bool = True) -> QueryMatrix:
    """
    This is a simple wrapper function for uproc-prot. If uproc exits with a non-zero error code, an exception is raised
    containing uprocs stderr output. If the binary cannot be found at the given path, ... (TODO).
    :param uproc_bin: Path to uproc-prot binary
    :param pfam_dir: Path to pfam database directory
    :param model_dir: Path to model directory
    :param infile: input file containing sequences in FASTA format
    :param merge: Merge pfam counts for sequences with the same name before the last underscore (e.g. `seqname_1` and
    `seqname_2`). This can be useful if the input file was generated by prodigal (which ist probably the case if you are
    using this function).
    :return: A count matrix containing the number of all pfams. Each row represents a sequence (or a set of sequences
    `seqname_*` in case or `merge=True`). They are ordered by first occurence in the input file. Each column represents
    a pfam. The column index corresponds to the internal representation (TODO: change this to PF00X index).
    """
    process = subprocess.Popen(
        [uproc_bin, "-p", "-n", "-F", "hf", pfam_dir, model_dir, infile],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    result, errors = process.communicate()

    if process.returncode != 0:
        raise Exception(errors)

    return QueryMatrix(_count_pfams(result.decode("utf8"), merge))


def _count_pfams(uproc_result: str, merge: bool = True):
    pfams = {}
    max_pfam = -1
    sequences = []

    for line in uproc_result.split("\n"):
        line = line.strip()
        if line == "":
            break

        seq, pfam = line.split(",")[:2]
        pfam = int(pfam)
        if merge:
            seq = seq.rpartition("_")[0]

        if seq not in pfams:
            pfams[seq] = []
            sequences.append(seq)

        pfams[seq].append(pfam)

        if pfam > max_pfam:
            max_pfam = pfam

    count_mat = np.zeros((len(sequences), max_pfam+1))
    for idx, seq in enumerate(sequences):
        for pfam in pfams[seq]:
            count_mat[idx, pfam] += 1

    return count_mat


def pfams_from_fasta():
    pass  # todo
    # Create a temporary directory, call prodigal and uproc and clean up at the end