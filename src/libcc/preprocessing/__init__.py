"""
This submodule contains several classes and functions that can be used to read the input data and construct count
matrices.

... A bit of informative text ...
"""

import subprocess


def prodigal(prodigal_bin: str, infile: str, outfile: str, show_progress=False) -> int:
    """
    A wrapper function for prodigal. If the binary is not found at the given path, ...(TODO).
    :param prodigal_bin: Path to prodigal binary
    :param infile: Input file in FASTA format
    :param outfile: Output file (because we want the results to be in FASTA format, the output has to be written to a
    file)
    :param show_progress: If this is set to true, the stderr output of prodigal is printed to stdout.
    (TODO: Replace this by a better progress indicator)
    :return: The error code of the prodigal process.
    """
    process = subprocess.Popen(
        [prodigal_bin, '-p', 'single', '-m', '-i', infile, '-a', outfile],  # TODO mode? -m?
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
    )

    if show_progress:
        while process.poll() is None:
            print(process.stderr.readline(), end="")
        print(process.stderr.read(), end="")

    returncode = process.wait()

    return returncode
