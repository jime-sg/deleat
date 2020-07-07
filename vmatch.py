#!/usr/bin/env python3
"""vmatch.py

    Interface for running Vmatch with command-line options and parsing
    results.

@author: Jimena Solana
"""

import subprocess
import os


def run(genome_file, repeat_length, out_dir):
    """Compute list of all exact repeats of length >= repeat_length in a
    genome, using Vmatch.

    First index the genome sequence by calling mkvtree (if it is not
    indexed yet). Then call vmatch to find all exact repeats of length
    >= repeat_length, both in the same and reverse complement strand.
    Args:
        genome_file (str): genome sequence file in FASTA format.
        repeat_length (int): min size of repeats to be found.
        out_dir (str): directory for output files.
    Returns:
        repeats_list (set of tuple(int, int)): found repeat locations.
    """
    index_files = (
        "index.al1", "index.bwt", "index.lcp", "index.ois", "index.sds",
        "index.sti1", "index.tis", "index.bck", "index.des", "index.llv",
        "index.prj", "index.skp", "index.suf"
    )
    if not all([os.path.isfile(os.path.join(out_dir, "vmatch", f))
                for f in index_files]):
        os.makedirs(os.path.join(out_dir, "vmatch"), exist_ok=True)
        # Vmatch: index genome
        subprocess.call(
            ["mkvtree", "-db", genome_file, "-dna",
             "-indexname", os.path.join(out_dir, "vmatch/index"),
             "-pl", "-allout"]
        )
    # Vmatch: find repeats
    repeats_file = os.path.join(out_dir, "vmatch/repeats.txt")
    with open(repeats_file, "w") as f:
        subprocess.call(
            ["vmatch", "-l", str(repeat_length), "-d", "-p",
             "-noidentity", "-nodist", "-noevalue", "-noscore",
             os.path.join(out_dir, "vmatch/index")],
            stdout=f
        )
    repeats_list = parse_results(repeats_file)
    return repeats_list


def parse_results(results_file):
    """Parse Vmatch output file into a set of repeat coordinates.

    Args:
        results_file (str): Vmatch output file path.
    Returns:
        results (set of tuple(int, int)): found repeat locations.
    """
    repeats_list = []
    with open(results_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                length = int(line.split()[0])
                repeat_s1 = int(line.split()[2]) + 1
                repeat_e1 = repeat_s1 + length - 1
                repeats_list.append((repeat_s1, repeat_e1))
                repeat_s2 = int(line.split()[6]) + 1
                repeat_e2 = repeat_s2 + length - 1
                repeats_list.append((repeat_s2, repeat_e2))
    results = set(repeats_list)
    return results
