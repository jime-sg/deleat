#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

import subprocess
from os import path, makedirs


def run(genome_file, repeat_length, out_path):
    """"""  # TODO
    index_files = (
        "index.al1", "index.bwt", "index.lcp", "index.ois", "index.sds",
        "index.sti1", "index.tis", "index.bck", "index.des", "index.llv",
        "index.prj", "index.skp", "index.suf"
    )
    if not all([path.isfile("%s/vmatch/%s" % (out_path, f))
                for f in index_files]):
        makedirs("%s/vmatch" % out_path, exist_ok=True)
        # Vmatch: index genome
        subprocess.call(
            ["mkvtree", "-db", genome_file, "-dna",
             "-indexname", "%s/vmatch/index" % out_path, "-pl", "-allout"]
        )
    # Vmatch: find repeats of length >= repeat_length
    with open("%s/vmatch/repeats.txt" % out_path, "w") as f:
        subprocess.call(
            ["vmatch", "-l", str(repeat_length), "-d", "-p",
             "-noidentity", "-nodist", "-noevalue", "-noscore",
             "%s/vmatch/index" % out_path],
            stdout=f
        )
    repeats_list = parse_results("%s/vmatch/repeats.txt" % out_path)
    return repeats_list


def parse_results(results_file):
    """"""  # TODO
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
    return set(repeats_list)