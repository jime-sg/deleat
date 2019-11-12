#!/usr/bin/env python3
"""
@author: Jimena Solana
"""

import subprocess
from os import makedirs


def run(genome_file, repeat_length, out_path):
    # Vmatch: index genome
    makedirs("%sindex" % out_path, exist_ok=True)
    #mkvtree = subprocess.check_output("which mkvtree", text=True, shell=True).strip()
    subprocess.call(
        ["mkvtree", "-db", genome_file,
         "-dna", "-indexname", "%sindex/index" % out_path, "-pl", "-allout"]
    )
    # Vmatch: find repeats of length >= repeat_length
    subprocess.call(
        ["vmatch", "-l", str(repeat_length), "-d", "-p",
         "-noidentity", "-nodist", "-noevalue", "-noscore",
         "%sindex/index" % out_path, ">", "%srepeats.txt" % out_path]
    )
    repeats_list = parse_results("%srepeats.txt" % out_path)
    return repeats_list


def parse_results(results_file):
    repeats_list = []
    with open(results_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                repeat_s = int(line.split()[2])-1
                repeat_e = repeat_s + int(line.split()[0])
                repeats_list.append((repeat_s, repeat_e))
    return repeats_list