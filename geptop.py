#!/usr/bin/env python3
"""
@author: Jimena Solana
"""

from Bio import SeqIO


def run(cds_file, cutoff, n_proc, out_path):
    seqs = SeqIO.parse(open(cds_file), "fasta")



