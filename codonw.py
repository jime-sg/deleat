#!/usr/bin/env python3
"""codonw.py
# TODO
@author: Jimena Solana
"""

import os
import subprocess

import pandas as pd


def run(fasta, feature_list, remove_blk=True):
    print("Running CodonW...")
    out_dir = os.path.split(fasta)[0]
    out = os.path.join(out_dir, "codonw.out")
    blk = os.path.join(out_dir, "codonw.blk")
    log = os.path.join(out_dir, "codonw_log.txt")
    with open(log, "w") as fo:
        subprocess.run(
            ["codonw", fasta, out, blk,
             "-nomenu", "-nowarn", "-silent", "-machine", "-t,",
             *feature_list],
            stdout=fo, stderr=fo
        )
    results = read_table(out)
    if remove_blk:
        os.remove(blk)
    print("Done.")
    return results


def read_table(csv):
    table = pd.read_csv(csv, index_col="title")
    table = table.iloc[:, :-1]  # Remove empty last column
    table.Nc = pd.to_numeric(table.Nc, errors="coerce")  # Missing values > NaN
    table.index = table.index.map(lambda x: x.split()[0])  # Index = locus_tag
    return table
