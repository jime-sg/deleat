#!/usr/bin/env python3
"""codonw.py
# TODO
@author: Jimena Solana
"""

import os
import subprocess

import pandas as pd


def run(fasta, feature_list, remove_blk=True):
    prefix = os.path.splitext(fasta)[0]
    out = prefix + ".out"
    blk = prefix + ".blk"
    subprocess.call(
        ["codonw", fasta,
         "-nomenu", "-nowarn", "-silent", "-machine", "-t,",
         *feature_list]
    )
    results = read_table(out, feature_list)
    if remove_blk:
        os.remove(blk)
    return results


def read_table(csv, columns):
    table = pd.read_csv(csv, usecols=["title", *columns], index_col="title")
    table.Nc = pd.to_numeric(table.Nc, errors="coerce")  # Missing values > NaN
    return table
