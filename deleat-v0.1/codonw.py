#!/usr/bin/env python3
"""codonw.py

    Interface for running CodonW with command-line options and parsing
    results.

@author: Jimena Solana
"""

import os
import subprocess

import pandas as pd


def run(fasta, feature_list, remove_blk=True):
    """Run CodonW for all genes in a FASTA file.

    Args:
        fasta (str): query genes FASTA (nt) file path.
        feature_list (list): list of features to calculate for each gene
            (as command-line arguments).
        remove_blk (bool): whether to remove output .blk file after
            completion.
    Returns:
        results (pd.DataFrame): table of results with each gene in a row
            and each calculated feature in a column.
    """
    print("  Running CodonW...")
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
    return results


def read_table(csv):
    """Parse CodonW .out file into a pd.DataFrame and clean data.

    Args:
        csv (str): CodonW .out file path.
    Returns:
        table (pd.DataFrame): CodonW results.
    """
    table = pd.read_csv(csv, index_col="title")
    table = table.iloc[:, :-1]  # Remove empty last column
    table.Nc = pd.to_numeric(table.Nc, errors="coerce")  # Missing values > NaN
    table.index = table.index.map(lambda x: x.split()[0])  # Index = locus_tag
    return table
