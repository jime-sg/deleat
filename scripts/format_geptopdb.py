#!/usr/bin/env python3
"""
@author: Jimena Solana
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

GEPTOPDB_DIR = "/home/jimena/Bartonella/Geptop2/datasets2"
GEPTOPDEG_FILE = "/home/jimena/Bartonella/Geptop2/DEG2"
OUT_DIR = "/home/jimena/Bartonella/Geptop2/datasets_f"

if __name__ == "__main__":
    with open(GEPTOPDEG_FILE, "r") as f:
        deg = [line.strip() for line in f]

    os.makedirs(OUT_DIR, exist_ok=True)
    orgs = os.listdir(GEPTOPDB_DIR)
    for org in orgs:
        with open(os.path.join(OUT_DIR, org), "w") as fo:
            for seq in SeqIO.parse(os.path.join(GEPTOPDB_DIR, org), "fasta"):
                if seq.id in deg:
                    desc = seq.description.replace("| ", "|DEG| ")
                else:
                    desc = seq.description.replace("| ", "|DNEG| ")

                SeqIO.write(
                    SeqRecord(seq.seq, id=desc, description=""),
                    fo, "fasta"
                )
