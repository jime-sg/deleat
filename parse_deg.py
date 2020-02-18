#!/usr/bin/env python3
"""
@author: Jimena Solana
"""

import pandas as pd

DEG_DIR = "/home/jimena/Bartonella/DEGdb/deg-p-15.2/"
ANNOTATION = DEG_DIR + "degannotation-p.dat"
FASTA = DEG_DIR + "degseq-p.dat"
ORGANISMS = DEG_DIR + "organisms.txt"

with open(ORGANISMS) as f:
    organisms = [org.strip() for org in f if not org.startswith("#")]

annotation = pd.read_table(ANNOTATION, sep="\t", header=0, dtype=str)
