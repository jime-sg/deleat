#!/usr/bin/env python3
"""parse_deg.py
@author: Jimena Solana
"""

from os import makedirs

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def make_header(sequence, annotation):
    gene_name = annotation.loc[sequence.id, "gene_name"]
    function = annotation.loc[sequence.id, "function"]
    header = "%s (%s)" % (function, gene_name)
    return header


if __name__ == "__main__":
    DEG_DIR = "/home/jimena/Bartonella/DEGdb/deg-p-15.2/"
    ANNOTATION = DEG_DIR + "degannotation-p.dat"
    FASTA = DEG_DIR + "degseq-p.dat"
    ORGANISMS = DEG_DIR + "organisms.txt"
    OUT_DIR = "/home/jimena/Bartonella/DEGdb/deg_byorg/"
    makedirs(OUT_DIR, exist_ok=True)

    with open(ORGANISMS) as f:
        organisms = {
            org.strip().split("\t")[0]: []
            for org in f if not org.startswith("#")
        }

    all_annot = pd.read_table(
        ANNOTATION, sep="\t", skiprows=1, index_col="deg_id",
        names=(
            "deg_org", "deg_id", "gene_name", "gene_ref", "cog",
            "class", "function", "organism", "refseq", "condition",
            "locus_tag", "go", "x"
        )
    )
    annot = all_annot.loc[all_annot["deg_org"].isin(organisms.keys())]

    deg_seqs = SeqIO.parse(FASTA, "fasta")
    for seq_ in deg_seqs:
        org = seq_.id[:7]
        if org in organisms.keys():
            new_seq = SeqRecord(
                seq=seq_.seq, id=seq_.id,
                description=make_header(seq_, annot)
            )
            organisms[org].append(new_seq)

    for organism in organisms:
        SeqIO.write(
            organisms[organism], OUT_DIR + organism + ".faa", "fasta"
        )
