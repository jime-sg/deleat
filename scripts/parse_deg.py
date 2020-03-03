#!/usr/bin/env python3
"""parse_deg.py
@author: Jimena Solana
"""

import os

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbimakeblastdbCommandline


def make_id(sequence, annotation):
    gi = annotation.loc[sequence.id, "gene_ref"]
    deg = sequence.id
    id_ = "gi|%s|deg|%s|" % (gi, deg)
    return id_


def make_description(sequence, annotation):
    gene_name = annotation.loc[sequence.id, "gene_name"]
    function = annotation.loc[sequence.id, "function"]
    description = "%s (%s)" % (function, gene_name)
    return description


def parse(annotation, fasta, organisms, prefix, out_dir):
    with open(organisms) as f:
        orgs = {
            org.strip().split("\t")[0]: []
            for org in f if not org.startswith("#")
        }

    if prefix == "DNEG":
        orgs = {k.replace("DEG", "DNEG"): v for k, v in orgs.items()}

    all_annot = pd.read_table(
        annotation,
        sep="\t", skiprows=1, low_memory=False, index_col="deg_id",
        names=(
            "deg_org", "deg_id", "gene_name", "gene_ref", "cog",
            "class", "function", "organism", "refseq", "condition",
            "locus_tag", "go", "x"
        )
    )
    annot = all_annot.loc[all_annot["deg_org"].isin(orgs.keys())]

    deg_seqs = SeqIO.parse(fasta, "fasta")
    for seq_ in deg_seqs:
        if prefix == "DEG":
            org = seq_.id[:7]
        elif prefix == "DNEG":
            org = seq_.id[:8]
        if org in orgs.keys():
            new_seq = SeqRecord(
                seq=seq_.seq,
                id=make_id(seq_, annot),
                description=make_description(seq_, annot)
            )
            orgs[org].append(new_seq)

    for org in orgs:
        SeqIO.write(
            orgs[org],
            os.path.join(out_dir, org + ".faa"), "fasta"
        )


def merge():
    orgs = os.listdir(os.path.join(OUT_DIR, "essential"))
    with open(os.path.join(OUT_DIR, "all", "deg.faa"), "w") as db:
        for org in orgs:
            with open(os.path.join(OUT_DIR, "all", org), "w") as f_all:
                for record in SeqIO.parse(
                        os.path.join(OUT_DIR, "essential", org), "fasta"
                ):
                    SeqIO.write(record, f_all, "fasta")
                    SeqIO.write(record, db, "fasta")

                org = org.replace("DEG", "DNEG")
                for record in SeqIO.parse(
                        os.path.join(OUT_DIR, "nonessential", org), "fasta"
                ):
                    SeqIO.write(record, f_all, "fasta")
                    SeqIO.write(record, db, "fasta")


if __name__ == "__main__":
    DEG_DIR = "/home/jimena/Bartonella/DEGdb/deg-p-15.2/"
    DNEG_DIR = "/home/jimena/Bartonella/DEGdb/deg-np-15.2/"
    DEG_ANNOTATION = DEG_DIR + "degannotation-p.dat"
    DNEG_ANNOTATION = DNEG_DIR + "degannotation-np.dat"
    DEG_FASTA = DEG_DIR + "degaa-p.dat"
    DNEG_FASTA = DNEG_DIR + "degaa-np.dat"
    ORGANISMS = "/home/jimena/Bartonella/DEGdb/organisms.txt"
    OUT_DIR = "/home/jimena/Bartonella/DEGdb/deg_byorg"
    os.makedirs(os.path.join(OUT_DIR, "essential"), exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "nonessential"), exist_ok=True)
    os.makedirs(os.path.join(OUT_DIR, "all"), exist_ok=True)

    # essential
    parse(
        DEG_ANNOTATION, DEG_FASTA, ORGANISMS, "DEG",
        os.path.join(OUT_DIR, "essential")
    )

    # nonessential
    parse(
        DNEG_ANNOTATION, DNEG_FASTA, ORGANISMS, "DNEG",
        os.path.join(OUT_DIR, "nonessential")
    )

    # all
    merge()

    # make blast db
    makeblastdb = NcbimakeblastdbCommandline(
        dbtype="prot",
        input_file=os.path.join(OUT_DIR, "all", "deg.faa")
    )
    out, err = makeblastdb()

