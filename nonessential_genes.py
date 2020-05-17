#!/usr/bin/env python3
"""nonessential_genes.py
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import seqfeatures
import geptop


def save_proteome(gb_file, out_file):
    proteins = []
    annot = SeqIO.read(gb_file, "genbank")
    for feature in annot.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            protein = SeqRecord(
                seq=Seq(feature.qualifiers["translation"][0]),
                id=make_id(feature),
                description=make_description(feature)
            )
            proteins.append(protein)
    SeqIO.write(proteins, out_file, "fasta")


def make_id(feature):
    if "locus_tag" in feature.qualifiers:
        id_ = feature.qualifiers["locus_tag"][0]
    else:
        id_ = "(no locus_tag)"
    return id_


def make_description(feature):
    if "product" in feature.qualifiers:
        product = feature.qualifiers["product"][0]
    else:
        product = "(no product)"
    location = "".join(str(part) for part in feature.location.parts)
    description = "%s %s" % (product, location)
    return description


def integrate_scores(seqfeature_s, geptop_s):  # FIXME
    a = int(seqfeature_s.split(";")[0].split(",")[1])
    b = int(seqfeature_s.split(";")[1].split(",")[1])
    c = float(seqfeature_s.split(";")[2].split(",")[1])
    d = geptop_s
    score = a + b + c + d
    return score


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="nonessential-genes",
        description=""  # FIXME
    )
    parser.add_argument(
        "-g0", dest="GB", required=True,
        help="original GenBank annotation file")
    parser.add_argument(
        "-d", dest="DEG", required=True,
        help="")  # TODO
    parser.add_argument(
        "-c", dest="CV", required=True,
        help="")  # TODO
    parser.add_argument(
        "-g", dest="GEPTOP_CUTOFF", required=True, type=float,  # FIXME: "geptop"?
        help="")  # TODO
    parser.add_argument(
        "-n", dest="NPROC", required=True, type=int,
        help="")  # TODO
    parser.add_argument(
        "-o", dest="OUT_DIR", required=True,
        help="directory for output files")
    args, unknown = parser.parse_known_args()
    GENBANK = args.GB
    OUT_DIR = args.OUT_DIR
    genbank_id = os.path.splitext(os.path.basename(GENBANK))[0]
    GENBANK_M1 = os.path.join(OUT_DIR, genbank_id + ".gbm1")
    PROTEOME = os.path.join(OUT_DIR, "proteome.faa")
    DEG = args.DEG
    CV = args.CV
    GEPTOP_CUTOFF = args.GEPTOP_CUTOFF
    NPROC = args.NPROC

    # Check input
    # TODO

    save_proteome(gb_file=GENBANK, out_file=PROTEOME)

    # Add seqfeatures results to annotation file
    genbank_m1 = seqfeatures.run(GENBANK)

    # Add Geptop results to annotation file
    geptop_results = geptop.run(
        query_file=PROTEOME,
        deg_path=DEG,
        cv_path=CV,
        cutoff=GEPTOP_CUTOFF,
        n_proc=NPROC,
        out_path=OUT_DIR
    )
    for gene in genbank_m1.features:
        if gene.type == "CDS" and "translation" in gene.qualifiers:
            locus_tag = gene.qualifiers["locus_tag"][0]
            geptop_score = geptop_results[locus_tag][0]
            gene.qualifiers["geptop"] = [geptop_score]

            # Add integrated score to annotation file
            seqfeature_scores = gene.qualifiers["seqfeatures"][0]
            gene.qualifiers["essentiality"] = [integrate_scores(
                seqfeature_scores,
                geptop_score
            )]
        elif gene.type in ("tRNA", "rRNA", "tmRNA", "ncRNA"):
            gene.qualifiers["essentiality"] = 1
        elif gene.type == "CDS":
            gene.qualifiers["essentiality"] = 0


    # Save modified-I GenBank file
    SeqIO.write(genbank_m1, GENBANK_M1, "genbank")
