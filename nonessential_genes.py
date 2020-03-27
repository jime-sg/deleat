#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser
from sys import argv  # FIXME

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


def integrate_scores(seqfeature_s, geptop_s):
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
    # GENBANK = "/home/jimena/Bartonella/NC_005955.gb"
    # PROTEOME = "/home/jimena/Escritorio/proteoma.faa"
    GENBANK = argv[1]
    PROTEOME = argv[2]
    OUT_GENBANK = argv[3]

    save_proteome(
        gb_file=GENBANK,
        out_file=PROTEOME
    )

    # Add seqfeatures results to annotation file
    mod_annotation = seqfeatures.run(GENBANK)

    # Add Geptop results to annotation file
    geptop_results = geptop.run(
        query_file=PROTEOME,
        deg_path="/home/jimena/Bartonella/DEGdb/deg_byorg/temp",
        cv_path="/home/jimena/Bartonella/DEGdb/cv",
        cutoff=0.24,
        n_proc=4,
        out_path="/home/jimena/Escritorio"
    )  # FIXME (leer argumentos de argparser)
    for gene in mod_annotation.features:
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
        # TODO elif -> otros tipos de genes

    # Write modified annotation file
    SeqIO.write(mod_annotation, OUT_GENBANK, "genbank")
