#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser

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


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="nonessential-genes",
        description=""  # FIXME
    )
    GENBANK = "/home/jimena/Bartonella/NC_005955.gb"
    PROTEOME = "/home/jimena/Escritorio/proteoma.faa"

    save_proteome(
        gb_file=GENBANK,
        out_file=PROTEOME
    )

    # seqfeatures_results = seqfeatures.run(GENBANK)
    # geptop_results = geptop.run(PROTEOME, )  # FIXME



