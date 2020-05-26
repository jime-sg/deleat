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
import pandas as pd

import seqfeatures
import geptop


def extract_cds(gb_file, out_aa, out_nt):
    proteins_aa = []
    proteins_nt = []
    annot = SeqIO.read(gb_file, "genbank")
    for feature in annot.features:
        if (feature.type == "CDS" and
                "translation" in feature.qualifiers and
                "locus_tag" in feature.qualifiers):
            protein_aa = SeqRecord(
                seq=Seq(feature.qualifiers["translation"][0]),
                id=feature.qualifiers["locus_tag"][0],
                description=make_description(feature)
            )
            protein_nt = SeqRecord(
                seq=feature.extract(annot.seq),  # FIXME
                id=feature.qualifiers["locus_tag"][0],
                description=make_description(feature)
            )
            proteins_aa.append(protein_aa)
            proteins_nt.append(protein_nt)
    SeqIO.write(proteins_aa, out_aa, "fasta")
    SeqIO.write(proteins_nt, out_nt, "fasta")


def make_description(feature):
    if "product" in feature.qualifiers:
        product = feature.qualifiers["product"][0]
    else:
        product = "(no product)"
    location = "".join(str(part) for part in feature.location.parts)
    description = "%s %s" % (product, location)
    return description


def strand(proteome, ori, ter):  # FIXME: comprobar
    proteins = list(SeqIO.read(proteome, "fasta"))
    strand_results = pd.DataFrame(
        index=[prot.id for prot in proteins],
        columns=["strand_lead"]
    )
    for prot in proteins:
        location = prot.description.split()[-1]
        start = int(location.split(":")[0][1:])
        sign = location[-2]
        if start > ori or start < ter:  # 0-180º replichore
            if sign == "+":  # leading strand
                strand_results.loc[prot.id, "strand_lead"] = 1
            else:  # lagging strand
                strand_results.loc[prot.id, "strand_lead"] = 0
        else:  # 180-360º replichore
            if sign == "+":  # lagging strand
                strand_results.loc[prot.id, "strand_lead"] = 0
            else:  # leading strand
                strand_results.loc[prot.id, "strand_lead"] = 1
    return strand_results


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
    DEG = args.DEG
    CV = args.CV
    GEPTOP_CUTOFF = args.GEPTOP_CUTOFF
    NPROC = args.NPROC

    # Check input
    # TODO

    # Get ori + ter coordinates
    # TODO
    ori = 0
    ter = 500000

    # Extract proteome
    proteome_aa = os.path.join(OUT_DIR, "proteome.faa")
    proteome_nt = os.path.join(OUT_DIR, "proteome.fna")
    extract_cds(gb_file=GENBANK, out_aa=proteome_aa, out_nt=proteome_nt)

    # Get strand data
    results = strand(proteome_aa, ori=ori, ter=ter)

    # Get Geptop scores
    geptop_results = geptop.run(
        query_file=proteome_aa,
        deg_path=DEG,
        cv_path=CV,
        cutoff=GEPTOP_CUTOFF,
        n_proc=NPROC,
        out_path=OUT_DIR
    )
    for gene, score in geptop_results.items():
        results.loc[gene, "geptop"] = score[0]


    # Add seqfeatures results to annotation file
    genbank_m1 = seqfeatures.run(GENBANK)

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
