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

import geptop
import codonw


CODONW_FEATURES = ["-enc", "-gc", "-sil_base", "-L_sym", "-L_aa", "-aro",
                   "-hyd"]


def get_feature_table(gb, out_dir, ori, ter, geptop_params, codonw_features):
    # Extract proteome
    proteome_aa = os.path.join(out_dir, "proteome.faa")
    proteome_nt = os.path.join(out_dir, "proteome.fna")
    extract_cds(gb_file=gb, out_aa=proteome_aa, out_nt=proteome_nt)

    # Strand data
    feature_table = strand(proteome_aa, ori=ori, ter=ter)

    # Geptop scores
    geptop_results = geptop.run(query_file=proteome_aa, **geptop_params)
    for gene, score in geptop_results.items():
        feature_table.loc[gene, "geptop"] = score[0]

    # CodonW features
    codonw_results = codonw.run(proteome_nt, codonw_features)
    feature_table.join(codonw_results)
    return results


def extract_cds(gb_file, out_aa, out_nt):
    proteins_aa = []
    proteins_nt = []
    annot = SeqIO.read(gb_file, "genbank")
    for feature in annot.features:
        if "essentiality" in feature.qualifiers:
            # Ignore already annotated genes
            continue
        if (feature.type == "CDS" and
                "translation" in feature.qualifiers and
                "locus_tag" in feature.qualifiers):
            protein_aa = SeqRecord(
                seq=Seq(feature.qualifiers["translation"][0]),
                id=feature.qualifiers["locus_tag"][0],
                description=make_description(feature)
            )
            protein_nt = SeqRecord(
                seq=feature.extract(annot.seq),
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
    # TODO: añadir ORI y TER (opcionales)
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
    ori = 0  # FIXME
    ter = 500000  # FIXME

    # Get table of all gene features
    geptop_params = {
        "deg_path": DEG,
        "cv_path": CV,
        "cutoff": GEPTOP_CUTOFF,
        "n_proc": NPROC,
        "out_path": OUT_DIR
    }
    results = get_feature_table(GENBANK, OUT_DIR, ori, ter,
                                geptop_params, CODONW_FEATURES)
    
    # TODO: predicción
    essentiality_scores = {}  # {locus_tag: p(essential)}

    # Create modified-I GenBank file
    annotation = SeqIO.read(GENBANK, "genbank")
    for gene in annotation.features:
        if (gene.type == "CDS" and
                "translation" in gene.qualifiers and
                "locus_tag" in gene.qualifiers):
            locus_tag = gene.qualifiers["locus_tag"][0]
            gene.qualifiers["essentiality"] = essentiality_scores[locus_tag]
        elif gene.type in ("tRNA", "rRNA", "tmRNA", "ncRNA"):
            gene.qualifiers["essentiality"] = 1
        elif gene.type == "CDS":  # pseudo-gene
            gene.qualifiers["essentiality"] = 0
    SeqIO.write(annotation, GENBANK_M1, "genbank")
