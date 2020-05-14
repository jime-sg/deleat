#!/usr/bin/env python3
"""define_deletions.py
@author: Jimena Solana
"""

from argparse import ArgumentParser

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pandas as pd


NONCODING_MARGIN = 200


def get_deletions(gb_m1, l, e):
    genbank_m1 = SeqIO.read(gb_m1, "genbank")
    essential_genes = [
        gene for gene in genbank_m1.features if is_essential(gene, e)
    ]
    essential_regions = CompoundLocation([
        FeatureLocation(gene.location.start - NONCODING_MARGIN,
                        gene.location.end + NONCODING_MARGIN)
        for gene in essential_genes
    ])
    nonessential_regions = complementary_compoundloc(
        0, len(genbank_m1), essential_regions
    )
    deletions = CompoundLocation([
        region for region in nonessential_regions if len(region) >= l
    ])
    return deletions


def is_essential(feature, threshold):
    if (feature.type == "CDS" and
            float(feature.qualifiers["essentiality"][0]) > threshold):  # FIXME: todos los CDS tendrán /essentiality??? comprobar
        essential = True
    else:
        essential = False
    return essential


def complementary_compoundloc(start, end, comploc):
    complementary = FeatureLocation(start, comploc.parts[0].start)
    for i in range(len(comploc.parts) - 1):
        print(comploc.parts[i], comploc.parts[i+1])
        complementary += FeatureLocation(
            comploc.parts[i].end + 1,
            comploc.parts[i+1].start - 1
        )
    else:
        complementary += FeatureLocation(comploc.parts[-1].end, end)
    return complementary


def save_genbank_m2(deletions, gb_m1, gb_m2):
    annot = SeqIO.read(gb_m1, "genbank")
    for n, deletion in enumerate(deletions.parts):
        deletion_feature = SeqFeature(
            location=deletion,
            type="misc_feature",
            id="deletion %s" % (n + 1)
        )
        annot.features.append(deletion_feature)
    SeqIO.write(annot, gb_m2, "genbank")


def make_table(deletions, gb_m1):
    annot = SeqIO.read(gb_m1, "genbank")
    table = pd.DataFrame(
        index=["D%d" % (n+1) for n in range(len(deletions))],
        columns=["start", "end", "length", "% of genome",
                 "contains (pseudo/hypot/non-hypot)"]
    )
    for n, deletion in enumerate(deletions.parts):
        name = "D%d" % (n+1)
        table.loc[name, "start"] = deletion.start
        table.loc[name, "end"] = deletion.end
        length = deletion.end - deletion.start + 1
        table.loc[name, "length"] = length
        table.loc[name, "% of genome"] = length / len(annot) * 100
        table.loc[name, "contains (pseudo/hypot/non-hypot)"] = gene_content(deletion, annot)
    return table


def gene_content(deletion, annot):
    genes = {"pseudo": 0, "hypot": 0, "non-hypot": 0}
    for gene in annot.features:
        if gene.location.start in deletion and gene.location.end in deletion:
            if gene.type in ("tRNA", "rRNA", "tmRNA", "ncRNA"):
                genes["non-hypot"] += 1
            elif gene.type == "CDS":
                if "pseudo" in gene.qualifiers:
                    genes["pseudo"] += 1
                elif ("product" in gene.qualifiers and
                      "hypothetical" in gene.qualifiers["product"][0]):
                    genes["hypot"] += 1
                else:
                    genes["non-hypot"] += 1
    result = "%d/%d/%d" % (genes["pseudo"], genes["hypot"], genes["non-hypot"])
    return result


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="define-deletions",
        description=""  # FIXME
    )
    GENBANK_M1 = "/home/jimena/Escritorio/NC_005955.gbm1"  # FIXME
    GENBANK_M2 = "/home/jimena/Escritorio/NC_005955.gbm2"  # FIXME
    OUT_TABLE = "/home/jimena/Escritorio/deleciones.csv"  # FIXME
    L = 10000  # FIXME
    E = 0.5  # FIXME

    proposed_deletions = get_deletions(GENBANK_M1, L, E)
    save_genbank_m2(proposed_deletions, GENBANK_M1, GENBANK_M2)

    deletions_table = make_table(proposed_deletions, GENBANK_M1)
    deletions_table.to_csv(OUT_TABLE)
