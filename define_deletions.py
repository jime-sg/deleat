#!/usr/bin/env python3
"""define_deletions.py
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser
import os

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
import pandas as pd


NONCODING_MARGIN = 200


def get_deletions(gb_m1, l, e):
    essential_genes = [gene for gene in gb_m1.features
                       if is_essential(gene, e)]
    essential_regions = CompoundLocation([
        FeatureLocation(gene.location.start - NONCODING_MARGIN,
                        gene.location.end + NONCODING_MARGIN)
        for gene in essential_genes
    ])
    essential_regions = merge_overlaps(essential_regions)
    nonessential_regions = complementary_compoundloc(0, len(gb_m1),
                                                     essential_regions)
    deletions = CompoundLocation([
        region for region in nonessential_regions.parts if len(region) >= l
    ])
    return deletions


def is_essential(feature, threshold):
    if (feature.type == "CDS" and
            float(feature.qualifiers["essentiality"][0]) > threshold):
        essential = True
    else:
        essential = False
    return essential


def merge_overlaps(comploc):
    parts = [comploc.parts[0]]
    for current in comploc.parts:
        previous = parts[-1]
        if previous.end >= current.start:
            parts[-1] = FeatureLocation(previous.start, current.end)
        else:
            parts.append(current)
    comploc_merged = CompoundLocation(parts)
    return comploc_merged


def complementary_compoundloc(start, end, comploc):
    complementary = FeatureLocation(start, comploc.parts[0].start)
    for i in range(len(comploc.parts) - 1):
        if comploc.parts[i+1].start - 1 < comploc.parts[i].end + 1:
            continue
        else:
            complementary += FeatureLocation(
                comploc.parts[i].end,
                comploc.parts[i+1].start
            )
    else:
        complementary += FeatureLocation(comploc.parts[-1].end, end)
    return complementary


def save_genbank_m2(deletions, gb_m1, gb_m2):
    for n, deletion in enumerate(deletions.parts):
        deletion_feature = SeqFeature(
            location=deletion,
            type="misc_feature"
        )
        deletion_feature.qualifiers["note"] = ["deletion D%s" % (n+1)]
        gb_m1.features.append(deletion_feature)
    SeqIO.write(gb_m1, gb_m2, "genbank")


def make_table(deletions, gb_m1):
    table = pd.DataFrame(
        index=["D%d" % (n+1) for n in range(len(deletions.parts))],
        columns=["start", "end", "length", "% of genome",
                 "contains (pseudo/hypot/non-hypot)"]
    )
    for n, deletion in enumerate(deletions.parts):
        table.iloc[n, 0] = deletion.start + 1  # 0-based -> 1-based
        table.iloc[n, 1] = deletion.end
        length = deletion.end - deletion.start
        table.iloc[n, 2] = length
        table.iloc[n, 3] = length / len(gb_m1) * 100
        table.iloc[n, 4] = gene_content(deletion, gb_m1)
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
        description=(
            "Compute table of proposed deletions. "
            "A deletion is defined as any region longer than (or equal to) "
            "DEL_LENGTH not containing any gene with essentiality score "
            "higher than ESS_THRESHOLD."
        )
    )
    parser.add_argument(
        "-g1", dest="GBM1", required=True,
        help="modified-I GenBank file")
    parser.add_argument(
        "-o", dest="OUT_DIR", required=True,
        help="directory for output files")
    parser.add_argument(
        "-l", dest="DEL_LENGTH", required=True, type=int,
        help="minimum deletion length")
    parser.add_argument(
        "-e", dest="ESS_THRESHOLD", required=True, type=float,
        help="gene essentiality threshold")
    args, unknown = parser.parse_known_args()
    # GENBANK_M1 = args.GBM1
    # OUT_DIR = args.OUT_DIR
    # genbank_id = os.path.splitext(os.path.basename(GENBANK_M1))[0]
    # GENBANK_M2 = os.path.join(OUT_DIR, genbank_id + ".gbm2")
    # OUT_TABLE = os.path.join(OUT_DIR, "proposed_deletions.csv")
    L = args.DEL_LENGTH
    E = args.ESS_THRESHOLD
    GENBANK_M1 = "/home/jimena/Escritorio/NC_005955.gbm1"  # FIXME
    GENBANK_M2 = "/home/jimena/Escritorio/NC_005955.gbm2"  # FIXME
    OUT_TABLE = "/home/jimena/Escritorio/deleciones.csv"  # FIXME

    # Check input
    # TODO

    annotation = SeqIO.read(GENBANK_M1, "genbank")

    proposed_deletions = get_deletions(annotation, L, E)
    save_genbank_m2(proposed_deletions, annotation, GENBANK_M2)

    deletions_table = make_table(proposed_deletions, annotation)
    deletions_table.to_csv(OUT_TABLE)
