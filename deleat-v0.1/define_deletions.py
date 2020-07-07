#!/usr/bin/env python3
"""define_deletions.py

    Compute a list of proposed deletions according to predicted
    essentiality scores.

@author: Jimena Solana
"""

from argparse import ArgumentParser
import os

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
import pandas as pd


NONCODING_MARGIN = 200


def get_deletions(gb_m1, l, e):
    """Define a list of deletions according to parameters L and E.

    A deletion is defined as any region longer than (or equal to) L not
    containing any gene with essentiality score higher than E.
    Args:
        gb_m1 (Bio.SeqRecord.SeqRecord): GenBank annotation.
        l (int): minimum deletion length.
        e (float): gene essentiality threshold.
    Returns:
        deletions (Bio.SeqFeature.CompoundLocation): list of proposed
            deletions.
    """
    essential_genes = [gene for gene in gb_m1.features
                       if is_essential(gene, e)]
    print("Essential gene count with threshold %.3f: %d"
          % (e, len(essential_genes)))
    essential_regions = CompoundLocation([
        FeatureLocation(gene.location.start - NONCODING_MARGIN,
                        gene.location.end + NONCODING_MARGIN)
        for gene in essential_genes
    ])
    essential_regions = merge_overlaps(essential_regions)
    nonessential_regions = complementary_compoundloc(0, len(gb_m1),
                                                     essential_regions)
    try:
        deletions = CompoundLocation([
            region for region in nonessential_regions.parts if len(region) >= l
        ])
    except ValueError:
        raise SystemExit("No deletions defined! Try readjusting L and/or E.")
    return deletions


def is_essential(feature, threshold):
    """Check whether a gene is essential."""
    if ("essentiality" in feature.qualifiers and
            float(feature.qualifiers["essentiality"][0]) > threshold):
        essential = True
    else:
        essential = False
    return essential


def merge_overlaps(comploc):
    """Merge overlapping regions in a CompoundLocation object.

    Args:
        comploc (Bio.SeqFeature.CompoundLocation): initial
            CompoundLocation.
    Returns:
        comploc_merged (Bio.SeqFeature.CompoundLocation):
            CompoundLocation without any overlapping parts.
    """
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
    """Get the complementary of a CompoundLocation (opposite parts).

    Args:
        start (int): start position of the genome.
        end (int): end position of the genome.
        comploc (Bio.SeqFeature.CompoundLocation): original
            CompoundLocation.
    Returns:
        complementary (Bio.SeqFeature.CompoundLocation): complementary
            CompoundLocation.
    """
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


def save_genbank_m2(deletions, gb_m1, gb_m2, l, e):
    """Generate modified-II GenBank file with annotated deletions.

    Args:
        deletions (Bio.SeqFeature.CompoundLocation): list of deletions.
        gb_m1 (Bio.SeqRecord.SeqRecord): modified-I GenBank annotation.
        gb_m2 (str): modified-II GenBank file path.
        l (int): minimum deletion length.
        e (float): gene essentiality threshold.
    """
    for n, deletion in enumerate(deletions.parts):
        deletion_feature = SeqFeature(
            location=deletion,
            type="misc_feature"
        )
        deletion_feature.qualifiers["note"] = ["deletion D%s" % (n+1),
                                               "L=%d" % l,
                                               "E=%.3f" %e]
        gb_m1.features.append(deletion_feature)
    SeqIO.write(gb_m1, gb_m2, "genbank")


def make_table(deletions, gb_m1):
    """Generate a summary table of deletion statistics.

    Args:
        deletions (Bio.SeqFeature.CompoundLocation): list of deletions.
        gb_m1 (Bio.SeqRecord.SeqRecord): modified-I GenBank annotation.
    Returns:
        table (pd.DataFrame): table of deletion statistics.
    """
    table = pd.DataFrame(
        index=["D%d" % (n+1) for n in range(len(deletions.parts))],
        columns=["start", "end", "length", "% of genome",
                 "contains (pseudo/hypot/non-hypot)"]
    )
    for n, deletion in enumerate(deletions.parts):
        table.iloc[n, 0] = deletion.start + 1  # 0-based > 1-based
        table.iloc[n, 1] = deletion.end
        length = deletion.end - deletion.start
        table.iloc[n, 2] = length
        table.iloc[n, 3] = length / len(gb_m1) * 100
        table.iloc[n, 4] = gene_content(deletion, gb_m1)
    return table


def gene_content(deletion, annot):
    """Describe the gene content of a deletion.

    For a specific deletions, returns the count of pseudo-genes, genes
    annotated as "hypothetical protein" and genes with functional
    annotation.
    Args:
        deletion (Bio.SeqFeature.FeatureLocation): a deletion.
        annot (Bio.SeqRecord.SeqRecord): modified-I GenBank annotation.
    Returns:
        result (str): gene contents (pseudo/hypot/non-hypot).
    """
    genes = {"pseudo": 0, "hypot": 0, "non-hypot": 0}
    for gene in annot.features:
        if gene.location.start in deletion or gene.location.end in deletion:
            if gene.type in ("tRNA", "rRNA", "tmRNA", "ncRNA"):
                genes["non-hypot"] += 1
            elif gene.type == "CDS" and "locus_tag" in gene.qualifiers:
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
    GENBANK_M1 = args.GBM1
    OUT_DIR = args.OUT_DIR
    os.makedirs(OUT_DIR, exist_ok=True)
    genbank_id = os.path.splitext(os.path.basename(GENBANK_M1))[0]
    GENBANK_M2 = os.path.join(OUT_DIR, genbank_id + ".gbm2")
    OUT_TABLE = os.path.join(OUT_DIR, "proposed_deletions.csv")
    L = args.DEL_LENGTH
    E = args.ESS_THRESHOLD

    # Check input
    try:
        annotation = SeqIO.read(GENBANK_M1, "genbank")
    except (FileNotFoundError, ValueError):
        raise SystemExit("\n\terror: could not read annotation file\n")
    try:
        for gene in annotation.features:
            if gene.type == "CDS":
                _ = gene.qualifiers["essentiality"][0]
                break
    except KeyError:
        raise SystemExit("\n\terror: invalid GenBank file (must be .gbm1)\n")
    if E < 0 or E > 1:
        raise SystemExit("\n\terror: invalid ESS_THRESHOLD")
    if L < 1 or L > len(annotation):
        raise SystemExit("\n\terror: invalid DEL_LENGTH")

    # Define deletions
    print("Computing deletion list...")
    proposed_deletions = get_deletions(annotation, L, E)
    save_genbank_m2(proposed_deletions, annotation, GENBANK_M2, L, E)
    deletions_table = make_table(proposed_deletions, annotation)
    deletions_table.to_csv(OUT_TABLE)
    print("Done. Results in %s and %s." % (OUT_TABLE, GENBANK_M2))
