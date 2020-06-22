#!/usr/bin/env python3
"""summarise.py
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser
import os

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord
import numpy as np

import circplot
from predict_essentiality import find_ori_ter
from define_deletions import complementary_compoundloc, is_essential
from revise_deletions import is_deletion


sep = "-" * 80 + "\n"


def save_genbank_m4(gb_m3, gb_m4):
    """
    # TODO
    Args:
        gb_m3:
        gb_m4:

    Returns:

    """
    deletions = []
    for feature in gb_m3.features:
        if is_deletion(feature):
            deletions.append(feature.location)
    deletions = CompoundLocation(deletions)
    non_deletions = complementary_compoundloc(0, len(gb_m3), deletions)
    reduced_annot = SeqRecord(
        seq=non_deletions.extract(gb_m3.seq),
        id=gb_m3.id,
        name=gb_m3.name,
        description=gb_m3.description,
        dbxrefs=gb_m3.dbxrefs,
        annotations=gb_m3.annotations
    )
    end = 0
    for nondel in non_deletions.parts:
        offset = nondel.start - end
        for feature in gb_m3.features:
            if (feature.location.start in nondel and
                    feature.location.end in nondel):
                feature.location = feature.location + (-offset)
                reduced_annot.features.append(feature)
        end = nondel.end - offset + 2
    SeqIO.write(reduced_annot, gb_m4, "genbank")


def get_stats(gb_m3):
    """
    # TODO
    Args:
        gb_m3:

    Returns:

    """
    l = 0
    e = 0
    deletions = []
    for feature in gb_m3.features:
        if is_deletion(feature):
            deletions.append(feature.location)
            if not (l and e):
                l = int(feature.qualifiers["note"][1].split("=")[-1])
                e = float(feature.qualifiers["note"][2].split("=")[-1])
    deletions = CompoundLocation(deletions)
    stats = {
        "org": gb_m3.annotations["organism"],
        "size": len(gb_m3),
        "l": l,
        "e": e,
        "del_n": 0,
        "del_kb": 0,
        "del_perc": 0,
        "del_genes_pseudo": 0,
        "del_genes_hypot": 0,
        "del_genes_nonhypot": 0,
        "ess_genes_rna": 0,
        "ess_genes_hypot": 0,
        "ess_genes_nonhypot": 0,
        "ness_genes_pseudo": 0,
        "ness_genes_hypot": 0,
        "ness_genes_nonhypot": 0
    }
    for feature in gb_m3.features:
        if is_essential(feature, e):  # essential
            if feature.type in ("tRNA", "rRNA", "tmRNA", "ncRNA"):
                stats["ess_genes_rna"] += 1
            elif feature.type == "CDS":
                if ("product" in feature.qualifiers and
                        "hypothetical" in feature.qualifiers["product"][0]):
                    stats["ess_genes_hypot"] += 1
                else:
                    stats["ess_genes_nonhypot"] += 1
        elif feature.type == "CDS":  # non-essential
            if "pseudo" in feature.qualifiers:
                stats["ness_genes_pseudo"] += 1
            elif ("product" in feature.qualifiers and
                  "hypothetical" in feature.qualifiers["product"][0]):
                stats["ness_genes_hypot"] += 1
            else:
                stats["ness_genes_nonhypot"] += 1
            if (feature.location.start in deletions or
                    feature.location.end in deletions):  # in a deletion
                if "pseudo" in feature.qualifiers:
                    stats["del_genes_pseudo"] += 1
                elif ("product" in feature.qualifiers and
                      "hypothetical" in feature.qualifiers["product"][0]):
                    stats["del_genes_hypot"] += 1
                else:
                    stats["del_genes_nonhypot"] += 1
        elif is_deletion(feature):
            stats["del_n"] += 1
            stats["del_kb"] += len(feature) / 1000
            stats["del_perc"] += len(feature) / len(gb_m3) * 100
    return stats


def write_stats(stats, log):
    with open(log, "w") as f:
        f.write(sep)
        f.write(
            "##### GENOME REDUCTION REPORT FOR %s (%d bp) #####\n"
            % (stats["org"], stats["size"])
        )
        f.write(sep)
        f.write(
            "Deletion design parameters:\n  L = %d\n  E = %.3f\n\n"
            % (stats["l"], stats["e"])
        )
        f.write(
            ("Predicted essential genes:\n"
             "  total %d\n  RNA %d\n  hypothetical %d\n  annotated %d\n\n")
            % (stats["ess_genes_rna"]
               + stats["ess_genes_hypot"]
               + stats["ess_genes_nonhypot"],
               stats["ess_genes_rna"],
               stats["ess_genes_hypot"],
               stats["ess_genes_nonhypot"])
        )
        f.write(
            ("Predicted non-essential genes:\n"
             "  total %d\n  pseudo %d\n  hypothetical %d\n  annotated %d\n\n")
            % (stats["ness_genes_pseudo"]
               + stats["ness_genes_hypot"]
               + stats["ness_genes_nonhypot"],
               stats["ness_genes_pseudo"],
               stats["ness_genes_hypot"],
               stats["ness_genes_nonhypot"])
        )
        f.write(
            "Designed deletions:\n  total %d\n  %d kb\n  %.2f %% of genome\n"
            % (stats["del_n"], stats["del_kb"], stats["del_perc"])
        )
        f.write(sep)


def best_deletion_order(gb_m3, ori, ter, log):
    """
    # TODO
    Args:
        gb_m3:
        ori:
        ter:
        log:

    Returns:

    """
    # Replichore imbalance = len(genome)/2 - len(0-180º replichore)
    if ter > ori:
        len_1 = ter - ori
    else:
        len_1 = len(gb_m3) - (ori - ter)
    imbalance = len(gb_m3)/2 - len_1
    # Calculate each deletion's contribution to the imbalance
    deletions = []
    for feature in gb_m3.features:
        if is_deletion(feature):
            name = feature.qualifiers["note"][0].split()[-1]
            midpoint = (feature.location.end + feature.location.start) / 2
            if midpoint > ori or midpoint < ter:  # 0-180º replichore
                contrib = len(feature)
            else:  # 180-360º replichore
                contrib = -len(feature)
            deletions.append((name, contrib))
    # Best next deletion is the one that minimises the imbalance
    with open(log, "a") as f:
        f.write("Initial replichore imbalance: %d\n" % imbalance)
        f.write("Best deletion order for minimising imbalance at each step:\n")
        n = 0
        while deletions:
            n += 1
            best = deletions.pop(
                np.argmin([abs(imbalance + deletion[1])
                           for deletion in deletions])
            )
            imbalance += best[1]
            f.write("  %d. %s -> imbalance = %d\n" % (n, best[0], imbalance))
        f.write(sep)


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="summarise",
        description=""  # FIXME
    )
    parser.add_argument(
        "-g3", dest="GBM3", required=True,
        help="modified-III GenBank file")
    parser.add_argument(
        "-o", dest="OUT_DIR", required=True,
        help="directory for output files")
    parser.add_argument(
        "-p1", dest="ORI", type=int,
        help="position of origin of replication")
    parser.add_argument(
        "-p2", dest="TER", type=int,
        help="position of terminus of replication")
    args, unknown = parser.parse_known_args()
    GENBANK_M3 = args.GBM3
    OUT_DIR = args.OUT_DIR
    os.makedirs(OUT_DIR, exist_ok=True)
    genbank_id = os.path.splitext(os.path.basename(GENBANK_M3))[0]
    GENBANK_M4 = os.path.join(OUT_DIR, genbank_id + ".gbm4")
    OUT_IMG = os.path.join(OUT_DIR, "genome_reduction")
    REPORT_LOG = os.path.join(OUT_DIR, "reduction_stats.txt")
    ORI = args.ORI
    TER = args.TER

    # Check input
    # TODO
    genbank_m3 = SeqIO.read(GENBANK_M3, "genbank")

    # Draw circular genome plot
    print("Drawing circular genome plot...")
    save_genbank_m4(genbank_m3, GENBANK_M4)
    circplot.plot(
        gb_outer=GENBANK_M3,
        gb_inner=GENBANK_M4,
        out_file=OUT_IMG,
        out_fmt="png"
    )
    os.remove(GENBANK_M4)

    # Final report
    print("Generating final reduction report...")
    reduction_stats = get_stats(genbank_m3)
    write_stats(reduction_stats, REPORT_LOG)
    if not (ORI and TER):
        ORI, TER = find_ori_ter(genbank_m3)
    best_deletion_order(genbank_m3, ORI, TER, REPORT_LOG)
    print("Done. Results in %s and %s." % (OUT_IMG + ".png", REPORT_LOG))
