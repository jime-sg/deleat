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
from define_deletions import complementary_compoundloc
from revise_deletions import is_deletion


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
    with open(log, "w") as f:
        f.write("Initial replichore imbalance: %d\n" % imbalance)
        f.write("Best deletion order for keeping the minimum possible "
                  "imbalance at each step:\n")
        while deletions:
            best = deletions.pop(
                np.argmin([abs(imbalance + deletion[1])
                for deletion in deletions])
            )
            imbalance += best[1]
            f.write("  %s -> imbalance = %d\n" % (best[0], imbalance))


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
        "-p1", dest="ORI",
        help="position of origin of replication")
    parser.add_argument(
        "-p2", dest="TER",
        help="position of terminus of replication")
    args, unknown = parser.parse_known_args()
    GENBANK_M3 = args.GBM3
    OUT_DIR = args.OUT_DIR
    os.makedirs(OUT_DIR, exist_ok=True)
    genbank_id = os.path.splitext(os.path.basename(GENBANK_M3))[0]
    GENBANK_M4 = os.path.join(OUT_DIR, genbank_id + ".gbm4")
    ORDER_LOG = os.path.join(OUT_DIR, "deletion_order.txt")
    OUT_IMG = os.path.join(OUT_DIR, "genome_reduction")
    ORI = args.ORI
    TER = args.TER

    # Check input
    # TODO
    genbank_m3 = SeqIO.read(GENBANK_M3, "genbank")

    # Determine best deletion order
    if not (ORI and TER):
        ORI, TER = find_ori_ter(genbank_m3)
    best_deletion_order(genbank_m3, ORI, TER, ORDER_LOG)

    # Draw circular genome plot
    save_genbank_m4(genbank_m3, GENBANK_M4)

    circplot.plot(
        gb_outer=GENBANK_M3,
        gb_inner=GENBANK_M4,
        out_file=OUT_IMG,
        out_fmt="png"
    )

    # TODO: remove reduced gb
