#!/usr/bin/env python3
"""summarize.py
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser
import os

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord

import circplot
from define_deletions import complementary_compoundloc


def is_deletion(feature):
    if (feature.type == "misc_feature" and
            "note" in feature.qualifiers and
            "deletion" in feature.qualifiers["note"][0]):
        return True
    else:
        return False


def save_genbank_m4(gb_m3, gb_m4):
    annot = SeqIO.read(gb_m3, "genbank")
    deletions = []
    for feature in annot.features:
        if is_deletion(feature):
            deletions.append(feature.location)
    deletions = CompoundLocation(deletions)
    non_deletions = complementary_compoundloc(0, len(annot), deletions)
    reduced_annot = SeqRecord(
        seq=non_deletions.extract(annot.seq),
        id=annot.id,
        name=annot.name,
        description=annot.description,
        dbxrefs=annot.dbxrefs,
        annotations=annot.annotations
    )
    end = 0
    for nondel in non_deletions.parts:
        offset = nondel.start - end
        for feature in annot.features:
            if (feature.location.start in nondel and
                    feature.location.end in nondel):
                feature.location = feature.location + (-offset)
                reduced_annot.features.append(feature)
        end = nondel.end - offset + 2
    SeqIO.write(reduced_annot, gb_m4, "genbank")


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="summarize",
        description=""  # FIXME
    )
    parser.add_argument(
        "-g3", dest="GBM3", required=True,
        help="modified-III GenBank file")
    parser.add_argument(
        "-o", dest="OUT_DIR", required=True,
        help="directory for output files")
    args, unknown = parser.parse_known_args()
    GENBANK_M3 = args.GBM3
    OUT_DIR = args.OUT_DIR
    genbank_id = os.path.splitext(os.path.basename(GENBANK_M3))[0]
    GENBANK_M4 = os.path.join(OUT_DIR, genbank_id + ".gbm4")
    OUT_IMG = os.path.join(OUT_DIR, "genome_reduction")

    # Check input
    # TODO

    save_genbank_m4(GENBANK_M3, GENBANK_M4)

    circplot.plot(
        gb_outer=GENBANK_M3,
        gb_inner=GENBANK_M4,
        out_file=OUT_IMG,
        out_fmt="png"
    )

    # TODO: remove reduced gb
