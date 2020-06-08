#!/usr/bin/env python3
"""revise_deletions.py

    Redefine deletion list after manual curation.

@author: Jimena Solana
"""

from argparse import ArgumentParser
import os

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import pandas as pd

from define_deletions import make_table


def is_deletion(feature):
    """Check whether a GenBank feature is an annotated deletion."""
    if (feature.type == "misc_feature" and
            "note" in feature.qualifiers and
            "deletion" in feature.qualifiers["note"][0]):
        return True
    else:
        return False


def parse_deletions(table):
    """Parse pd.DataFrame of deletion coordinates into a
    CompoundLocation.

    Args:
        table (pd.DataFrame): table of deletion coordinates.
    Returns:
        dels_loc (Bio.SeqFeature.CompoundLocation): deletion list.
    """
    dels = []
    for row in range(len(table.index)):
        s = int(table.iloc[row, 0]) - 1  # 1-based > 0-based
        e = int(table.iloc[row, 1])
        dels.append(FeatureLocation(s, e))
    dels_loc = CompoundLocation(dels)
    return dels_loc


def save_genbank_m3(deletions, names, gb_m2, gb_m3):
    """Generate modified-III GenBank file with annotated deletions.

    Args:
        deletions (Bio.SeqFeature.CompoundLocation): list of deletions.
        names (list of str): deletion names.
        gb_m2 (Bio.SeqRecord.SeqRecord): modified-II GenBank annotation.
        gb_m3 (str): modified-III GenBank file path.
    """
    # Delete previous proposed deletions
    features_nodel = []
    for feature in gb_m2.features:
        if not is_deletion(feature):
            features_nodel.append(feature)
    gb_m2.features = features_nodel
    # Add revised deletions
    for n, deletion in enumerate(deletions.parts):
        deletion_feature = SeqFeature(
            location=deletion,
            type="misc_feature"
        )
        deletion_feature.qualifiers["note"] = ["deletion %s" % (names[n])]
        gb_m2.features.append(deletion_feature)
    SeqIO.write(gb_m2, gb_m3, "genbank")


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="revise-deletions",
        description="Redefine deletion list after manual curation."
    )
    parser.add_argument(
        "-g2", dest="GBM2", required=True,
        help="modified-II GenBank file")
    parser.add_argument(
        "-t", dest="TABLE", required=True,
        help="revised table of proposed deletions (csv)")
    parser.add_argument(
        "-o", dest="OUT_DIR", required=True,
        help="directory for output files")
    args, unknown = parser.parse_known_args()
    GENBANK_M2 = args.GBM2
    IN_TABLE = args.TABLE
    OUT_DIR = args.OUT_DIR
    os.makedirs(OUT_DIR, exist_ok=True)
    genbank_id = os.path.splitext(os.path.basename(GENBANK_M2))[0]
    GENBANK_M3 = os.path.join(OUT_DIR, genbank_id + ".gbm3")
    OUT_TABLE = os.path.join(OUT_DIR, "revised_deletions.csv")

    # Check input
    try:
        genbank_m2 = SeqIO.read(GENBANK_M2, "genbank")
    except (FileNotFoundError, ValueError):
        raise SystemExit("\n\terror: could not read .gbm2 file\n")
    contains_dels = False
    for feature in genbank_m2.features:
        if is_deletion(feature):
            contains_dels = True
            break
    if not contains_dels:
        raise SystemExit("\n\terror: invalid GenBank file (must be .gbm2)\n")
    try:
        in_table = pd.read_csv(IN_TABLE,
                               index_col=0, usecols=range(3), comment="#")
    except (FileNotFoundError, ValueError):
        raise SystemExit("\n\terror: could not read deletion table\n")
    if not in_table.values.all() in range(1, len(genbank_m2)+1):
        raise SystemExit("\n\terror: invalid deletion coordinates\n")
    if any(in_table.index.duplicated()):
        duplicates = in_table.index[in_table.index.duplicated()].to_list()
        raise SystemExit(
            "\n\terror: duplicate deletion names: %s\n" % duplicates
        )

    # Redefine deletions
    revised_deletions = parse_deletions(in_table)
    names = in_table.index.to_list()
    save_genbank_m3(revised_deletions, names,  genbank_m2, GENBANK_M3)
    deletions_table = make_table(revised_deletions, genbank_m2)
    deletions_table.index = in_table.index
    deletions_table.to_csv(OUT_TABLE)
