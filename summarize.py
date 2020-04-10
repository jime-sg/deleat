#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

import circplot


def is_neregion(feature):  # FIXME
    if (feature.type == "misc_feature" and
            "note" in feature.qualifiers and
            feature.qualifiers["note"][0] == "non-essential region"):
        return True
    else:
        return False


def complementary_compoundloc(start, end, comploc):
    complementary = FeatureLocation(start, comploc.parts[0].start)
    for i in range(len(comploc.parts) - 1):
        complementary += FeatureLocation(
            comploc.parts[i].end + 1,
            comploc.parts[i+1].start - 1
        )
    else:
        complementary += FeatureLocation(comploc.parts[-1].end, end)
    return complementary


def reduced_genbank(gb_in, gb_out):
    original = SeqIO.read(gb_in, "genbank")
    deletions = []
    for feature in original.features:
        if is_neregion(feature):
            deletions.append(feature.location)
    deletions = CompoundLocation(deletions)
    non_deletions = complementary_compoundloc(0, len(original), deletions)
    reduced_annot = SeqRecord(
        seq=non_deletions.extract(original.seq),
        id=original.id,
        name=original.name,
        description=original.description,
        dbxrefs=original.dbxrefs,
        annotations=original.annotations
    )
    end = 0
    for nondel in non_deletions.parts:
        offset = nondel.start - end
        for feature in original.features:
            if (feature.location.start in nondel and
                    feature.location.end in nondel):
                feature.location = feature.location + (-offset)
                reduced_annot.features.append(feature)
        end = nondel.end - offset + 2
    SeqIO.write(reduced_annot, gb_out, "genbank")


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="summarize",
        description=""  # FIXME
    )
    MODIFIED_GB = "/home/jimena/Bartonella/NC_005955_wregions.gb"  # FIXME
    REDUCED_GB = "/home/jimena/Escritorio/NC_005955_wregions_reduced.gb"  # FIXME

    reduced_genbank(MODIFIED_GB, REDUCED_GB)

    circplot.plot(
        gb_outer=MODIFIED_GB,
        gb_inner=REDUCED_GB,
        out_file="/home/jimena/Dropbox/TFM/paso_circos/test4",
        out_fmt="png"
    )

    # TODO: remove reduced gb
