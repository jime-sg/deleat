#!/usr/bin/env python3
"""design_all_primers.py
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser
import os
import subprocess

from Bio import SeqIO

from summarize import is_deletion


if __name__ == "__main__":
    # Parse command-line arguments and init constants
    parser = ArgumentParser(
        prog="design-all-primers",
        description="Design primers for large genome deletions by megapriming."
    )
    parser.add_argument(
        "-g3", dest="GENBANK_M3", required=True,
        help="modified-III GenBank file")
    parser.add_argument(
        "-o", dest="OUT_DIR", required=True,
        help="directory for output files")
    parser.add_argument(
        "-e", dest="ENZYME", required=True,
        help="restriction enzyme used in the experiment (must cut the "
             "megapriming product at both ends only)")
    parser.add_argument(
        "-L", dest="HR_LENGTH", metavar="HR_LENGTH", type=int, default=20,
        choices=range(15, 101),
        help=("min substrate length required for homologous recombination "
              "events (optional, default %(default)s bp)"))
    args, unknown = parser.parse_known_args()
    GENBANK_M3 = args.GENBANK_M3
    OUT_DIR = args.OUT_DIR
    os.makedirs(OUT_DIR, exist_ok=True)
    HR_LENGTH = args.HR_LENGTH
    ENZYME = args.ENZYME

    # Save genome sequence in FASTA format
    genbank_m3 = SeqIO.read(GENBANK_M3, "genbank")
    fasta = os.path.join(OUT_DIR, "genome.fna")
    SeqIO.write(genbank_m3, fasta, "fasta")

    # Design primers
    for feature in genbank_m3.features:
        if is_deletion(feature):
            name = feature.qualifiers["note"][0].split()[-1]
            start = feature.location.start
            end = feature.location.end
            print(
                "Designing primers for deletion %s (%s-%s)..."
                % (name, start, end)
            )
            subprocess.call(
                ["python", "design_primers.py",
                 "-g", fasta, "-o", OUT_DIR,
                 "-n", name, "-d1", str(start), "-d2", str(end),
                 "-e", ENZYME, "-L", str(HR_LENGTH)]
            )

    # Delete temporary FASTA file
    os.remove(fasta)
