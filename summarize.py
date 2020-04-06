#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser

from Bio import SeqIO

from pycircos import Genome


NEREGION_REQ = (lambda x: "note" in x.qualifiers and
                          x.qualifiers["note"][0] == "non-essential region")  # FIXME


def circular_plot(genbank, out_file, out_fmt):
    genome = Genome()
    # Outer circle
    genome.read_locus(
        SeqIO.parse(genbank, "genbank"),
        interspace=0, bottom=904, height=3,
        color_list=["#7c00ff"],
    )
    chrom = list(genome.locus_dict.keys())[0]

    # Ticks
    genome.plot_ticks(bottom=909, space=100000)

    # + strand genes
    genome.plot_feature(
        feat_type="gene",
        bottom=850, height=50,
        requirement=lambda x: x.strand == 1,
        color="#74AFB9"
    )
    # - strand genes
    genome.plot_feature(
        feat_type="gene",
        bottom=800, height=50,
        requirement=lambda x: x.strand == -1,
        color="#C481A1"
    )

    # Non-essential regions
    genome.plot_feature(
        feat_type="misc_feature",
        bottom=720, height=40,
        requirement=NEREGION_REQ,
        color="#eb0000"
    )
    # TODO: features que mostrar en el gráfico
    genome.save(file_name=out_file, format_type=out_fmt)


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="summarize",
        description=""  # FIXME
    )

    circular_plot(
        genbank="/home/jimena/Bartonella/NC_005955.gb",
        out_file="/home/jimena/Dropbox/TFM/paso_circos/test4",
        out_fmt="png"
    )
