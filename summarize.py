#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser

from Bio import SeqIO

from pycircos import Genome


def circular_plot(genbank, out_file, out_fmt):
    genome = Genome()
    # Outer circle
    genome.read_locus(
        SeqIO.parse(genbank, "genbank"),
        interspace=0, bottom=900, height=10,
        color_list=["#7c00ff"],
    )
    chrom = list(genome.locus_dict.keys())[0]
    # Leading strand genes
    genome.plot_feature(
        feat_type="gene",
        bottom=740,
        requirement=lambda x: x.strand == 1,
        color="#74AFB9"
    )
    # Lagging strand genes
    genome.plot_feature(
        feat_type="gene",
        bottom=820,
        requirement=lambda x: x.strand == -1,
        color="#C481A1"
    )
    # TODO: features que mostrar en el gráfico
    genome.save(file_name=out_file, format_type=out_fmt)


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="summarize",
        description=""  # FIXME
    )
