#!/usr/bin/env python3
"""circplot.py

    Create a circular genome plot comparing a genome and its reduced
    version.

@author: Jimena Solana
"""

from Bio import SeqIO

import pycircos
from design_all_primers import is_deletion


COLORS = {
    "circle": "#7c00ff",
    "pos": "#c9ebff",
    "neg": "#fff08c",
    "gp": "#b174f7",
    "gm": "#91f76f",
    "rt": "#000000",
    "rr": "#517c8c",
    "ner": "#eb0000"
}
GENE_H = 60
RNA_H = 20
GC_H = 60
SPACE = 80
DELETION_H = 50


def plot(gb_outer, gb_inner, out_file, out_fmt):
    """Plot the original and reduced genomes.

    Args:
        gb_outer (str): modified-III GenBank file path.
        gb_inner (str): modified-IV GenBank file path.
        out_file (str): output file path.
        out_fmt (str): output file format.
    """
    y = 900
    genome = pycircos.Genome()
    print("Drawing original genome...")
    # OUTER: circle
    h = 2
    y += h
    genome.read_locus(
        SeqIO.parse(gb_outer, "genbank"),
        interspace=0, bottom=y + 4, height=h,
        color_list=[COLORS["circle"]]
    )
    # OUTER: ticks
    genome.plot_ticks(bottom=y + 6, height=20, space=100000, labels=True)
    # + strand genes
    print("  Genes...")
    h = GENE_H
    y -= h
    genome.plot_feature(
        feat_type="gene",
        bottom=y, height=h,
        requirement=lambda x: x.strand == 1,
        color=COLORS["gp"]
    )
    # OUTER: - strand genes
    y -= h
    genome.plot_feature(
        feat_type="gene",
        bottom=y, height=h,
        requirement=lambda x: x.strand == -1,
        color=COLORS["gm"]
    )
    # OUTER: RNA genes
    h = RNA_H
    y -= h + 5
    genome.plot_feature(
        feat_type="tRNA",
        bottom=y, height=h, thicken=5,
        color=COLORS["rt"]
    )
    genome.plot_feature(
        feat_type="rRNA",
        bottom=y, height=h,
        color=COLORS["rr"]
    )
    # OUTER: GC skew
    print("  GC skew graph...")
    h = GC_H
    y -= 2 * h
    chrom = list(genome.locus_dict.keys())[0]
    genome.calc_gcskew(chrom, window_size=1000, slide_size=1000)
    genome.plot_data(
        chrom, pycircos.np.array(genome.locus_dict[chrom]["gc_skew"]),
        bottom=y + h, height=h,
        xaxes=True, plot_style="fill",
        color="#e8e8e8", color1=COLORS["pos"], color2=COLORS["neg"]
    )
    # OUTER: deletions
    h = DELETION_H
    y -= h
    genome.plot_feature(
        feat_type="misc_feature",
        bottom=y, height=h,
        requirement=is_deletion,
        color=COLORS["ner"]
    )
    print("Drawing reduced genome...")
    # INNER: circle
    h = 2
    y -= h + SPACE
    genome.read_locus(
        SeqIO.parse(gb_inner, "genbank"),
        interspace=0, bottom=y + 4, height=h,
        color_list=[COLORS["circle"]],
        reset=True
    )
    # INNER: ticks
    genome.plot_ticks(bottom=y + 6, height=20, space=100000)
    # INNER: + strand genes
    print("  Genes...")
    h = GENE_H
    y -= h
    genome.plot_feature(
        feat_type="gene",
        bottom=y, height=h,
        requirement=lambda x: x.strand == 1,
        color=COLORS["gp"]
    )
    # INNER: - strand genes
    y -= h
    genome.plot_feature(
        feat_type="gene",
        bottom=y, height=h,
        requirement=lambda x: x.strand == -1,
        color=COLORS["gm"]
    )
    # INNER: RNA genes
    h = RNA_H
    y -= h + 5
    genome.plot_feature(
        feat_type="tRNA",
        bottom=y, height=h, thicken=5,
        color=COLORS["rt"]
    )
    genome.plot_feature(
        feat_type="rRNA",
        bottom=y, height=h,
        color=COLORS["rr"]
    )
    # OUTER: GC skew
    print("  GC skew graph...")
    h = GC_H
    y -= 2 * h
    chrom = list(genome.locus_dict.keys())[0]
    genome.calc_gcskew(chrom, window_size=1000, slide_size=1000)
    genome.plot_data(
        chrom, pycircos.np.array(genome.locus_dict[chrom]["gc_skew"]),
        bottom=y + h, height=h,
        xaxes=True, plot_style="fill",
        color="#e8e8e8", color1=COLORS["pos"], color2=COLORS["neg"]
    )
    print("Saving image...")
    genome.save(file_name=out_file, format_type=out_fmt)
    print("Done.")
