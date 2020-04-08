#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from Bio import SeqIO

import pycircos


COLORS = {
    "circle": "#7c00ff",
    "pos": "#a9ff8a",
    "neg": "#ff8a8a",
    "gp": "#9d78ff",
    "gm": "#ffd94d",
    "rt": "#000000",
    "rr": "#517c8c",
    "ner": "#eb0000"
}
GENE_H = 60
RNA_H = 20
SPACE = 70
NEREGION_H = 50
# NEREGION_REQ = (lambda x: "note" in x.qualifiers and
#                           x.qualifiers["note"][0] == "non-essential region")  # FIXME
NEREGION_REQ = (lambda x: "colour" in x.qualifiers)  # FIXME


def plot(genbank, out_file, out_fmt):
    y = 900
    genome = pycircos.Genome()
    # Outer circle
    h = 2
    y += h
    genome.read_locus(
        SeqIO.parse(genbank, "genbank"),
        interspace=0, bottom=y+h, height=h,
        color_list=[COLORS["circle"]],
    )

    # Ticks
    genome.plot_ticks(bottom=y+6, height=20, space=100000)

    # + strand genes
    h = GENE_H
    y -= h
    genome.plot_feature(
        feat_type="gene",
        bottom=y, height=h,
        requirement=lambda x: x.strand == 1,
        color=COLORS["gp"]
    )
    # - strand genes
    y -= h
    genome.plot_feature(
        feat_type="gene",
        bottom=y, height=h,
        requirement=lambda x: x.strand == -1,
        color=COLORS["gm"]
    )
    # GC skew
    h = GENE_H
    chrom = list(genome.locus_dict.keys())[0]
    genome.calc_gcskew(chrom, window_size=1000, slide_size=1000)
    genome.plot_data(
        chrom, pycircos.np.array(genome.locus_dict[chrom]["gc_skew"]),
        bottom=y+h, height=h,
        xaxes=True, plot_style="fill",
        color="#e8e8e8", color1=COLORS["pos"], color2=COLORS["neg"]
    )
    # RNA genes
    h = RNA_H
    y -= h + 5
    genome.plot_feature(
        feat_type="tRNA",
        bottom=y, height=h,
        color=COLORS["rt"]
    )
    genome.plot_feature(
        feat_type="rRNA",
        bottom=y, height=h,
        color=COLORS["rr"]
    )

    # Non-essential regions
    h = NEREGION_H
    y -= h + 40
    genome.plot_feature(
        feat_type="misc_feature",
        bottom=y, height=h,
        requirement=NEREGION_REQ,
        color=COLORS["ner"]
    )

    # + strand genes
    h = GENE_H
    y -= h + SPACE
    genome.plot_feature(
        feat_type="gene",
        bottom=y, height=h,
        requirement=lambda x: x.strand == 1,
        color=COLORS["gp"]
    )
    genome.plot_feature(
        feat_type="CDS",
        bottom=y, height=h,
        requirement=lambda x: x.strand == 1 and "pseudo" in x.qualifiers,
        color=COLORS["gp"]
    )
    # - strand genes
    y -= h
    genome.plot_feature(
        feat_type="gene",
        bottom=y, height=h,
        requirement=lambda x: x.strand == -1,
        color=COLORS["gm"]
    )
    genome.plot_feature(
        feat_type="CDS",
        bottom=y, height=h,
        requirement=lambda x: x.strand == -1 and "pseudo" in x.qualifiers,
        color=COLORS["gm"]
    )
    # RNA genes
    h = RNA_H
    y -= h + 5
    genome.plot_feature(
        feat_type="tRNA",
        bottom=y, height=h,
        color=COLORS["rt"]
    )
    genome.plot_feature(
        feat_type="rRNA",
        bottom=y, height=h,
        color=COLORS["rr"]
    )
    # TODO: features que mostrar en el gráfico
    genome.save(file_name=out_file, format_type=out_fmt)
