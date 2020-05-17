#!/usr/bin/env python3
"""circplot.py
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
GC_H = 60
SPACE = 80
DELETION_H = 50
DELETION_REQ = (lambda x: "note" in x.qualifiers and
                          "deletion" in x.qualifiers["note"][0])


def plot(gb_outer, gb_inner, out_file, out_fmt):
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
    genome.plot_ticks(bottom=y + 6, height=20, space=100000)
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
        requirement=DELETION_REQ,
        color=COLORS["ner"]
    )
    # TODO: más features que mostrar
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
