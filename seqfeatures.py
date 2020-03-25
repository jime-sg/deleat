#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from Bio import SeqIO
from Bio.SeqUtils import GC as gc_content


def run(gb_file):
    annotation = SeqIO.read(gb_file, "genbank")
    genomeseq = annotation.seq
    for gene in annotation.features:
        if gene.type == "CDS":
            seqfeatures = get_features(gene, genomeseq)
            seqfeatures = ";".join(
                str(k) + "," + str(v) for k, v in seqfeatures.items()
            )
            gene.qualifiers["seqfeatures"] = seqfeatures
    return annotation


def get_features(gene, genome):
    features = {
        "length": len(gene),
        "strand": gene.strand,
        "gc": gc_content(gene.extract(genome))/100
    }
    return features


if __name__ == "__main__":
    # gb = run("/home/jimena/Bartonella/NC_005955.gb")
    # SeqIO.write(gb, "/home/jimena/Escritorio/newannot.gbm", "genbank")
    pass
