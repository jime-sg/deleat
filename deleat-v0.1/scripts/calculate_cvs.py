#!/usr/bin/env python3
"""calculate_cvs.py

    Calculate proteome composition vector for all DEG reference
    organisms, according to the algorithm by Qi et al. 2004 'Whole
    Proteome Prokaryote Phylogeny Without Sequence Alignment: A K-String
    Composition Approach'.

@author: Jimena Solana
"""

import os
import json

from ..geptop import composition_vector


if __name__ == "__main__":
    REFSEQS_DIR = "/home/jimena/Bartonella/DEGdb/deg_byorg/all"
    OUT_DIR = "/home/jimena/Bartonella/DEGdb/cv"

    os.makedirs(OUT_DIR, exist_ok=True)
    reference_genomes = os.listdir(REFSEQS_DIR)

    for genome in reference_genomes:
        id_ = os.path.splitext(genome)[0]
        print(id_)
        cv = composition_vector(os.path.join(REFSEQS_DIR, genome))
        with open(os.path.join(OUT_DIR, id_ + ".json"), "w") as fo:
            json.dump(cv, fo)
