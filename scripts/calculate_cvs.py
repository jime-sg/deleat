#!/usr/bin/env python3
"""
@author: Jimena Solana
"""

import os
import json

from ..geptop import composition_vector, CV_K


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
