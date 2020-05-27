#!/usr/bin/env python3
"""train_classifier.py
# TODO
@author: Jimena Solana
"""

from sys import argv
import os

from nonessential_genes import (find_ori_ter, get_feature_table, extract_cds,
                                strand)
import geptop
import codonw


ORGANISMS = "/home/jimena/Bartonella/DEGdb/organisms2.txt"


if __name__ == "__main__":
    RESULTS_DIR = ""  # TODO
    GENBANK = argv[1]
    OUT_DIR = os.path.join(RESULTS_DIR, DEG_ID)
    os.makedirs(OUT_DIR, exist_ok=True)
    DEG = "/home/jimena/Bartonella/DEGdb/deg_byorg/all"
    CV = "/home/jimena/Bartonella/DEGdb/cv"
    GEPTOP_CUTOFF = 0.24
    NPROC = 2

    with open(ORGANISMS, "r") as fi:
        organisms = {
            org.strip().split("\t")[-1]: org.strip().split("\t")[0]
            for org in fi if not org.startswith("#")
        }

    geptop_params = {
        "deg_path": DEG,  # hay que hacer un directorio auxiliar con todos los organismos de deg menos el que se está analizando
        "cv_path": CV,
        "cutoff": GEPTOP_CUTOFF,
        "n_proc": NPROC,
        "out_path": OUT_DIR
    }

