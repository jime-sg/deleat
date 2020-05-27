#!/usr/bin/env python3
"""extract_features.py
# TODO
@author: Jimena Solana
"""

from sys import argv
import os

from Bio import SeqIO

from ..nonessential_genes import (find_ori_ter, get_feature_table,
                                CODONW_FEATURES)

DEG = "/home/jimena/Bartonella/DEGdb/deg_byorg/all"
CV = "/home/jimena/Bartonella/DEGdb/cv"
RESULTS_DIR = ""  # TODO


if __name__ == "__main__":
    GENBANK = argv[1]
    NPROC = argv[2]
    DEG_ID = os.path.splitext(os.path.basename(GENBANK))[0]
    OUT_DIR = os.path.join(RESULTS_DIR, DEG_ID)
    os.makedirs(OUT_DIR, exist_ok=True)

    annot = SeqIO.read(GENBANK, "genbank")
    ori, ter = find_ori_ter(annot.seq)

    geptop_params = {
        "deg_path": DEG,  # hay que hacer un directorio auxiliar con todos los organismos de deg menos el que se está analizando
        "cv_path": CV,
        "n_proc": NPROC,
        "out_path": OUT_DIR
    }

    results = get_feature_table(GENBANK, OUT_DIR, ori, ter,
                                geptop_params, CODONW_FEATURES)
    results.to_csv(os.path.join(RESULTS_DIR, DEG_ID + "_feature_table.csv"))
