#!/usr/bin/env python3
"""extract_features.py
# TODO
@author: Jimena Solana
"""

from sys import argv
import os
import shutil

from Bio import SeqIO

from predict_essentiality import (find_ori_ter, get_feature_table,
                                  CODONW_FEATURES)
from geptop import ORG_NAMES

DEG = "/home/jimena/Bartonella/DEGdb/deg_byorg/all"  # FIXME
CV = "/home/jimena/Bartonella/DEGdb/cv"  # FIXME
RESULTS_DIR = "/home/jimena/Bartonella/ge_classifier/features"  # FIXME


if __name__ == "__main__":
    GENBANK = argv[1]
    NPROC = argv[2]
    DEG_ID = os.path.splitext(os.path.basename(GENBANK))[0]
    OUT_DIR = os.path.join(RESULTS_DIR, DEG_ID)
    os.makedirs(OUT_DIR, exist_ok=True)

    if DEG_ID.endswith(("a", "b")):
        print("Extracting features for organism:", DEG_ID, ORG_NAMES[DEG_ID[:-1]])
    else:
        print("Extracting features for organism:", DEG_ID, ORG_NAMES[DEG_ID])

    annot = SeqIO.read(GENBANK, "genbank")
    ori, ter = find_ori_ter(annot)

    # create dir with all DEG organisms except the one being analysed
    deg2 = os.path.join(RESULTS_DIR, "temp")
    shutil.copytree(DEG, deg2)
    if DEG_ID.endswith(("a", "b")):
        os.remove(os.path.join(deg2, DEG_ID[:-1] + ".faa"))
    else:
        os.remove(os.path.join(deg2, DEG_ID + ".faa"))

    geptop_params = {
        "deg_path": deg2,
        "cv_path": CV,
        "n_proc": NPROC,
        "out_path": OUT_DIR
    }

    results = get_feature_table(annot, OUT_DIR, ori, ter,
                                geptop_params, CODONW_FEATURES)
    results.to_csv(os.path.join(RESULTS_DIR, DEG_ID + "_feature_table.csv"))

    shutil.rmtree(deg2)
