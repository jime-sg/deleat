#!/usr/bin/env python3
"""fetch_deg_genbank.py

    Download GenBank annotation files for each DEG reference organism.

@author: Jimena Solana
"""

import os
import string

from Bio import Entrez


ORGANISMS = "/home/jimena/Bartonella/DEGdb/organisms2.txt"
OUT_DIR = "/home/jimena/Bartonella/DEGdb/genbank"
EMAIL = "jisogon@alumni.uv.es"


def fetch(accesion, out_file):
    with Entrez.efetch(
            db="nuccore",
            id=accesion,
            rettype="gbwithparts",
            retmode="text"
    ) as handle:
        with open(out_file, "w") as fo:
            fo.write(handle.read())


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    Entrez.email = EMAIL

    with open(ORGANISMS, "r") as fi:
        organisms = {
            org.strip().split("\t")[-1]: org.strip().split("\t")[0]
            for org in fi if not org.startswith("#")
        }

    for organism in organisms:
        print("Fetching %s..." % organism)
        if ";" not in organism:
            fetch(organism, os.path.join(OUT_DIR, organisms[organism] + ".gb"))
        else:
            n = 0
            for acc in organism.split(";"):
                suf = string.ascii_lowercase[n]
                fetch(acc, os.path.join(OUT_DIR, organisms[organism] + suf + ".gb"))
                n += 1
