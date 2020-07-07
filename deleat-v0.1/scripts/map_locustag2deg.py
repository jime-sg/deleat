#!/usr/bin/env python3
"""map_locustag2deg.py

    Map locus_tag protein ids to their corresponding DEG ids by BLAST
    search.

@author: Jimena Solana
"""

import os
import shutil

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import pandas as pd

from ..geptop import make_blastdb, remove_blastdb

DEG = "/home/jimena/Bartonella/DEGdb/deg_byorg/all"
ORGANISMS = "/home/jimena/Bartonella/DEGdb/organisms.txt"
FASTA_LOCUSTAG_DIR = "/home/jimena/Escritorio/pruebas/features"
NPROC = 4
OUT_FILE = "/home/jimena/Escritorio/pruebas/locustag2deg.csv"


def blast(query, subject, out_dir, n_proc):
    query_fa = os.path.join(out_dir, "temp.fa")
    SeqIO.write(query, query_fa, "fasta")
    results_file = os.path.join(out_dir, "blast_results.xml")
    command = NcbiblastpCommandline(
        query=query_fa, db=subject,
        evalue=1, outfmt=5, num_threads=n_proc,
        out=results_file
    )
    out, err = command()
    os.remove(query_fa)
    return results_file


def same(query, hit):
    t = 0.95  # threshold
    L = len(query)
    if hit.align_length / L > t and hit.identities / L > t:
        return True
    else:
        return False


if __name__ == "__main__":
    with open(ORGANISMS) as f:
        orgs = [org.strip().split("\t")[0]
                for org in f if not org.startswith("#")]

    matches = []
    for deg_id in orgs:
        print(deg_id)
        fasta_locustag = os.path.join(FASTA_LOCUSTAG_DIR,
                                      deg_id, "proteome.faa")
        fasta_deg = os.path.join(DEG, deg_id + ".faa")
        if deg_id.endswith(("a", "b")):
            fasta_deg = os.path.join(DEG, deg_id[:-1] + ".faa")

        make_blastdb(fasta_deg)
        blast_dir = os.path.join(FASTA_LOCUSTAG_DIR, deg_id, "blast")
        os.makedirs(blast_dir, exist_ok=True)

        for protein in SeqIO.parse(fasta_locustag, "fasta"):
            results = blast(protein, fasta_deg, blast_dir, NPROC)
            with open(results) as f:
                record = NCBIXML.read(f)
            if record.alignments:
                hsp = record.alignments[0].hsps[0]
                if same(protein, hsp):
                    hit_name = record.alignments[0].hit_def
                    hit_deg_id = hit_name.split()[0].split("|")[3]
                    pair = (protein.id, hit_deg_id)
                    matches.append(pair)
                else:
                    matches.append((protein.id, "-"))
            else:
                matches.append((protein.id, "-"))

        shutil.rmtree(blast_dir)
        remove_blastdb(fasta_deg)

    map_results = pd.DataFrame(matches, columns=['locus_tag', 'deg_id'])
    map_results.to_csv(OUT_FILE)
