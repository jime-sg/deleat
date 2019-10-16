#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

import primer3
from Bio import SeqIO
from Bio.Restriction import BamHI



def set_subsequence(genome, coords):
    s = coords[0]
    e = coords[1]
    subseq = genome[s:e]

    return s, e, subseq


def BamHI_cuts(seq):
    cuts = False
    if BamHI.search(seq):
        cuts = True

    return cuts


GENOME = "/home/jimena/Bartonella/NC_005955.fna"
DEL_COORDS = (1307339, 1327913 + 1)

genome = SeqIO.read(GENOME, "fasta")

# desired deletion
del_start, del_end, del_seq = set_subsequence(genome, DEL_COORDS)

# deletion + margins (PCR1/2)
total_coords = (del_start - 850, del_end + 850)
total_start, total_end, total_seq = set_subsequence(genome, total_coords)

# PCR1
margin1_coords = (total_start, del_start)
margin1_start, margin1_end, margin1_seq = set_subsequence(genome, margin1_coords)

# PCR2
margin2_coords = (del_end, total_end)
margin2_start, margin2_end, margin2_seq = set_subsequence(genome, margin2_coords)


margins_ok = False
Bam_ok = False
repeat_ok = False
while not margins_ok:
    while not Bam_ok
        margin1_Bam_ok = False
        margin2_Bam_ok = False
        if BamHI_cuts(margin1_seq):
            displ1 = BamHI.search(margin1_seq)
            s1 = margin1_start + displ1
            e1 = margin1_end + displ1
            margin1_start, margin1_end, margin1_seq = set_subsequence(genome, (s1, e1))
        else:
            margin1_Bam_ok = True
        if BamHI_cuts(margin2_seq):
            displ2 = len(margin2_seq) - BamHI.search(margin2_seq)
            s2 = margin1_start - displ2
            e2 = margin1_end - displ2
            margin2_start, margin2_end, margin2_seq = set_subsequence(genome, (s2, e2))
        else:
            margin2_Bam_ok = True

        if margin1_Bam_ok and margin2_Bam_ok:
            Bam_ok = True

    while not repeat_ok:
        # ...
        # if repeats: reset_subsequence, Bam_ok = False, break
        # else: repeat_ok
    

    if Bam_ok and repeat_ok:
        margins_ok = True

p3_seqargs = {
    "SEQUENCE_ID": "deletion_2",
    "SEQUENCE_TEMPLATE": str(p3_seq.seq),
    "SEQUENCE_INCLUDED_REGION": [1, 850]
}

p3_globalargs = {
    "PRIMER_TASK": "generic",
    "PRIMER_PRODUCT_SIZE_RANGE": [750, 850]
}

results = primer3.bindings.designPrimers(p3_seqargs, p3_globalargs)

print("N\tsize\tpos_fragment\tpos_genome\tseq\tTM")
for n in range(0, results["PRIMER_PAIR_NUM_RETURNED"]):
    p_size = results["PRIMER_PAIR_%d_PRODUCT_SIZE" %n]
    print("%s\t%s" %(n+1, p_size))

    for m in ("LEFT", "RIGHT"):
        o = m + "_" + str(n)
        seq = results["PRIMER_%s_SEQUENCE" %o]
        pos_f = (results["PRIMER_%s" %o][0],
                 results["PRIMER_%s" %o][0]+results["PRIMER_%s" %o][1])
        pos_g = (pos_f[0]+p3_start, pos_f[1]+p3_start)
        tm = results["PRIMER_%s_TM" %o]
        print("\t\t%10s\t%s\t%s\t%.1f ºC"
              %(pos_f, pos_g, seq, tm))
