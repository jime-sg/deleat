#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

import primer3
from Bio import SeqIO
from Bio.Restriction import BamHI
from regions import Region



def BamHI_cuts(seq):
    cuts = False
    if BamHI.search(seq):
        cuts = True

    return cuts


GENOME = "/home/jimena/Bartonella/NC_005955.fna"
DEL_COORDS = (1307339, 1327913 + 1)

genome = SeqIO.read(GENOME, "fasta")

# desired deletion
del_region = Region(genome, DEL_COORDS)

# deletion + margins (PCR1/2)
del_plusmargins_coords = (del_region.s() - 850, del_region.e() + 850)
del_plusmargins = Region(genome, del_plusmargins_coords)

# PCR1
margin1_coords = (del_plusmargins.s(), del_region.s())
margin1 = Region(genome, margin1_coords)

# PCR2
margin2_coords = (del_region.e(), del_plusmargins.e())
margin2 = Region(genome, margin2_coords)


margins_ok = False
Bam_ok = False
repeat_ok = False
while not margins_ok:
    while not Bam_ok
        margin1_Bam_ok = False
        margin2_Bam_ok = False
        if BamHI_cuts(margin1.subseq()):
            displ1 = BamHI.search(margin1.subseq())
            s1 = margin1.s() + displ1
            e1 = margin1.e() + displ1
            margin1 = Region(genome, (s1, e1))
        else:
            margin1_Bam_ok = True
        if BamHI_cuts(margin2.subseq()):
            displ2 = len(margin2.subseq()) - BamHI.search(margin2.subseq())
            s2 = margin1.s() - displ2
            e2 = margin1.e() - displ2
            margin2 = Region(genome, (s2, e2))
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
