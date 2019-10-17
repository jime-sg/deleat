#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

import primer3
from Bio import SeqIO
from Bio.Restriction import BamHI
from regions import Region


GENOME = "/home/jimena/Bartonella/NC_005955.fna"
DEL_COORDS = (1307339, 1327913 + 1)
MARGIN_SIZE = 850

genome = SeqIO.read(GENOME, "fasta")


# define initial regions
## desired deletion
del_region = Region(genome, DEL_COORDS)

## deletion + margins (PCR1/2)
del_plusmargins_coords = (del_region.s() - MARGIN_SIZE,
                          del_region.e() + MARGIN_SIZE)
del_plusmargins = Region(genome, del_plusmargins_coords)

## PCR1
margin1_coords = (del_plusmargins.s(),
                  del_region.s())
margin1 = Region(genome, margin1_coords)

## PCR2
margin2_coords = (del_region.e(),
                  del_plusmargins.e())
margin2 = Region(genome, margin2_coords)


# check PCR1/2 region quality
margin1_ok = False
margin2_ok = False

## check PCR1
while not margin1_ok:
    margin1_Bam_ok = False
    margin1_repeat_ok = False
    # check PCR1 for BamHI target sites
    while not margin1_Bam_ok:
        pos1 = BamHI.search(margin1.subseq())
        if pos1:
            margin1.displace_past(pos1, "right")
        else:
            margin1_Bam_ok = True

    # check PCR1 for repetitive regions
    # while not margin1_repeat_ok:
    # ...
    # if repeats:
    # margin1.displace_past(x, "right")
    # else:
    # margin1_repeat_ok = True

    # if margin1_Bam_ok and margin1_repeat_ok:
    #     margin1_ok = True

## check PCR2
while not margin2_ok:
    margin2_Bam_ok = False
    margin2_repeat_ok = False
    # check PCR2 for BamHI target sites
    while not margin2_Bam_ok:
        pos2 = BamHI.search(margin2.subseq())
        if pos2:
            margin2.displace_past(pos2, "left")
        else:
            margin2_Bam_ok = True

    # check PCR2 for repetitive regions
    # while not margin2_repeat_ok:
    # ...
    # if repeats:
    # margin2.displace_past(x, "left")
    # else:
    # margin2_repeat_ok = True

    # if margin2_Bam_ok and margin1_repeat_ok:
    #     margin2_ok = True




# primer design
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
