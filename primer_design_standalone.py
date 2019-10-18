#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

import primer3
from Bio import SeqIO
from Bio.Restriction import BamHI
from regions import Region


def check_BamHItargets_and_repeats(seq, seq_name, dir):
    seq_ok = False
    log.write("\nChecking %s...\n" % seq_name)
    while not seq_ok:
        bam_ok = False
        repeat_ok = True  # TEMP

        # check seq for BamHI target sites
        while not bam_ok:
            pos = BamHI.search(seq.subseq())
            if pos:
                pos = pos[0]
                globalpos = pos + seq.s()
                seq.displace_past(pos, dir)
                log.write("!! BamHI target site found at %d. Location reset to %d - %d.\n"
                          % (globalpos, seq.s(), seq.e()))
            else:
                bam_ok = True

        # check seq for repetitive regions
        # while not repeat_ok:
        # ...
        # if repeats:
        # seq.displace_past(x, dir)
        # log.write
        # else:
        # repeat_ok = True

        if bam_ok and repeat_ok:
            seq_ok = True


GENOME = "/home/jimena/Bartonella/NC_005955.fna"
LOG_DIR = "/home/jimena/Bartonella/deletions/deletion2/"
DEL_COORDS = (1459510, 1467913 + 1)
MARGIN_SIZE = 850

genome = SeqIO.read(GENOME, "fasta")
log = open(LOG_DIR + "primer_design.txt", "w")


# define initial regions
## desired deletion
del_region = Region(genome, DEL_COORDS)

## deletion + margins (for PCR1/2)
del_plusmargins_coords = (del_region.s() - MARGIN_SIZE,
                          del_region.e() + MARGIN_SIZE)
del_plusmargins = Region(genome, del_plusmargins_coords)

## margin 1
margin1_coords = (del_plusmargins.s(),
                  del_region.s())
margin1 = Region(genome, margin1_coords)

## margin 2
margin2_coords = (del_region.e(),
                  del_plusmargins.e())
margin2 = Region(genome, margin2_coords)

# check margin region quality
log.write("-"*30 + "\nChecking margins...\n")
log.write("Initial locations:\n\tMargin 1 (left):  %d - %d\n\tMargin 2 (right): %d - %d\n"
          %(margin1.s(), margin1.e(), margin2.s(), margin2.e()))
check_BamHItargets_and_repeats(margin1, "margin 1", dir="right")
check_BamHItargets_and_repeats(margin2, "margin 2", dir="left")
log.write("\nMargin quality check: done.")
log.write("\nFinal locations:\n\tMargin 1 (left):  %d - %d\n\tMargin 2 (right): %d - %d\n"
      %(margin1.s(), margin1.e(), margin2.s(), margin2.e()))
log.write("-"*30 + "\n")

#log.close()
#exit()


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


log.close()