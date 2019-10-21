#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

import primer3
from Bio import SeqIO
from Bio.Restriction import BamHI
from Bio.Seq import Seq
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


def design_primers(start, end):
    p3_seqargs = {
        "SEQUENCE_TEMPLATE": str(genome.seq),
        "SEQUENCE_INCLUDED_REGION": [start, end - start]
    }

    p3_globalargs = {
        "PRIMER_TASK": "generic",
        "PRIMER_PRODUCT_SIZE_RANGE": [750, 850],
        "PRIMER_PRODUCT_OPT_SIZE": 800,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": 0.2,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 0.2,
        "PRIMER_NUM_RETURN": 10
    }

    results = primer3.bindings.designPrimers(p3_seqargs, p3_globalargs)

    return results


def parse_p3results(results_dict):
    results_str = "N\tsize\tF_start-F_end..R_start-R_end\tseq_F\tseq_R\tTM_F (ºC)\tTM_R (ºC)\n"
    for n in range(0, results_dict["PRIMER_PAIR_NUM_RETURNED"]):
        prod_size = results_dict["PRIMER_PAIR_%d_PRODUCT_SIZE" % n]
        lp_start = results_dict["PRIMER_LEFT_%d" % n][0]
        lp_end = results_dict["PRIMER_LEFT_%d" % n][0] + results_dict["PRIMER_LEFT_%d" % n][1] - 1
        rp_start = results_dict["PRIMER_RIGHT_%d" % n][0]
        rp_end = results_dict["PRIMER_RIGHT_%d" % n][0] + results_dict["PRIMER_RIGHT_%d" % n][1] - 1
        lp_seq = results_dict["PRIMER_LEFT_%d_SEQUENCE" % n]
        rp_seq = results_dict["PRIMER_RIGHT_%d_SEQUENCE" % n]
        lp_tm = results_dict["PRIMER_LEFT_%d_TM" % n]
        rp_tm = results_dict["PRIMER_RIGHT_%d_TM" % n]

        results_str += "%d\t%d\t%d-%d..%d-%d\t%-24s\t%-24s\t%.1f\t%.1f\n" \
                       % (n + 1, prod_size, lp_start, lp_end, rp_start, rp_end, lp_seq, rp_seq, lp_tm, rp_tm)

    return results_str


GENOME = "/home/jimena/Bartonella/NC_005955.fna"
LOG_DIR = "/home/jimena/Bartonella/deletions/deletion2/"
DEL_COORDS = (1307339, 1327913)
MARGIN_SIZE = 850

sep = "-" * 80 + "\n"

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
log.write(sep + "Checking margins...\n")
log.write("Initial locations:\n\tMargin 1 (left):  %d - %d\n\tMargin 2 (right): %d - %d\n"
          % (margin1.s(), margin1.e(), margin2.s(), margin2.e()))
check_BamHItargets_and_repeats(margin1, "margin 1", dir="right")
check_BamHItargets_and_repeats(margin2, "margin 2", dir="left")
log.write("\nMargin quality check: done.")
log.write("\nFinal locations:\n\tMargin 1 (left):  %d - %d\n\tMargin 2 (right): %d - %d\n"
          % (margin1.s(), margin1.e(), margin2.s(), margin2.e()))

# log.close()
# exit()


# primer design
log.write(sep + "Designing primers...\n")
results_margin1 = design_primers(margin1.s(), margin1.e())
results_margin2 = design_primers(margin2.s(), margin2.e())
log.write("\nBest primer pairs for margin 1 (left):\n")
log.write(parse_p3results(results_margin1))
log.write("\nBest primer pairs for margin 2 (right):\n")
log.write(parse_p3results(results_margin2))

primers_raw = {
    "PCR1F": results_margin1["PRIMER_LEFT_0_SEQUENCE"].upper(),
    "PCR1R": results_margin1["PRIMER_RIGHT_0_SEQUENCE"].upper(),
    "PCR2F": results_margin2["PRIMER_LEFT_0_SEQUENCE"].upper(),
    "PCR2R": results_margin2["PRIMER_RIGHT_0_SEQUENCE"].upper()
}

primers_tailed = {
    "PCR1F": "gcacggatcc" + primers_raw["PCR1F"],
    "PCR1R": primers_raw["PCR1R"],
    "PCR2F": Seq(primers_raw["PCR1R"]).reverse_complement().lower() + primers_raw["PCR2F"],
    "PCR2R": "gcacggatcc" + primers_raw["PCR2R"]
}

log.write("\nSelected first pairs of primers:\n")
for a, b in primers_raw.items():
    log.write("%s: %s\n" % (a, b))

log.write("\nAdded tails:\n")
for a, b in primers_tailed.items():
    log.write("%s: %s\n" % (a, b))

# PCR1/2 region design
# check for Bam targets

log.write(sep)
log.close()
