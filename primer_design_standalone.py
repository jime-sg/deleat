#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

import primer3
from Bio import SeqIO
from Bio.Restriction import BamHI
from Bio.SeqRecord import SeqRecord
from regions import Region, Primer, PrimerSet


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
        "PRIMER_NUM_RETURN": 10,
        "PRIMER_FIRST_BASE_INDEX": 1
    }

    results = primer3.bindings.designPrimers(p3_seqargs, p3_globalargs)

    return results


def write_primer_pairs(primer_dict):
    results_str = "N\tsize\tF_start-F_end..R_start-R_end\tseq_F\tseq_R\tTM_F (ºC)\tTM_R (ºC)\n"
    n_pairs = len(primer_dict)//2
    for n in range(n_pairs):
        lp = primer_dict["LEFT_%d" % n]
        rp = primer_dict["RIGHT_%d" % n]
        prod_size = rp.e() - lp.s() + 1

        results_str += "%d\t%d\t%d-%d..%d-%d\t%-24s\t%-24s\t%.1f\t%.1f\n" \
                       % (n+1, prod_size, lp.s(), lp.e(), rp.s(), rp.e(), lp.seq(), rp.seq(), lp.tm(), rp.tm())

    return results_str


def save_pcr_regions(regions_dict, dir):
    with open(dir + "PCR_regions.fna", "w") as f:
        pcr1 = regions_dict["PCR1"]
        pcr1_r = SeqRecord(pcr1.subseq(), id="PCR1", description="%d:%d" % (pcr1.s(), pcr1.e()))
        pcr2 = regions_dict["PCR2"]
        pcr2_r = SeqRecord(pcr2.subseq(), id="PCR2", description="%d:%d" % (pcr2.s(), pcr2.e()))
        prod = pcr1.subseq() + pcr2.subseq()
        prod_r = SeqRecord(prod, id="total_product", description="")
        SeqIO.write((pcr1_r, pcr2_r, prod_r), f, "fasta")



GENOME = "/home/jimena/Bartonella/NC_005955.fna"
LOG_DIR = "/home/jimena/Bartonella/deletions/deletion2/"
DEL_COORDS = (1307339, 1327913)
MARGIN_SIZE = 850

sep = "-" * 80 + "\n"


if __name__ == "__main__":
    genome = SeqIO.read(GENOME, "fasta")
    log = open(LOG_DIR + "primer_design.txt", "w")


    # define initial regions
    del_region = Region(DEL_COORDS, genome)
    del_plusmargins_coords = (del_region.s() - MARGIN_SIZE,
                              del_region.e() + MARGIN_SIZE)
    del_plusmargins = Region(del_plusmargins_coords, genome)
    margin1_coords = (del_plusmargins.s(),
                      del_region.s())
    margin1 = Region(margin1_coords, genome)
    margin2_coords = (del_region.e(),
                      del_plusmargins.e())
    margin2 = Region(margin2_coords, genome)


    # check margin region quality
    log.write(sep + "Checking margins...\n")
    log.write("Initial locations:\n\tMargin 1 (left):  %d - %d\n\tMargin 2 (right): %d - %d\n"
              % (margin1.s(), margin1.e(), margin2.s(), margin2.e()))
    check_BamHItargets_and_repeats(margin1, "margin 1", dir="right")
    check_BamHItargets_and_repeats(margin2, "margin 2", dir="left")
    log.write("\nMargin quality check: done.")
    log.write("\nFinal locations:\n\tMargin 1 (left):  %d - %d\n\tMargin 2 (right): %d - %d\n"
              % (margin1.s(), margin1.e(), margin2.s(), margin2.e()))


    # primer design
    log.write(sep + "Designing primers...\n")
    all_primers = {1: {}, 2: {}}
    for i, margin in enumerate((margin1, margin2)):
        p3_results = design_primers(margin.s(), margin.e())
        n_pairs = p3_results["PRIMER_PAIR_NUM_RETURNED"]
        for j in range(n_pairs):
            for k in ("LEFT", "RIGHT"):
                primer = Primer(p3_results, j, k)
                all_primers[i+1][primer.id] = primer

    log.write("\nBest primer pairs for margin 1 (left):\n")
    log.write(write_primer_pairs(all_primers[1]))
    log.write("\nBest primer pairs for margin 2 (right):\n")
    log.write(write_primer_pairs(all_primers[2]))

    chosen_ids = ["LEFT_0", "RIGHT_0", "LEFT_0", "RIGHT_0"]
    chosen_primers = {
        "1F": all_primers[1][chosen_ids[0]],
        "1R": all_primers[1][chosen_ids[1]],
        "2F": all_primers[2][chosen_ids[2]],
        "2R": all_primers[2][chosen_ids[3]]
    }

    megapriming = PrimerSet(chosen_primers)
    log.write("\nSelected pairs of primers:\n")
    for name, primer in megapriming.primers_raw_dict.items():
        log.write("%s: %s\n" % (name, primer.seq()))

    megapriming.add_tails()
    log.write("\nAdded tails:\n")
    for name, primer in megapriming.primers_tailed_dict.items():
        log.write("%s: %s\n" % (name, primer))

    pcr_regions = megapriming.get_PCR_regions(genome)
    log.write("\nDetermined PCR regions for megapriming:\n")
    log.write("\tPCR1 (left):  %d - %d\n\tPCR2 (right): %d - %d\n"
              % (pcr_regions["PCR1"].s(), pcr_regions["PCR1"].e(),
                 pcr_regions["PCR2"].s(), pcr_regions["PCR2"].e()))

    save_pcr_regions(pcr_regions, LOG_DIR)
    log.write("\nPCR region sequences saved at %s.\n" % LOG_DIR)


    # PCR1/2 region design
    # check for Bam targets


    log.write(sep)
    log.close()
