#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

from Bio import SeqIO
from Bio.Restriction import BamHI
from regions import Region
from primers import PrimerSet, design_primers, write_primer_pairs, save_pcr_regions


def check_BamHItargets_and_repeats(seq, seq_name, direction):
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
                seq.displace_past(pos, direction)
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
DEL_COORDS = (1307339, 1327913)
MARGIN_SIZE = 850

sep = "-" * 80 + "\n"


if __name__ == "__main__":
    genome = SeqIO.read(GENOME, "fasta")
    log = open(LOG_DIR + "primer_design.txt", "w")

    # define initial regions
    del_region = Region(DEL_COORDS, genome)
    margin1_coords = (del_region.s() - MARGIN_SIZE,
                      del_region.s())
    margin1 = Region(margin1_coords, genome)
    margin2_coords = (del_region.e(),
                      del_region.e() + MARGIN_SIZE)
    margin2 = Region(margin2_coords, genome)

    # check margin region quality
    log.write(sep + "Checking margins...\n")
    log.write("Initial locations:\n\tMargin 1 (left):  %d - %d\n\tMargin 2 (right): %d - %d\n"
              % (margin1.s(), margin1.e(), margin2.s(), margin2.e()))
    check_BamHItargets_and_repeats(margin1, "margin 1", direction="right")
    check_BamHItargets_and_repeats(margin2, "margin 2", direction="left")
    log.write("\nMargin quality check: done.")
    log.write("\nFinal locations:\n\tMargin 1 (left):  %d - %d\n\tMargin 2 (right): %d - %d\n"
              % (margin1.s(), margin1.e(), margin2.s(), margin2.e()))

    # primer design
    log.write(sep + "Designing primers...\n")
    all_primers = {}
    all_primers[1] = design_primers(margin1, genome)
    log.write("\nBest primer pairs for margin 1 (left):\n")
    log.write(write_primer_pairs(all_primers[1]))
    all_primers[2] = design_primers(margin2, genome)
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

    # check for Bam targets

    log.write(sep)
    log.close()
