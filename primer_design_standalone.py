#!/usr/bin/env python3
"""primer_design_standalone.py
# TODO
@author: Jimena Solana
"""

from Bio import SeqIO
from Bio.Restriction import BamHI

from regions import Region
import vmatch
from primers import (design_primers, write_primer_pairs, choose_primers,
                     save_pcr_regions)


GENOME = "/home/jimena/Bartonella/NC_005955.fna"
LOG_DIR = "/home/jimena/Bartonella/deletions/deletion2/"
DEL_COORDS = (33900, 80000)
MARGIN_SIZE = 1000
INTERNAL_MARGIN = 200
HR_LENGTH = 20

sep = "-"*80 + "\n"


def check_BamHItargets_and_repeats(seq, repeats_list, L, direction):
    """"""  # TODO
    seq_ok = False
    bam_ok = False
    repeat_ok = False
    while not seq_ok:
        # Check seq for BamHI target sites
        while not bam_ok:
            pos = BamHI.search(seq.subseq())
            if pos:
                pos = pos[0]
                globalpos = pos + seq.s()
                seq.shift_past(pos, direction)
                repeat_ok = False
                log.write(
                    "!! BamHI target site found at %d. "
                    "Location reset to %d - %d.\n"
                    % (globalpos, seq.s(), seq.e())
                )
            else:
                bam_ok = True

        # Check seq for repetitive regions
        while not repeat_ok:
            repeat_ok = True
            for repeat in repeats_list:
                if seq.overlap(repeat) > L:
                    if direction == "right":
                        seq.shift_past(repeat[1], direction)
                    elif direction == "left":
                        seq.shift_past(repeat[0], direction)
                    repeat_ok = False
                    bam_ok = False
                    log.write(
                        "!! Repeated sequence found at (%d-%d). "
                        "Location reset to %d - %d.\n"
                        % (repeat[0], repeat[1], seq.s(), seq.e())
                    )
                    break

        if bam_ok and repeat_ok:
            seq_ok = True


if __name__ == "__main__":
    genome = SeqIO.read(GENOME, "fasta")
    log = open(LOG_DIR + "primer_design.txt", "w")

    # Define initial regions
    del_region = Region(DEL_COORDS, genome)
    margin1_coords = (del_region.s() - MARGIN_SIZE + INTERNAL_MARGIN,
                      del_region.s() + INTERNAL_MARGIN)
    margin1 = Region(margin1_coords, genome)
    margin2_coords = (del_region.e() - INTERNAL_MARGIN,
                      del_region.e() + MARGIN_SIZE - INTERNAL_MARGIN)
    margin2 = Region(margin2_coords, genome)
    log.write(
        sep + "##### PRIMER DESIGN FOR DELETION %d - %d (%.1f kb) #####\n"
        % (del_region.s(), del_region.e(),
           (del_region.e()-del_region.s())/1000)
    )

    # Calculate list of all repeats of length >= HR_LENGTH
    # (possible substrates for homologous recombination)
    repeats_list = vmatch.run(GENOME, HR_LENGTH, LOG_DIR)

    # Check margin region quality
    log.write(sep + "Checking margins...\n")
    log.write(
        "Initial locations:\n\tMargin 1 (left):  %d - %d\n\t"
        "Margin 2 (right): %d - %d\n"
        % (margin1.s(), margin1.e(), margin2.s(), margin2.e())
    )
    log.write("\nChecking margin 1...\n")
    check_BamHItargets_and_repeats(margin1,
                                   repeats_list, L=HR_LENGTH,
                                   direction="right")
    log.write("\nChecking margin 2...\n")
    check_BamHItargets_and_repeats(margin2,
                                   repeats_list, L=HR_LENGTH,
                                   direction="left")
    log.write("\nDone.")
    log.write(
        "\nFinal locations:\n\tMargin 1 (left):  %d - %d\n\t"
        "Margin 2 (right): %d - %d\n"
        % (margin1.s(), margin1.e(), margin2.s(), margin2.e())
    )

    # Design primers
    log.write(sep + "Designing primers...\n")
    all_primers = {
        1: design_primers(margin1, crit_pos=(del_region.s(), "L")),
        2: design_primers(margin2, crit_pos=(del_region.e(), "R"))
    }
    log.write("\nBest primer pairs for margin 1 (left):\n")
    log.write(write_primer_pairs(all_primers[1]))
    log.write("\nBest primer pairs for margin 2 (right):\n")
    log.write(write_primer_pairs(all_primers[2]))

    # Choose primer pairs
    log.write("\nChoosing primers...")
    megapriming = choose_primers(all_primers, genome)
    log.write("\nDone.\n")
    log.write(sep)
    log.write("\nSelected pairs of primers:\n")
    for name, primer in megapriming.primers_raw_dict.items():
        log.write("%s: %s\n" % (name, primer.seq()))
    log.write("\nAdded tails:\n")
    for name, primer in megapriming.primers_tailed_dict.items():
        log.write("%s: %s\n" % (name, primer))

    # Define PCR regions
    pcr1 = megapriming.PCR_dict["PCR1"]
    pcr2 = megapriming.PCR_dict["PCR2"]
    log.write("\nDetermined PCR regions for megapriming:\n")
    log.write(
        "\tPCR1 (left):  %d - %d (%d bp)\n"
        "\tPCR2 (right): %d - %d (%d bp)\n"
        % (pcr1.s(), pcr1.e(), pcr1.e()-pcr1.s(),
           pcr2.s(), pcr2.e(), pcr2.e()-pcr2.s())
    )
    save_pcr_regions(megapriming, LOG_DIR)
    log.write("\nPCR region sequences saved at %sPCR_regions.fna.\n" % LOG_DIR)

    log.write(sep)
    log.close()
