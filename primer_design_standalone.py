#!/usr/bin/env python3
"""primer_design_standalone.py
# TODO
@author: Jimena Solana
"""

from Bio import SeqIO
from Bio.Restriction import BamHI

from regions import Region
from primers import (design_primers, write_primer_pairs, choose_primers,
                     save_pcr_regions)


def check_BamHItargets_and_repeats(seq, direction):
    """"""  # TODO
    seq_ok = False
    while not seq_ok:
        bam_ok = False
        repeat_ok = True  # FIXME

        # Check seq for BamHI target sites
        while not bam_ok:
            pos = BamHI.search(seq.subseq())
            if pos:
                pos = pos[0]
                globalpos = pos + seq.s()
                seq.displace_past(pos, direction)
                log.write(
                    "!! BamHI target site found at %d. "
                    "Location reset to %d - %d.\n"
                    % (globalpos, seq.s(), seq.e())
                )
            else:
                bam_ok = True

        # Check seq for repetitive regions
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

    # Define initial regions
    del_region = Region(DEL_COORDS, genome)
    margin1_coords = (del_region.s() - MARGIN_SIZE,
                      del_region.s())
    margin1 = Region(margin1_coords, genome)
    margin2_coords = (del_region.e(),
                      del_region.e() + MARGIN_SIZE)
    margin2 = Region(margin2_coords, genome)

    # Check margin region quality
    log.write(sep + "Checking margins...\n")
    log.write(
        "Initial locations:\n\tMargin 1 (left):  %d - %d\n\t"
        "Margin 2 (right): %d - %d\n"
        % (margin1.s(), margin1.e(), margin2.s(), margin2.e())
    )
    log.write("\nChecking margin 1...\n")
    check_BamHItargets_and_repeats(margin1, direction="right")
    log.write("\nChecking margin 2...\n")
    check_BamHItargets_and_repeats(margin2, direction="left")
    log.write("\nMargin quality check: done.")
    log.write(
        "\nFinal locations:\n\tMargin 1 (left):  %d - %d\n\t"
        "Margin 2 (right): %d - %d\n"
        % (margin1.s(), margin1.e(), margin2.s(), margin2.e())
    )

    # Design primers
    log.write(sep + "Designing primers...\n")
    all_primers = {
        1: design_primers(margin1),
        2: design_primers(margin2)
    }
    log.write("\nBest primer pairs for margin 1 (left):\n")
    log.write(write_primer_pairs(all_primers[1]))
    log.write("\nBest primer pairs for margin 2 (right):\n")
    log.write(write_primer_pairs(all_primers[2]))

    # Choose primer pairs
    log.write("\nChecking for BamHI targets in the megapriming product...\n")
    megapriming = choose_primers(all_primers, genome)
    log.write(sep)
    log.write("Selected pairs of primers:\n")
    for name, primer in megapriming.primers_raw_dict.items():
        log.write("%s: %s\n" % (name, primer.seq()))
    log.write("\nAdded tails:\n")
    for name, primer in megapriming.primers_tailed_dict.items():
        log.write("%s: %s\n" % (name, primer))

    # Define PCR regions
    pcr1 = megapriming.PCR_dict["PCR1"]
    pcr2 = megapriming.PCR_dict["PCR2"]
    log.write("\nDetermined PCR regions for megapriming:\n")
    log.write("\tPCR1 (left):  %d - %d\n\tPCR2 (right): %d - %d\n"
              % (pcr1.s(), pcr2.e(), pcr2.s(), pcr2.e()))
    save_pcr_regions(megapriming, LOG_DIR)
    log.write("\nPCR region sequences saved at %sPCR_regions.fna.\n" % LOG_DIR)

    log.write(sep)
    log.close()
