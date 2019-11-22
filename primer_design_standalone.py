#!/usr/bin/env python3
"""primer_design_standalone.py
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser

from Bio import SeqIO
from Bio.Restriction import BamHI

from regions import Region
import vmatch
import primers


MARGIN_SIZE = 1000
INTERNAL_MARGIN = 200
sep = "-"*80 + "\n"


def check_BamHItargets_and_repeats(seq, repeats_list, L, direction):
    """Ensure that a region does not contain BamHI targets or repeats.
    
    Check whether a region on a genome contains BamHI targets or any
    sequence that may be substrate for homologous recombination (repeat
    equal or longer than L). If it does, shift it (right or left) until
    it does not.
    Args:
        seq (regions.Region): a region on the reference genome.
        repeats_list (set): list of possible HR substrate locations.
        L (int): length considered sufficient for HR events.
        direction (str, [left|right]): direction in which the region
            should be shifted.
    """
    seq_ok = False
    bam_ok = False
    repeat_ok = False
    while not seq_ok:
        # Check seq for BamHI target sites
        while not bam_ok:
            pos = BamHI.search(seq.subseq())
            if pos:
                offset = pos[0]
                seq.shift_past(offset, direction)
                repeat_ok = False
                globalpos = seq.s() + offset
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
                        offset = repeat[1] - seq.s()
                        seq.shift_past(offset, direction)
                    elif direction == "left":
                        offset = repeat[0] - seq.s()
                        seq.shift_past(offset, direction)
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
    # Parse command-line arguments and init constants
    parser = ArgumentParser(
        description="Design primers for a large genomic deletion by megapriming."
    )
    parser.add_argument(
        "-g", dest="GENOME", required=True,
        help="file containing genome sequence in FASTA format")
    parser.add_argument(
        "-l", dest="LOG_DIR", required=True,
        help="directory for output files")
    parser.add_argument(
        "-n", dest="DEL_NAME", required=True,
        help="name of desired deletion (e.g. D2)")
    parser.add_argument(
        "-d1", dest="DEL_START", required=True, type=int,
        help="start position of desired deletion")
    parser.add_argument(
        "-d2", dest="DEL_END", required=True, type=int,
        help="end position of desired deletion")
    parser.add_argument(
        "-L", dest="HR_LENGTH", metavar="HR_LENGTH", type=int, default=20,
        choices=range(15, 100),
        help=("min substrate length required for homologous recombination "
              "events (optional, default %(default)s bp)"))
    parser.add_argument(
        "-t", dest="LOCUS_TAG", default="",
        help=("locus tag of first deleted gene (optional, for primer naming "
              "purposes only)"))
    args = parser.parse_args()
    GENOME = args.GENOME
    LOG_DIR = args.LOG_DIR
    DEL_COORDS = (args.DEL_START, args.DEL_END)
    DEL_NAME = args.DEL_NAME
    HR_LENGTH = args.HR_LENGTH
    LOCUS_TAG = args.LOCUS_TAG

    # Check input
    try:
        genome = SeqIO.read(GENOME, "fasta")
    except (FileNotFoundError, ValueError):
        raise SystemExit("\n\terror: could not read genome file\n")
    try:
        log = open("%s/%s_primer_design.txt" % (LOG_DIR, DEL_NAME), "w")
    except FileNotFoundError:
        raise SystemExit("\n\terror: could not find output directory\n")
    if (not all([coord in range(1, len(genome)+1) for coord in DEL_COORDS])
            or DEL_COORDS[0] >= DEL_COORDS[1]):
        log.close()
        raise SystemExit("\n\terror: invalid deletion coordinates\n")

    # Define initial regions
    del_region = Region(DEL_COORDS, genome)
    margin1_coords = (del_region.s() - MARGIN_SIZE + INTERNAL_MARGIN,
                      del_region.s() + INTERNAL_MARGIN)
    margin1 = Region(margin1_coords, genome)
    margin2_coords = (del_region.e() - INTERNAL_MARGIN,
                      del_region.e() + MARGIN_SIZE - INTERNAL_MARGIN)
    margin2 = Region(margin2_coords, genome)
    log.write(
        sep + "##### PRIMER DESIGN FOR DELETION '%s' %d - %d (%.1f kb) #####\n"
        % (DEL_NAME, del_region.s(), del_region.e(),
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
        1: primers.design(margin1, crit_pos=(del_region.s(), "L")),
        2: primers.design(margin2, crit_pos=(del_region.e(), "R"))
    }
    log.write("\nBest primer pairs for margin 1 (left):\n")
    log.write(primers.write_pairs(all_primers[1]))
    log.write("\nBest primer pairs for margin 2 (right):\n")
    log.write(primers.write_pairs(all_primers[2]))

    # Choose primer pairs
    log.write("\nChoosing primers...")
    megapriming = primers.choose(all_primers, genome)
    log.write("\nDone.\n")
    log.write(sep)
    log.write("\nSelected pairs of primers:\n")
    for name, primer in megapriming.primers_raw_dict.items():
        log.write("%s: %s\n" % (name, primer.seq()))
    log.write("\nAdded tails:\n")
    for name, primer in megapriming.primers_tailed_dict.items():
        syst_name = primers.get_name(DEL_NAME, name, 1, LOCUS_TAG, 1)  # FIXME
        log.write("%s: %s\n" % (syst_name, primer))

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
    primers.save_pcr_regions(megapriming, DEL_NAME, LOG_DIR)
    log.write(
        "\nPCR region sequences saved at %s/%s_PCR_regions.fna.\n"
        % (LOG_DIR, DEL_NAME)
    )

    log.write(sep)
    log.close()
