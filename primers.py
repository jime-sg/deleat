#!/usr/bin/env python3
"""primers.py
# TODO
@author: Jimena Solana
"""

from itertools import product

import primer3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import BamHI

from regions import Region


# Primer3 settings
PCR_REGION_SIZE_RANGE = [600, 1000]
PCR_REGION_OPT_SIZE = 800
LT_SIZE_WT = 0.005
GT_SIZE_WT = 0.005
N_PRIMERS = 20

PCR_REGIONS_MAX_SIZEDIFF = 150


class Primer:
    """
    # TODO
    Attributes:
        id (str, [LEFT|RIGHT]_[0-9]+): the primer's unique name.
        start (int): position of the primer's leftmost edge on the
            template sequence.
        end (int): position of the primer's rightmost edge on the
            template sequence.
        sequence (Bio.Seq.Seq): the primer's sequence.
        temp (float): the primer's melting temperature.
        penalty (float): penalty calculated by primer3.
        gc_perc (float): GC% content.
        self_any_th:  # TODO
        self_end_th:  # TODO
        hairpin_th:  # TODO
        end_stability:  # TODO
    """

    def __init__(self, p3results_dict, primer_n, primer_dir):
        """Collect all information returned by primer3 about a primer.
        
        Parse results dictionary returned by primer3 in order to extract
        all information about a primer, and save it to instance
        attributes: inits Primer with id, start, end, sequence, temp,
        penalty, gc_perc, self_any_th, self_end_th, hairpin_th and
        end_stability.
        Args:
            p3results_dict (dict): results returned by
                primer3.bindings.designPrimers(), with all information
                about designed primer pairs.
            primer_n (int): index of primer pair to which the primer of
                interest belongs.
            primer_dir (str, [LEFT|RIGHT]): direction of the primer of
                interest.
        """
        self.id = primer_dir + "_" + str(primer_n)
        if primer_dir == "LEFT":  # Position returned by primer3 = start
            self.start = p3results_dict["PRIMER_%s" % self.id][0]
            self.end = self.start + p3results_dict["PRIMER_%s" % self.id][1]-1
        elif primer_dir == "RIGHT":  # Position returned by primer3 = end
            self.end = p3results_dict["PRIMER_%s" % self.id][0]
            self.start = self.end - p3results_dict["PRIMER_%s" % self.id][1]+1
        self.sequence = Seq(p3results_dict["PRIMER_%s_SEQUENCE" % self.id].upper())
        self.temp = p3results_dict["PRIMER_%s_TM" % self.id]
        self.penalty = p3results_dict["PRIMER_%s_PENALTY" % self.id]
        self.gc_perc = p3results_dict["PRIMER_%s_GC_PERCENT" % self.id]
        self.self_any_th = p3results_dict["PRIMER_%s_SELF_ANY_TH" % self.id]
        self.self_end_th = p3results_dict["PRIMER_%s_SELF_END_TH" % self.id]
        self.hairpin_th = p3results_dict["PRIMER_%s_HAIRPIN_TH" % self.id]
        self.end_stability = p3results_dict["PRIMER_%s_END_STABILITY" % self.id]

    def seq(self):
        """Return the primer's sequence."""
        return self.sequence

    def s(self):
        """Return the primer's leftmost position on the template
        sequence.
        """
        return self.start

    def e(self):
        """Return the primer's rightmost position on the template
        sequence.
        """
        return self.end

    def tm(self):
        """Return the primer's melting temperature."""
        return self.temp


class PrimerSet:
    """A set of four primers that defines two PCR regions for megapriming.
    
    Attributes:
        PCR1F, PCR1R, PCR2F, PCR2R (regions.Region): forward and reverse
            primers for PCR1 (left), and forward and reverse primers for PCR2
            (right), respectively.
        primers_raw_dict (dict): dictionary of untailed (raw) primers.
        PCR1Ft, PCR1Rt, PCR2Ft, PCR2Rt (Bio.Seq.Seq): primers with added tails
            necessary for megapriming.
        primers_tailed_dict (dict): dictionary of tailed primers.
        PCR_dict (dict): dictionary of PCR regions defined by the primer set.
        PCR1_region (regions.Region): PCR1 (left) region defined by primers
            PCR1F and PCR1R.
        PCR2_region (regions.Region): PCR2 (left) region defined by primers
            PCR2F and PCR2R.
    """

    def __init__(self, primer_dict, global_seq):
        """Define the set of four primers, add tails and delimit PCR regions.

        Inits PrimerSet with PCR1F, PCR1R, PCR2F, PCR2R, primers_raw_dict,
        primers_tailed_dict, PCR1Ft, PCR1Rt, PCR2Ft, PCR2Rt, PCR_dict,
        PCR1_region and PCR2_region.
        Args:
             primer_dict (dict): dictionary of 1F, 1R, 2F and 2R primers.
             global_seq (Bio.Seq.Seq): template sequence.
        """
        self.PCR1F = primer_dict["1F"]
        self.PCR1R = primer_dict["1R"]
        self.PCR2F = primer_dict["2F"]
        self.PCR2R = primer_dict["2R"]
        self.primers_raw_dict = {
            "PCR1_F": self.PCR1F,
            "PCR1_R": self.PCR1R,
            "PCR2_F": self.PCR2F,
            "PCR2_R": self.PCR2R
        }
        self.primers_tailed_dict = self.add_tails()
        self.PCR1Ft = self.primers_tailed_dict["PCR1_Ft"]
        self.PCR1Rt = self.primers_tailed_dict["PCR1_Rt"]
        self.PCR2Ft = self.primers_tailed_dict["PCR2_Ft"]
        self.PCR2Rt = self.primers_tailed_dict["PCR2_Rt"]
        self.PCR_dict = self.get_PCR_regions(global_seq)
        self.PCR1_region = self.PCR_dict["PCR1"]
        self.PCR2_region = self.PCR_dict["PCR2"]

    def add_tails(self):
        """Add tails necessary for megapriming to corresponding primers.

        - PCR1_Bam-F: add BamHI target site at 5'
        - PCR1_R: unchanged
        - PCR2_F: add reverse complement of PCR1_R at 5'
        - PCR2_Bam-R: add BamHI target site at 5'
        Returns:
             primers_tailed_dict (dict): dictionary of tailed primers.
        """
        primers_tailed_dict = {
            "PCR1_Ft": Seq("gcacggatcc") + self.PCR1F.seq(),
            "PCR1_Rt": self.PCR1R.seq(),
            "PCR2_Ft": self.PCR1R.seq().reverse_complement().lower()
                       + self.PCR2F.seq(),
            "PCR2_Rt": Seq("gcacggatcc") + self.PCR2R.seq()
        }
        return primers_tailed_dict

    def get_PCR_regions(self, global_seq):
        """Define PCR regions (left+right) delimited by the primer set.
        
        Args:
            global_seq: 

        Returns:

        """  # TODO
        pcr1_start = self.PCR1F.s()
        pcr1_end = self.PCR1R.e()
        pcr2_start = self.PCR2F.s()
        pcr2_end = self.PCR2R.e()
        pcr1_region = Region((pcr1_start, pcr1_end), global_seq)
        pcr2_region = Region((pcr2_start, pcr2_end), global_seq)
        pcr_dict = {"PCR1": pcr1_region,
                    "PCR2": pcr2_region}
        return pcr_dict

    def get_product(self):
        """"""  # TODO
        return self.PCR1_region.subseq() + self.PCR2_region.subseq()


def design_primers(region, crit_pos):
    """"""  # TODO
    primer_dict = {}
    p3_results = p3_design(region, crit_pos)
    n_pairs = p3_results["PRIMER_PAIR_NUM_RETURNED"]
    for i in range(n_pairs):
        for j in ("LEFT", "RIGHT"):
            primer = Primer(p3_results, i, j)
            primer_dict[primer.id] = primer
    return primer_dict


def p3_design(region, crit_pos):
    """"""  # TODO
    start = region.s()
    end = region.e()
    x = crit_pos[0]
    d = crit_pos[1]
    if d == "L":  # PCR1 (left)
        # Reverse primer must be between critical position and margin end
        y = end-x
        region_list = [-1, -1, x, y]
    elif d == "R":  # PCR2 (right)
        # Forward primer must be between margin start and critical position
        y = x-start
        region_list = [x-y, y, -1, -1]

    p3_seqargs = {
        "SEQUENCE_TEMPLATE": str(region.global_seq().seq),
        "SEQUENCE_INCLUDED_REGION": [start, end-start],
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": region_list
    }
    p3_globalargs = {
        "PRIMER_TASK": "generic",
        "PRIMER_PRODUCT_SIZE_RANGE": PCR_REGION_SIZE_RANGE,
        "PRIMER_PRODUCT_OPT_SIZE": PCR_REGION_OPT_SIZE,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": LT_SIZE_WT,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": GT_SIZE_WT,
        "PRIMER_NUM_RETURN": N_PRIMERS,
        "PRIMER_FIRST_BASE_INDEX": 1  # 1-based indexing
    }
    results = primer3.bindings.designPrimers(p3_seqargs, p3_globalargs)
    return results


def write_primer_pairs(primer_dict):
    """Convert dictionary of primer design results into a readable table.

    Args:
        primer_dict (dict): primer design results for one region.
    Returns:
        results_str (str): same results, parsed and readable.
    """
    results_str = ("N\tsize\tF_start-F_end..R_start-R_end\tseq_F\tseq_R\t"
                   "TM_F (ºC)\tTM_R (ºC)\n")
    n_pairs = len(primer_dict)//2
    for n in range(n_pairs):
        lp = primer_dict["LEFT_%d" % n]
        rp = primer_dict["RIGHT_%d" % n]
        prod_size = rp.e() - lp.s() + 1
        with open("/home/jimena/Escritorio/product_lengths.txt", "a") as f:  # FIXME
            f.write(str(prod_size) + "\n")

        results_str += (
                "%d\t%d\t%d-%d..%d-%d\t%-24s\t%-24s\t%.1f\t%.1f\n"
                % (n+1, prod_size, lp.s(), lp.e(), rp.s(), rp.e(),
                   lp.seq(), rp.seq(), lp.tm(), rp.tm())
        )
    return results_str


def choose_primers(primer_dict, global_seq):
    """"""  # TODO
    n_pairs = len(primer_dict[1])//2
    combinations = list(product(range(n_pairs), repeat=2))
    # Prioritize primer pair quality = minimize sum of primer pair indexes
    combinations.sort(key=lambda x: x[0]+x[1])
    i = 0
    sizediff_ok = False
    bam_ok = False
    while not sizediff_ok and not bam_ok:
        chosen_primers = {
            "1F": primer_dict[1]["LEFT_%d" % combinations[i][0]],
            "1R": primer_dict[1]["RIGHT_%d" % combinations[i][0]],
            "2F": primer_dict[2]["LEFT_%d" % combinations[i][1]],
            "2R": primer_dict[2]["RIGHT_%d" % combinations[i][1]]
        }
        primer_set = PrimerSet(chosen_primers, global_seq)
        mp_product = primer_set.get_product()
        size_diff = abs(len(primer_set.PCR1_region)
                        - len(primer_set.PCR2_region))

        if size_diff > PCR_REGIONS_MAX_SIZEDIFF:
            i += 1
            bam_ok = False
            continue
        else:
            sizediff_ok = True
        if BamHI.search(mp_product):
            i += 1
            sizediff_ok = False
        else:
            bam_ok = True

        with open("/home/jimena/Escritorio/product_diffs.txt", "a") as f:  # FIXME
            f.write(str(size_diff) + "\n")

    return primer_set


def save_pcr_regions(primer_set, path):
    """Save defined PCR regions (and total product) to a FASTA file.

    Args:
        primer_set (primers.PrimerSet): chosen primers.
        path (str): output file path.
    """
    with open("%s/PCR_regions.fna" % path, "w") as f:
        pcr1 = primer_set.PCR1_region
        pcr2 = primer_set.PCR2_region
        prod = primer_set.get_product()
        pcr1_r = (
            SeqRecord(pcr1.subseq(), id="PCR1", description="%d:%d"
            % (pcr1.s(), pcr1.e()))
        )
        pcr2_r = (
            SeqRecord(pcr2.subseq(), id="PCR2", description="%d:%d"
            % (pcr2.s(), pcr2.e()))
        )
        prod_r = SeqRecord(prod, id="total_product", description="")
        SeqIO.write((pcr1_r, pcr2_r, prod_r), f, "fasta")
