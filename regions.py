#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq


class Region:
    """
    Manage subsequence basic properties.
    """
    def __init__(self, coords, global_seq):
        """

        :rtype: Bio.SeqRecord
        """
        self.start = coords[0]
        self.end = coords[1]
        self.globalseq = global_seq
        self.seqfeature = SeqFeature(FeatureLocation(self.start, self.end))
        self.subsequence = self.seqfeature.extract(global_seq).seq

    def subseq(self):
        return self.subsequence
    def s(self):
        return self.start
    def e(self):
        return self.end

    def displace_past(self, pos, dir):
        if dir == "right":
            displ = pos
        elif dir == "left":
            displ = -(len(self.subsequence) - pos)

        self.__init__((self.start + displ, self.end + displ),
                      self.globalseq)


class Primer:
    def __init__(self, p3results_dict, primer_n, primer_dir):
        self.id = primer_dir + "_" + str(primer_n)
        self.start = p3results_dict["PRIMER_%s" % self.id][0]
        if primer_dir == "LEFT":
            self.end = self.start + p3results_dict["PRIMER_%s" % self.id][1] - 1
        elif primer_dir == "RIGHT":
            self.end = self.start - p3results_dict["PRIMER_%s" % self.id][1] + 1
        self.sequence = Seq(p3results_dict["PRIMER_%s_SEQUENCE" % self.id].upper())
        self.temp = p3results_dict["PRIMER_%s_TM" % self.id]
        self.penalty = p3results_dict["PRIMER_%s_PENALTY" % self.id]
        self.gc = p3results_dict["PRIMER_%s_GC_PERCENT" % self.id]
        self.self_any_th = p3results_dict["PRIMER_%s_SELF_ANY_TH" % self.id]
        self.self_end_th = p3results_dict["PRIMER_%s_SELF_END_TH" % self.id]
        self.hairpin_th = p3results_dict["PRIMER_%s_HAIRPIN_TH" % self.id]
        self.end_stability = p3results_dict["PRIMER_%s_END_STABILITY" % self.id]

    def seq(self):
        return self.sequence

    def s(self):
        return self.start

    def e(self):
        return self.end

    def tm(self):
        return self.temp


class PrimerSet:
    def __init__(self, primer_dict):
        self.PCR1F = primer_dict["1F"]
        self.PCR1R = primer_dict["1R"]
        self.PCR2F = primer_dict["2F"]
        self.PCR2R = primer_dict["2R"]

        self.primers_raw_dict = {
            "PCR1_F": self.PCR1F,
            "PCR1_R": self.PCR1R,
            "PCR2_F": self.PCR2R,
            "PCR2_R": self.PCR2R
        }

    def add_tails(self):
        # PCR1_Bam-F: add BamHI target site at 5'
        self.PCR1Ft = Seq("gcacggatcc") + self.PCR1F.seq()
        # PCR1_R: unchanged
        self.PCR1Rt = self.PCR1R.seq()
        # PCR2_F: add revcomp of PCR1_R at 5'
        self.PCR2Ft = self.PCR2F.seq().reverse_complement().lower() + self.PCR2F.seq()
        # PCR2_Bam-R: add BamHI target site at 5'
        self.PCR2Rt = Seq("gcacggatcc") + self.PCR2R.seq()

        self.primers_tailed_dict = {
            "PCR1_Ft": self.PCR1Ft,
            "PCR1_Rt": self.PCR1Rt,
            "PCR2_Ft": self.PCR2Ft,
            "PCR2_Rt": self.PCR2Rt
        }

    def get_PCR_regions(self, global_seq):
        pcr1_start = self.PCR1F.s()
        pcr1_end = self.PCR1R.e()
        pcr2_start = self.PCR2F.s()
        pcr2_end = self.PCR2R.e()
        self.PCR1_region = Region((pcr1_start, pcr1_end), global_seq)
        self.PCR2_region = Region((pcr2_start, pcr2_end), global_seq)
        pcr_dict = {"PCR1": self.PCR1_region,
                    "PCR2": self.PCR2_region}

        return pcr_dict

