#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

import primer3
from regions import Region
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Primer:
    def __init__(self, p3results_dict, primer_n, primer_dir):
        self.id = primer_dir + "_" + str(primer_n)
        if primer_dir == "LEFT":
            self.start = p3results_dict["PRIMER_%s" % self.id][0]
            self.end = self.start + p3results_dict["PRIMER_%s" % self.id][1] - 1
        elif primer_dir == "RIGHT":
            self.end = p3results_dict["PRIMER_%s" % self.id][0]
            self.start = self.end - p3results_dict["PRIMER_%s" % self.id][1] + 1
        self.sequence = Seq(p3results_dict["PRIMER_%s_SEQUENCE" % self.id].upper())
        self.temp = p3results_dict["PRIMER_%s_TM" % self.id]
        self.penalty = p3results_dict["PRIMER_%s_PENALTY" % self.id]
        self.gc_perc = p3results_dict["PRIMER_%s_GC_PERCENT" % self.id]
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
            "PCR2_F": self.PCR2F,
            "PCR2_R": self.PCR2R
        }

    def add_tails(self):
        # PCR1_Bam-F: add BamHI target site at 5'
        self.PCR1Ft = Seq("gcacggatcc") + self.PCR1F.seq()
        # PCR1_R: unchanged
        self.PCR1Rt = self.PCR1R.seq()
        # PCR2_F: add revcomp of PCR1_R at 5'
        self.PCR2Ft = self.PCR1R.seq().reverse_complement().lower() + self.PCR2F.seq()
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


def design_primers(region, global_seq):
    primer_dict = {}
    p3_results = p3_design((region.s(), region.e()),
                           global_seq)
    n_pairs = p3_results["PRIMER_PAIR_NUM_RETURNED"]
    for i in range(n_pairs):
        for j in ("LEFT", "RIGHT"):
            primer = Primer(p3_results, i, j)
            primer_dict[primer.id] = primer

    return primer_dict


def p3_design(coords, template):
    start = coords[0]
    end = coords[1]
    p3_seqargs = {
        "SEQUENCE_TEMPLATE": str(template.seq),
        "SEQUENCE_INCLUDED_REGION": [start, end - start]
    }
    p3_globalargs = {
        "PRIMER_TASK": "generic",
        "PRIMER_PRODUCT_SIZE_RANGE": [750, 850],
        #"PRIMER_PRODUCT_OPT_SIZE": 800,
        #"PRIMER_PAIR_WT_PRODUCT_SIZE_LT": 0.2,
        #"PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 0.2,
        "PRIMER_NUM_RETURN": 20,
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


def save_pcr_regions(regions_dict, path):
    with open(path + "PCR_regions.fna", "w") as f:
        pcr1 = regions_dict["PCR1"]
        pcr1_r = SeqRecord(pcr1.subseq(), id="PCR1", description="%d:%d" % (pcr1.s(), pcr1.e()))
        pcr2 = regions_dict["PCR2"]
        pcr2_r = SeqRecord(pcr2.subseq(), id="PCR2", description="%d:%d" % (pcr2.s(), pcr2.e()))
        prod = pcr1.subseq() + pcr2.subseq()
        prod_r = SeqRecord(prod, id="total_product", description="")
        SeqIO.write((pcr1_r, pcr2_r, prod_r), f, "fasta")
