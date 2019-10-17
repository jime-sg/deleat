#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

from Bio.SeqFeature import SeqFeature, FeatureLocation


class Region(seq, coords):
    def __init__(self):
        self.seqfeature = SeqFeature(FeatureLocation(start, end))
        self.subseq = self.seqfeature.extract(seq)
        self.s = coords[0]
        self.e = coords[1]

    def subseq(self):
        return self.subseq
    def s(self):
        return self.s
    def e(self):
        return self.e
