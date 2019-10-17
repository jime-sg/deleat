#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

from Bio.SeqFeature import SeqFeature, FeatureLocation


class Region:
    """
    Manage subsequences basic properties.
    """
    def __init__(self, seq, coords):
        """

        :rtype: Bio.SeqRecord
        """
        self.origseq = seq
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

    def displace_past(self, pos, dir):
        if dir == "right":
            displ = pos
        if dir == "left":
            displ = -(len(self.subseq) - pos)

        self.__init__(self.origseq, (self.s + displ, self.e + displ))