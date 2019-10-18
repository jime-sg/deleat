#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Jimena Solana
"""

from Bio.SeqFeature import SeqFeature, FeatureLocation


class Region:
    """
    Manage subsequence basic properties.
    """
    def __init__(self, seq, coords):
        """

        :rtype: Bio.SeqRecord
        """
        self.origseq = seq
        self.start = coords[0]
        self.end = coords[1]
        self.seqfeature = SeqFeature(FeatureLocation(self.start, self.end))
        self.subsequence = self.seqfeature.extract(seq).seq

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

        self.__init__(self.origseq,
                      (self.start + displ, self.end + displ))