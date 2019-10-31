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

    def __init__(self, coords, global_seq):
        """

        :rtype: Bio.SeqRecord
        """
        self.start = coords[0]
        self.end = coords[1]
        self.globalseq = global_seq
        self.seqfeature = SeqFeature(FeatureLocation(self.start - 1, self.end))
        self.subsequence = self.seqfeature.extract(global_seq).seq

    def subseq(self):
        return self.subsequence

    def s(self):
        return self.start

    def e(self):
        return self.end

    def global_seq(self):
        return self.globalseq

    def displace_past(self, position, direction):
        if direction == "right":
            displ = position
        elif direction == "left":
            displ = -(len(self.subsequence) - position)

        self.__init__((self.start + displ, self.end + displ),
                      self.globalseq)
