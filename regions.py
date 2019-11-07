#!/usr/bin/env python3
"""regions.py
# TODO
@author: Jimena Solana
"""

from Bio.SeqFeature import SeqFeature, FeatureLocation


class Region:
    """
    Manage subsequence basic properties. # TODO
    """

    def __init__(self, coords, global_seq):
        """"""  # TODO
        self.start = coords[0]
        self.end = coords[1]
        self.globalseq = global_seq
        # Biopython is 1-based but FeatureLocation takes Python slicing-style
        # positions: [20:30] -> 19..30
        self.seqfeature = SeqFeature(FeatureLocation(self.start - 1, self.end))
        self.subsequence = self.seqfeature.extract(global_seq).seq

    def subseq(self):
        """"""  # TODO
        return self.subsequence

    def s(self):
        """"""  # TODO
        return self.start

    def e(self):
        """"""  # TODO
        return self.end

    def global_seq(self):
        """"""  # TODO
        return self.globalseq

    def shift_past(self, position, direction):
        """"""  # TODO
        if direction == "right":
            displ = position
        elif direction == "left":
            displ = -(len(self.subsequence) - position)

        self.__init__((self.start + displ, self.end + displ),
                      self.globalseq)
