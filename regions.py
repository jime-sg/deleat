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

    def __len__(self):
        """"""  # TODO
        return len(self.subseq())

    def global_seq(self):
        """"""  # TODO
        return self.globalseq

    def shift_past(self, position, direction):
        """"""  # TODO
        if direction == "right":
            offset = position
        elif direction == "left":
            offset = -(len(self.subsequence) - position)

        self.__init__((self.start + offset, self.end + offset),
                      self.globalseq)

    def overlap(self, coords):
        """Calculate length of overlap with a given location.

        Args:
            coords (tuple): start and end of location.
        Returns:
            i (int): length of overlap between region and given location.
        """
        i = min(self.e(), coords[1]) - max(self.s(), coords[0]) + 1
        if i < 0:
            i = 0
        return i
