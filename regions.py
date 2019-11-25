#!/usr/bin/env python3
"""regions.py
# TODO
@author: Jimena Solana
"""

from Bio.SeqFeature import SeqFeature, FeatureLocation


class Region:
    """Manage subsequence basic properties.

    Attributes:
        start (int): the region's start position on the global sequence.
        end (int): the region's end position on the global sequence.
        globalsequence (Bio.Seq.Seq): the sequence to which the region
            belongs.
        subsequence (Bio.Seq.Seq): the region's nucleotide sequence.
    """

    def __init__(self, coords, global_seq):
        """Init a region's basic information: start, end,
        globalsequence, subsequence.

        Args:
            coords (tuple(int, int)): start and end coordinates of the
                region.
            global_seq (Bio.Seq.Seq): the original sequence to which the
                region belongs.
        """
        self.start = coords[0]
        self.end = coords[1]
        self.globalsequence = global_seq
        # Biopython is 1-based but FeatureLocation takes Python
        # slicing-style positions: [20:30] -> 19..30
        feature = SeqFeature(FeatureLocation(self.start - 1, self.end))
        self.subsequence = feature.extract(global_seq).seq

    def subseq(self):
        """Return the region's nucleotide sequence (Bio.Seq.Seq)."""
        return self.subsequence

    def s(self):
        """Return the region's start position (int)."""
        return self.start

    def e(self):
        """Return the region's end position (int)."""
        return self.end

    def __len__(self):
        """Return the region's length (int)."""
        return len(self.subseq())

    def global_seq(self):
        """Return the original sequence to which the region belongs
        (Bio.Seq.Seq)."""
        return self.globalsequence

    def shift_past(self, position, direction):
        """Shift the region's location beyond a given position.
        
        The new location can be placed exactly to the left or to the
        right of the given position, so that the region does not overlap
        it.
        Args:
            position (int): boundary for the region's new location.
            direction (str): to which side of the position the region
                should be shifted ([left|right]).
        """
        if direction == "right":
            offset = position
        elif direction == "left":
            offset = -(len(self.subsequence) - position)

        self.__init__((self.start + offset, self.end + offset),
                      self.globalsequence)

    def overlap(self, coords):
        """Calculate length of overlap with a given location.

        Args:
            coords (tuple(int, int)): start and end of location.
        Returns:
            i (int): length of overlap between region and given location.
        """
        i = min(self.e(), coords[1]) - max(self.s(), coords[0]) + 1
        if i < 0:
            i = 0
        return i
