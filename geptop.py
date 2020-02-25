#!/usr/bin/env python3
"""
@author: Jimena Solana
"""

from collections import Counter
from itertools import product

from Bio import SeqIO
from more_itertools import windowed


CV_K = 6

def run(cds_file, cutoff, n_proc, out_path):
    seqs = SeqIO.parse(open(cds_file), "fasta")


def word2num(word):
    num = 0
    for i, n in enumerate(word):
        num += n * (26**i)
    return num


def num2word(num):
    word = []
    while num:
        word.append(num % 26)
        num = num // 26
    return tuple(word)


def composition_vector(species_fasta, K):
    freqs_k = Counter()
    freqs_km1 = Counter()
    freqs_km2 = Counter()
    total_words = 0

    for protein in SeqIO.parse(species_fasta, "fasta"):
        total_words += len(protein) - K + 1
        protein_n = tuple(
            ord(c)-64 for c in protein.seq.upper()
            if ord(c)-64 in range(1, 27)
        )
        for kword in windowed(protein_n, K):
            kword_n = word2num(kword)
            km1word_n = word2num(kword[:-1])
            km2word_n = word2num(kword[:-2])
            freqs_k[kword_n] += 1
            freqs_km1[km1word_n] += 1
            freqs_km2[km2word_n] += 1
        # last words of length K-1 and K-2
        freqs_km1[word2num(protein_n[-K+1:])] += 1
        freqs_km2[word2num(protein_n[-K+1:-1])] += 1
        freqs_km2[word2num(protein_n[-K+2:])] += 1

    probs_k = {w: f/total_words for w, f in freqs_k.items()}
    probs_km1 = {w: f/total_words for w, f in freqs_km1.items()}
    probs_km2 = {w: f/total_words for w, f in freqs_km2.items()}

    probs_0 = dict()
    for kword_n in probs_k:
        # kword_n="ABCDE" -> a="ABCD"; b="BCDE"; c="BCD"
        a = word2num(num2word(kword_n)[:-1])
        b = word2num(num2word(kword_n)[1:])
        c = word2num(num2word(kword_n)[1:-1])
        probs_0[kword_n] = probs_km1[a] * probs_km1[b] / probs_km2[c]

    cv = dict()
    for kword_n in probs_k:
        if probs_0[kword_n] == 0:
            cv[kword_n] = 0
        else:
            cv[kword_n] = ((probs_k[kword_n] - probs_0[kword_n])
                           / probs_0[kword_n])
    return cv


if __name__ == "__main__":
    pass
    #comp_vector = composition_vector("/home/jimena/Bartonella/CDSa/all.faa", 6)
    #for x, y in comp_vector.items():
        #print(x, y)
    #print(len(comp_vector))
