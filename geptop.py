#!/usr/bin/env python3
"""
@author: Jimena Solana
"""

from collections import Counter
from math import sqrt
from time import time  # FIXME

from Bio import SeqIO
# from Bio.Blast.Applications import NcbiblastxCommandline
from more_itertools import windowed


CV_K = 6

def run(cds_file, cutoff, n_proc, out_path):
    pass
    # seqs = SeqIO.parse(open(cds_file), "fasta")


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
    freqs_k = Counter()  # K-words
    freqs_km1 = Counter()  # (K-1)-words
    freqs_km2 = Counter()  # (K-2)-words
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


def distance(cv1, cv2):
    a = sum(
        {kword: cv1[kword] * cv2[kword]
         for kword in cv1.keys() & cv2.keys()}.values()
    )
    b = sum(map(lambda x: x*x, cv1.values()))
    c = sum(map(lambda x: x*x, cv2.values()))

    corr = a / sqrt(b * c)
    dist = (1 - corr) / 2
    return dist


"""
def StrToNum(string):
    aa_dict = {'A':0,  'C':1,  'D':2,  'E':3,  'F':4,  'G':5,
               'H':6,  'I':7,  'K':8,  'L':9,  'M':10, 'N':11,
               'P':12, 'Q':13, 'R':14, 'S':15, 'T':16, 'V':17,
               'W':18, 'Y':19,
               'B':2,  'U':1,  'X':5,  'Z':3,  'J':7}
    number = 0
    for i in range(len(string)):
        if string[i] not in aa_dict:
            return -1
        bit = aa_dict[string[i]]
        power = len(string)-i-1
        number += bit*(20**power)
    return number


def CompositionVector(species):
    kstring = 6
    k = {}
    k0 = {}
    k1 = {}
    k2 = {}
    try:
        for seqrecord in SeqIO.parse(species, "fasta"):
            seq = str(seqrecord.seq)
            len0 = len(seq) - kstring + 3
            for s in range(len0):
                start = s
                end = kstring + s - 2
                num = StrToNum(seq[start:end])
                if num not in k2:
                    k2[num] = 1
                else:
                    k2[num] += 1
                if s < len0 - 2:
                    num = StrToNum(seq[start:end+2])
                    if num not in k0:
                        k0[num] = 1
                    else:
                        k0[num] += 1
                if s < len0 - 1:
                    num = StrToNum(seq[start:end+1])
                    if num not in k1:
                        k1[num] = 1
                    else:
                        k1[num] += 1
            if -1 in k0:
                del k0[-1]
            if -1 in k1:
                del k1[-1]
            if -1 in k2:
                del k2[-1]

            string0 = sum(k0.values())
            string1 = sum(k1.values())
            string2 = sum(k2.values())

        for n1 in k1:
            for aa in range(20):
                n0 = n1 * 20 + aa
                n2 = n1 % (20**(kstring-2)) * 20 + aa
                if n2 in k1:
                    if n0 in k0:
                        n3 = n1 % (20**(kstring-2))
                        p0 = 1.0 * k1[n1] * k1[n2] * string0 * string2/k2[n3]/string1/string1
                        k[n0] = (k0[n0]-p0)/p0
                    else:
                        k[n0] =- 1
        return k

    except (Exception, IOError):
        print ('Compute CV failed to: %s'%IOError)


def Distance(CV1, CV2):
    O = 0
    P = 0
    Q = 0

    for value in CV1:
        if value in CV2:
            O += CV1[value]*CV2[value]
            P += CV1[value]*CV1[value]
            Q += CV2[value]*CV2[value]
        else:
            P += CV1[value]*CV1[value]

    for value in CV2:
        if value not in CV1:
            Q += CV2[value]*CV2[value]

    return (1-O/sqrt(P*Q))/2
"""


if __name__ == "__main__":
    t0 = time()
    comp_vector1 = composition_vector("/home/jimena/Escritorio/bq1.faa", 6)
    comp_vector2 = composition_vector("/home/jimena/Escritorio/bq2.faa", 6)
    t1 = time()
    d = distance(comp_vector1, comp_vector2)
    t2 = time()
    print(d, "%.2f s" % (t1-t0), "%.2f s" % (t2-t1))

    # t0 = time()
    # comp_vector1 = CompositionVector("/home/jimena/Escritorio/bq1.faa")
    # comp_vector2 = CompositionVector("/home/jimena/Escritorio/bq2.faa")
    # t1 = time()
    # d = Distance(comp_vector1, comp_vector2)
    # t2 = time()
    # print(d, "%.2f s" % (t1-t0), "%.2f s" % (t2-t1))
    # for x, y in comp_vector.items():
        # print(x, y)
    # print(len(comp_vector))
