#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from collections import Counter
from math import sqrt
from time import time  # FIXME
import json
import os
from glob import glob

from Bio import SeqIO
from Bio.Blast.Applications import (NcbimakeblastdbCommandline,
                                    NcbiblastpCommandline)
from Bio.Blast import NCBIXML
from more_itertools import windowed


CV_K = 6
BLAST_EVALUE = 1
ORG_NAMES = {
    "DEG1001": "Bacillus subtilis 168",
    "DEG1002": "Staphylococcus aureus N315",
    "DEG1003": "Vibrio cholerae N16961",
    "DEG1005": "Haemophilus influenzae Rd KW20",
    "DEG1006": "Mycoplasma genitalium G37",
    "DEG1007": "Streptococcus pneumoniae",
    "DEG1008": "Helicobacter pylori 26695",
    "DEG1010": "Mycobacterium tuberculosis H37Rv",
    "DEG1011": "Salmonella typhimurium LT2",
    "DEG1012": "Francisella novicida U112",
    "DEG1013": "Acinetobacter baylyi ADP1",
    "DEG1014": "Mycoplasma pulmonis UAB CTIP",
    "DEG1015": "Pseudomonas aeruginosa UCBPP-PA14",
    "DEG1016": "Salmonella enterica serovar Typhi",
    "DEG1017": "Staphylococcus aureus NCTC 8325",
    "DEG1018": "Escherichia coli MG1655 I",
    "DEG1019": "Escherichia coli MG1655 II",
    "DEG1020": "Caulobacter crescentus",
    "DEG1021": "Streptococcus sanguinis",
    "DEG1022": "Porphyromonas gingivalis ATCC 33277",
    "DEG1023": "Bacteroides thetaiotaomicron VPI-5482",
    "DEG1024": "Burkholderia thailandensis E264",
    "DEG1025": "Mycobacterium tuberculosis H37Rv II",
    "DEG1026": "Salmonella enterica subsp. enterica serovar Typhimurium str. 14028S",
    "DEG1027": "Mycobacterium tuberculosis H37Rv III",
    "DEG1028": "Sphingomonas wittichii RW1",
    "DEG1029": "Shewanella oneidensis MR-1",
    "DEG1030": "Pseudomonas aeruginosa PAO1",
    "DEG1031": "Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819",
    "DEG1032": "Salmonella enterica serovar Typhimurium SL1344",
    "DEG1033": "Salmonella enterica serovar Typhi Ty2",
    "DEG1034": "Bacteroides fragilis 638R",
    "DEG1035": "Burkholderia pseudomallei K96243",
    "DEG1036": "Pseudomonas aeruginosa PAO1",
    "DEG1037": "Streptococcus pyogenes MGAS5448",
    "DEG1038": "Streptococcus pyogenes NZ131",
    "DEG1039": "Porphyromonas gingivalis ATCC 33277",
    "DEG1040": "Synechococcus elongatus PCC 7942",
    "DEG1041": "Rhodopseudomonas palustris CGA009",
    "DEG1042": "Streptococcus agalactiae A909",
    "DEG1043": "Acinetobacter baumannii ATCC 17978",
    "DEG1044": "Acinetobacter baumannii ATCC 17978",
    "DEG1045": "Agrobacterium fabrum str. C58",
    "DEG1046": "Brevundimonas subvibrioides ATCC 15264",
    "DEG1047": "Bacillus thuringiensis BMB171",
    "DEG1048": "Escherichia coli ST131 strain EC958",
    "DEG1049": "Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819",
    "DEG1050": "Campylobacter jejuni subsp. jejuni 81-176"
}


def run(query_file, deg_path, cv_path, cutoff, n_proc, out_path):
    """
    Get essentiality scores and classification for each protein in a
    file, according to the Geptop algorithm (prediction based on
    orthology and phylogeny) using the organisms in the Database of
    Essential Genes as reference.
    Args:
        query_file (str): query proteins in FASTA format.
        deg_path (str): directory containing DEG reference proteomes.
        cv_path (str): directory containing pre-computed composition
            vectors for DEG reference proteomes.
        cutoff (float): cutoff score for essentiality classification.
        n_proc (int): number of processes for parallelization.  # FIXME
        out_path (str): output directory.
    Returns:
        results (dict of str: (float, bool)): essentiality score and
            classification (True/False) for each protein.
    """
    blast_path = os.path.join(out_path, "blast_results")
    os.makedirs(blast_path, exist_ok=True)

    print("Reading reference proteome...")
    genes = SeqIO.parse(query_file, "fasta")
    scores = dict.fromkeys((gene.description for gene in genes), 0)
    make_blastdb(query_file)

    # Calculate query CV
    query_cv = composition_vector(query_file)

    print("Reading DEG database...")
    deg_organisms = [(os.path.splitext(os.path.basename(file))[0],
                      os.path.join(deg_path, file))
                     for file in os.listdir(deg_path)]

    print("Finding essential orthologs in:")
    for deg_id, deg_file in deg_organisms:
        t0 = time()
        print(
            "  %s %s" % (deg_id, ORG_NAMES[deg_id]),
            end=" ", flush=True
        )

        # Find ortholog pairs
        make_blastdb(deg_file)
        blast_all(query_file, deg_file, blast_path, n_proc)
        orthologs = rbh(
            os.path.join(blast_path, deg_id + "_f.xml"),
            os.path.join(blast_path, deg_id + "_r.xml")
        )
        os.remove(os.path.join(blast_path, deg_id + "_f.xml"))
        os.remove(os.path.join(blast_path, deg_id + "_r.xml"))
        remove_blastdb(deg_file)

        # Calculate species distance
        dist = get_distance(query_cv, deg_id, cv_path)
        if dist == 0:
            dist = 0.01

        # Get essentality score
        for pair in orthologs:
            scores[pair[0]] += is_essential(pair[1])/dist

        t1 = time()
        print("(%.2f s)" % (t1 - t0))

    os.rmdir(blast_path)
    remove_blastdb(query_file)
    print("Done. Normalizing essentiality scores...")
    scores = normalize(scores)

    # Classify genes
    results = {gene.split()[0]: (score, score > cutoff)
               for gene, score in scores.items()}
    return results


def composition_vector(species_fasta):
    """
    Calculate the composition vector of a given proteome, according to
    the algorithm by Qi et al. 2004 'Whole Proteome Prokaryote Phylogeny
    Without Sequence Alignment: A K-String Composition Approach'.
    Args:
        species_fasta (str): species proteome in FASTA format.
    Returns:
        cv (dict of int: float): composition vector.
    """
    freqs_k = Counter()  # K-words
    freqs_km1 = Counter()  # (K-1)-words
    freqs_km2 = Counter()  # (K-2)-words
    total_words = 0

    for protein in SeqIO.parse(species_fasta, "fasta"):
        if len(protein.seq) < CV_K:
            continue
        total_words += len(protein) - CV_K + 1  # total possible words
        protein_n = tuple(
            ord(c)-64 for c in protein.seq.upper()
            if ord(c)-64 in range(1, 27)
        )
        # Count words in the protein
        for kword in windowed(protein_n, CV_K):
            kword_n = word2num(kword)
            km1word_n = word2num(kword[:-1])
            km2word_n = word2num(kword[:-2])
            freqs_k[kword_n] += 1
            freqs_km1[km1word_n] += 1
            freqs_km2[km2word_n] += 1
        # Last (K-1)- and (K-2)-words
        freqs_km1[word2num(protein_n[-CV_K+1:])] += 1
        freqs_km2[word2num(protein_n[-CV_K+1:-1])] += 1
        freqs_km2[word2num(protein_n[-CV_K+2:])] += 1

    # Calculate relative frequencies
    probs_k = {w: f/total_words for w, f in freqs_k.items()}
    probs_km1 = {w: f/total_words for w, f in freqs_km1.items()}
    probs_km2 = {w: f/total_words for w, f in freqs_km2.items()}

    # Calculate predicted probabilities by Markov model
    probs_0 = dict()
    for kword_n in probs_k:
        # kword_n="ABCDE" -> a="ABCD"; b="BCDE"; c="BCD"
        a = word2num(num2word(kword_n)[:-1])
        b = word2num(num2word(kword_n)[1:])
        c = word2num(num2word(kword_n)[1:-1])
        probs_0[kword_n] = probs_km1[a] * probs_km1[b] / probs_km2[c]

    # Subtract random mutation background
    cv = dict()
    for kword_n in probs_k:
        if probs_0[kword_n] == 0:
            cv[kword_n] = 0
        else:
            cv[kword_n] = ((probs_k[kword_n] - probs_0[kword_n])
                           / probs_0[kword_n])
    return cv


def word2num(word):
    """Transform a letter string into a number in base 26."""
    num = 0
    for i, n in enumerate(word):
        num += n * (26**i)
    return num


def num2word(num):
    """Transform a number in base 26 into a letter string."""
    word = []
    while num:
        word.append(num % 26)
        num = num // 26
    return tuple(word)


def make_blastdb(seqs_file):
    """Make a BLAST database from a protein FASTA file.
    Args:
        seqs_file (str): input protein sequences in FASTA format.
    """
    makeblastdb = NcbimakeblastdbCommandline(
        dbtype="prot",
        input_file=seqs_file
    )
    out, err = makeblastdb()


def remove_blastdb(in_file):
    """Remove BLAST database files.
    Args:
        in_file (str): input file for construction of the database.
    """
    for db_file in glob(in_file + ".p*"):
        os.remove(db_file)


def blast_all(query, subject, xml_path, n_proc):
    """Perform a reciprocal all-vs-all BLAST search between two
    proteomes.
    Args:
        query (str): filename of first proteome in FASTA format.
        subject (str): filename of second proteome in FASTA format.
        xml_path (str): directory for BLAST output in XML format.
        n_proc (int): number of threads for BLAST search.
    """
    subj_id = os.path.splitext(os.path.basename(subject))[0]

    forward = NcbiblastpCommandline(
        query=query, db=subject,
        evalue=BLAST_EVALUE, outfmt=5, num_threads=n_proc,
        out=os.path.join(xml_path, "%s_f.xml" % subj_id)
    )
    out, err = forward()

    reverse = NcbiblastpCommandline(
        query=subject, db=query,
        evalue=BLAST_EVALUE, outfmt=5, num_threads=n_proc,
        out=os.path.join(xml_path, "%s_r.xml" % subj_id)
    )
    out, err = reverse()


def rbh(f_xml, r_xml):
    """Parse the output from a reciprocal all-vs-all BLAST search into a
    list of ortholog pairs.
    Args:
        f_xml (str): results of forward BLAST search in XML format.
        r_xml (str): results of reverse BLAST search in XML format.
    Returns:
        orthologs (list of (str, str)): ortholog pairs (as gene
            descriptions).
    """
    orthologs = []

    with open(f_xml) as f:
        records_f = list(NCBIXML.parse(f))
    with open(r_xml) as r:
        records_r = list(NCBIXML.parse(r))

    for record_f in records_f:
        if record_f.alignments:
            q = record_f.query  # query title
            s = record_f.alignments[0].hit_def  # best hit title
            pair_f = (q, s)
            for record_r in records_r:
                if record_r.alignments:
                    s = record_r.query
                    q = record_r.alignments[0].hit_def
                    pair_r = (q, s)
                    if pair_f == pair_r:
                        orthologs.append(pair_f)
                        break
    return orthologs


def get_distance(query_cv, ref_org, cv_path):
    """Calculate phylogenetic distance between two proteomes by the CV
    method.
    Args:
        query_cv (dict of int: float): composition vector of query
            proteome.
        ref_org (str): DEG id of subject proteome.
        cv_path (str): directory containing pre-computed composition
            vectors for DEG reference proteomes.
    Returns:
        dist (float): phylogenetic distance.
    """
    with open(os.path.join(cv_path, ref_org + ".json")) as f:
        ref_cv = {int(a): b for a, b in json.load(f).items()}

    dist = distance(query_cv, ref_cv)
    return dist


def distance(cv1, cv2):
    """Calculate distance between two composition vectors.
    Args:
        cv1 (dict of int: float): first composition vector.
        cv2 (dict of int: float): second composition vector.
    Returns:
        dist (float): distance.
    """
    # Calculate correlation between vectors as their cosine similarity
    a = sum(
        {kword: cv1[kword] * cv2[kword]
         for kword in cv1.keys() & cv2.keys()}.values()
    )
    b = sum(map(lambda x: x*x, cv1.values()))
    c = sum(map(lambda x: x*x, cv2.values()))
    corr = a / sqrt(b * c)

    # Get distance by normalizing to (0, 1)
    dist = (1 - corr) / 2
    return dist


def get_all_distances(query_org, ref_path, deg_path):  # FIXME
    reference_cvs = os.listdir(ref_path)
    for genome in reference_cvs:
        id_ = os.path.splitext(genome)[0]
        print(id_)
        t0 = time()
        d = get_distance(query_org, id_, deg_path)
        t1 = time()
        print(d, "%.2f s" % (t1 - t0))


def is_essential(hit_id):
    """Determine whether a DEG hit is essential (1) or not (0)."""
    deg_id = hit_id.split("|")[3]
    if "N" in deg_id:  # "DNEG" in id
        essential = 0
    else:  # "DEG" in id
        essential = 1
    return essential


def normalize(raw):
    """Normalize essentiality scores.
    Args:
        raw (dict of str: float): raw essentiality scores.
    Returns:
        norm (dict of str: float): normalized essentiality scores.
    """
    norm = {}
    s_max = max(raw.values())
    s_min = min(raw.values())
    for gene in raw:
        norm[gene] = (raw[gene] - s_min) / (s_max - s_min)
    return norm


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
    res = run(
        query_file="/home/jimena/Bartonella/CDSa/all.faa",
        deg_path="/home/jimena/Bartonella/DEGdb/deg_byorg/all",
        cv_path="/home/jimena/Bartonella/DEGdb/cv",
        cutoff=0.24,
        n_proc=4,
        out_path="/home/jimena/Escritorio"
    )
    for j, k in res.items():
        print(str(j) + "\t" + str(k))

    # t0 = time()
    # with open(os.path.join(REFSEQS_CV_DIR, "DEG1011.json"), "r") as fi:
    #     comp_vector1 = {int(a): b for a, b in json.load(fi).items()}
    # with open(os.path.join(REFSEQS_CV_DIR, "DEG1011.json"), "r") as fi:
    #     comp_vector2 = {int(a): b for a, b in json.load(fi).items()}
    # comp_vector1 = composition_vector("/home/jimena/Bartonella/DEGdb/deg_byorg/all/DEG1019.faa", 6)
    # comp_vector2 = composition_vector("/home/jimena/Bartonella/CDSa/all.faa", 6)
    # t1 = time()
    # d = distance(comp_vector1, comp_vector2)
    # t2 = time()
    # print(d, "%.2f s" % (t1-t0), "%.2f s" % (t2-t1))

    # t0 = time()
    # with open(os.path.join(REFSEQS_CV_DIR, "DEG1020.json"), "r") as fi:
    #     comp_vector1 = {int(a): b for a, b in json.load(fi).items()}
    # comp_vector1 = CompositionVector("/home/jimena/Bartonella/DEGdb/deg_byorg/all/DEG1019.faa")
    # comp_vector2 = CompositionVector("/home/jimena/Bartonella/CDSa/all.faa")
    # t1 = time()
    # d = Distance(comp_vector1, comp_vector2)
    # t2 = time()
    # print(d, "%.2f s" % (t1-t0), "%.2f s" % (t2-t1))
    # for x, y in comp_vector.items():
    #     print(x, y)
    # print(len(comp_vector))
