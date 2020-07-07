#!/usr/bin/env python3
"""predict_essentiality.py

    Get predicted essentiality scores for all genes in a bacterial
    genome.

@author: Jimena Solana
"""

from argparse import ArgumentParser
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC_skew
import numpy as np
import pandas as pd
from joblib import load
from sklearn.pipeline import Pipeline

import geptop
import codonw


CODONW_FEATURES = ["-enc", "-gc", "-L_aa", "-hyd"]
FEATURES = ["strand_lead", "geptop", "Nc", "GC", "L_aa", "Gravy"]
DEG = os.path.join(os.path.dirname(__file__), "data/deg")
CV = os.path.join(os.path.dirname(__file__), "data/cv")
CLASSIFIER = os.path.join(os.path.dirname(__file__),
                          "classifier/classifier.joblib")


def find_ori_ter(gb):
    """Find positions of origin and terminus of replication in a genome,
    using the cumulative GC skew method.

    Args:
        gb (Bio.SeqRecord.SeqRecord): GenBank annotation.
    Returns:
        ori, ter (int): origin and terminus positions.
    """
    dna = gb.seq
    w = 100
    cum_gcskew = np.cumsum(GC_skew(dna, window=w))
    ori = np.argmin(cum_gcskew) * w
    ter = np.argmax(cum_gcskew) * w
    print("Found ori and ter positions: %d, %d" % (ori, ter))
    return ori, ter


def get_feature_table(gb, out_dir, ori, ter, geptop_params, codonw_features):
    """Calculate features for essentiality prediction for every gene.

    Args:
        gb (Bio.SeqRecord.SeqRecord): GenBank annotation.
        out_dir (str): directory for output files.
        ori (int): position of origin of replication.
        ter (int): position of terminus of replication.
        geptop_params (dict of str:str): parameters for essential
            ortholog mapping.
        codonw_features (list of str): list of CodonW features to
            calculate for each gene (as command-line arguments).
    Returns:
        feature_table (pd.DataFrame): table of results with each gene in
            a row and each calculated feature in a column.
    """
    # Extract proteome
    proteome_aa = os.path.join(out_dir, "proteome.faa")
    proteome_nt = os.path.join(out_dir, "proteome.fna")
    extract_cds(gb=gb, out_aa=proteome_aa, out_nt=proteome_nt)

    # Strand data
    feat_table = strand(proteome_aa, ori=ori, ter=ter)

    # Geptop scores
    geptop_results = geptop.run(query_file=proteome_aa, **geptop_params)
    for gene, score in geptop_results.items():
        feat_table.loc[gene, "geptop"] = score[0]

    # CodonW features
    codonw_results = codonw.run(proteome_nt, codonw_features)
    feat_table = feat_table.join(codonw_results)
    return feat_table


def extract_cds(gb, out_aa, out_nt):
    """Save CDS sequences as multiFASTA, both nucleotide and amino acid.

    For each GenBank feature annotated as a CDS with translation (a
    non-pseudogenised protein-coding gene), extract the nt and aa
    sequence and save each in a multiFASTA file, cumulatively. These
    files are created for downstream analyses.
    Args:
        gb (Bio.SeqRecord.SeqRecord): GenBank annotation.
        out_aa (str): protein sequence multiFASTA file.
        out_nt (str): nucleotide sequence multiFASTA file.
    """
    proteins_aa = []
    proteins_nt = []
    for feature in gb.features:
        if "essentiality" in feature.qualifiers:
            # Ignore already annotated genes
            continue
        if (feature.type == "CDS" and
                "translation" in feature.qualifiers and
                "locus_tag" in feature.qualifiers):
            protein_aa = SeqRecord(
                seq=Seq(feature.qualifiers["translation"][0]),
                id=feature.qualifiers["locus_tag"][0],
                description=make_description(feature)
            )
            protein_nt = SeqRecord(
                seq=feature.extract(gb.seq),
                id=feature.qualifiers["locus_tag"][0],
                description=make_description(feature)
            )
            proteins_aa.append(protein_aa)
            proteins_nt.append(protein_nt)
    SeqIO.write(proteins_aa, out_aa, "fasta")
    SeqIO.write(proteins_nt, out_nt, "fasta")


def make_description(feature):
    """Generate a FASTA header for a sequence."""
    if "product" in feature.qualifiers:
        product = feature.qualifiers["product"][0]
    else:
        product = "(no product)"
    location = "".join(str(part) for part in feature.location.parts)
    description = "%s %s" % (product, location)
    return description


def strand(proteome, ori, ter):
    """Determine strand directionality (leading or lagging) of each gene.

    Args:
        proteome (str): protein sequence multiFASTA file path.
        ori (int): position of origin of replication.
        ter (int): position of terminus of replication.
    Returns:
        strand_results (pd.DataFrame): table of results with each gene
            in a row and strand directionality in the only column
            (encoded as 1 = leading, 0 = lagging).
    """
    proteins = list(SeqIO.parse(proteome, "fasta"))
    strand_results = pd.DataFrame(
        index=[prot.id for prot in proteins],
        columns=["strand_lead"]
    )
    for prot in proteins:
        location = prot.description.split()[-1]
        start = int(location.split(":")[0][1:])
        sign = location[-2]
        if start > ori or start < ter:  # 0-180º replichore
            if sign == "+":  # leading strand
                strand_results.loc[prot.id, "strand_lead"] = 1
            else:  # lagging strand
                strand_results.loc[prot.id, "strand_lead"] = 0
        else:  # 180-360º replichore
            if sign == "+":  # lagging strand
                strand_results.loc[prot.id, "strand_lead"] = 0
            else:  # leading strand
                strand_results.loc[prot.id, "strand_lead"] = 1
    return strand_results


def get_essentiality_scores(feat_table, features, classifier):
    """
    # TODO
    Args:
        feat_table:
        features:
        classifier:

    Returns:

    """
    classifier = load(classifier)
    X_target = feat_table[features].values
    preprocess = Pipeline([
        ("imputation", classifier.named_steps["imputation"]),
        ("scaling", classifier.named_steps["scaling"])
    ])
    X_target = preprocess.fit(X_target).transform(X_target)  # Impute + scale
    y_probs = classifier.predict_proba(X_target)[:, 0]
    ess_scores = dict(zip(feat_table.index, y_probs))
    return ess_scores


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="predict-essentiality",
        description="Get predicted essentiality scores for all genes in a "
                    "bacterial genome"
    )
    parser.add_argument(
        "-g", dest="GB", required=True,
        help="original GenBank annotation file")
    parser.add_argument(
        "-o", dest="OUT_DIR", required=True,
        help="directory for output files")
    parser.add_argument(
        "-p1", dest="ORI", type=int,
        help="position of origin of replication")
    parser.add_argument(
        "-p2", dest="TER", type=int,
        help="position of terminus of replication")
    parser.add_argument(
        "-n", dest="NPROC", required=True, type=int,
        help="number of CPUs to use for execution")
    args, unknown = parser.parse_known_args()
    GENBANK = args.GB
    OUT_DIR = args.OUT_DIR
    os.makedirs(OUT_DIR, exist_ok=True)
    OUT_FEATTABLE = os.path.join(OUT_DIR, "feature_table.csv")
    OUT_ESSTABLE = os.path.join(OUT_DIR, "essentiality_table.csv")
    genbank_id = os.path.splitext(os.path.basename(GENBANK))[0]
    GENBANK_M1 = os.path.join(OUT_DIR, genbank_id + ".gbm1")
    NPROC = args.NPROC
    ORI = args.ORI
    TER = args.TER
    GEPTOP_PARAMS = {
        "deg_dir": DEG,
        "cv_dir": CV,
        "n_proc": NPROC,
        "out_dir": OUT_DIR
    }

    # Check input
    try:
        annotation = SeqIO.read(GENBANK, "genbank")
    except (FileNotFoundError, ValueError):
        raise SystemExit("\n\terror: could not read annotation file\n")
    if ORI and TER:
        if not all([coord in range(1, len(annotation)+1)
                    for coord in (ORI, TER)]):
            raise SystemExit("\n\terror: invalid ori/ter coordinates\n")

    # Get ori + ter coordinates
    if not (ORI and TER):
        ORI, TER = find_ori_ter(annotation)

    # Get table of all gene features
    print("Extracting gene features...")
    feature_table = get_feature_table(annotation, OUT_DIR, ORI, TER,
                                      GEPTOP_PARAMS, CODONW_FEATURES)
    feature_table.to_csv(OUT_FEATTABLE)

    # Load classifier and get essentiality scores
    print("Calculating essentiality scores...")
    essentiality_scores = get_essentiality_scores(feature_table,
                                                  FEATURES, CLASSIFIER)
    essentiality_table = pd.DataFrame(list(essentiality_scores.items()),
                                      columns=["locus_tag", "ess_score"])
    essentiality_table.to_csv(OUT_ESSTABLE, index=False)

    # Create modified-I GenBank file
    for gene in annotation.features:
        if (gene.type == "CDS" and
                "translation" in gene.qualifiers and
                "locus_tag" in gene.qualifiers):
            if "essentiality" in gene.qualifiers:
                # Ignore already annotated genes
                continue
            locus_tag = gene.qualifiers["locus_tag"][0]
            gene.qualifiers["essentiality"] = essentiality_scores[locus_tag]
        elif gene.type in ("tRNA", "rRNA", "tmRNA", "ncRNA"):
            gene.qualifiers["essentiality"] = 1
        elif gene.type == "CDS":  # pseudo-gene
            gene.qualifiers["essentiality"] = 0
    SeqIO.write(annotation, GENBANK_M1, "genbank")
    print(
        "Done. Results in %s, %s and %s."
        % (OUT_FEATTABLE, OUT_ESSTABLE, GENBANK_M1)
    )
