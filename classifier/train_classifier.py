#!/usr/bin/env python3
"""train_classifier.py

    Train a logistic regression classifier with gene feature data from
    reference organisms.

@author: Jimena Solana
"""

import os
import warnings

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import LeaveOneGroupOut
from joblib import dump


DATA_DIR = os.path.join(os.path.dirname(__file__), "features")
ID_MAP = os.path.join(DATA_DIR, "locustag2deg.csv")
FEATURES = ["strand_lead", "geptop", "Nc", "GC", "L_aa", "Gravy"]


if __name__ == "__main__":
    warnings.filterwarnings("ignore")

    # Data preprocessing
    print("Reading data...")
    data_files = [f for f in os.listdir(DATA_DIR)
                  if f.endswith("_feature_table.csv")]
    data_files.sort()
    all_genes = pd.DataFrame()
    for file in data_files:
        df = pd.read_csv(os.path.join(DATA_DIR, file), index_col=0)
        df["organism"] = file.split("_")[0]
        all_genes = pd.concat([all_genes, df])
    all_genes = all_genes[~all_genes.index.duplicated(keep="last")]
    locustag2deg = pd.read_csv(ID_MAP,
                               usecols=[1, 2], index_col="locus_tag",
                               na_values="-")
    locustag2deg = locustag2deg[~locustag2deg.index.duplicated(keep="last")]
    all_genes = all_genes.join(locustag2deg, how="inner")
    deg_genes = all_genes.dropna(subset=["deg_id"])
    deg_genes["essential"] = deg_genes["deg_id"].map(
        lambda x: "non-essential" if "N" in x else "essential"
    )
    X = deg_genes[FEATURES].values  # feature array
    y = deg_genes["essential"].values  # target array
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        test_size=0.4,
                                                        random_state=123)

    # Hyperparameter tuning
    print("Tuning classifier...")
    params = {
        "logistic_regression__C": np.logspace(-4, 4),
    }
    search_pipeline = Pipeline([
        ("imputation", SimpleImputer()),  # impute missing data (Nc values)
        ("scaling", MinMaxScaler()),  # scale data to [0, 1]
        ("logistic_regression", LogisticRegression(solver="saga",
                                                   class_weight="balanced",
                                                   random_state=123))
    ])
    search = GridSearchCV(search_pipeline, params, scoring="roc_auc", n_jobs=-1)
    search.fit(X_train, y_train)
    best_C = search.best_params_["logistic_regression__C"]
    pipeline = Pipeline([
        ("imputation", SimpleImputer()),
        ("scaling", MinMaxScaler()),
        ("logistic_regression", LogisticRegression(C=best_C,
                                                   solver="saga",
                                                   class_weight="balanced",
                                                   random_state=123))
    ])
    pipeline.fit(X_train, y_train)

    # Model evaluation
    y_pred_prob = pipeline.predict_proba(X_test)[:, 1]
    print("  Classifer ROC AUC: %.4f" % roc_auc_score(y_test, y_pred_prob))
    cv_scores = cross_val_score(pipeline, X, y,
                                cv=5, scoring="roc_auc", n_jobs=-1)
    print("  ROC AUC with 5-fold cross-validation: %.4f +- %.4f"
          % (cv_scores.mean(), cv_scores.std()))
    le = LabelEncoder().fit(deg_genes["organism"].unique())
    groups = le.transform(deg_genes["organism"])
    cross_validator = LeaveOneGroupOut().get_n_splits(groups=groups)
    cv_scores2 = cross_val_score(pipeline, X, y,
                                 cv=cross_validator, scoring="roc_auc",
                                 n_jobs=-1)
    print("  ROC AUC with leave-one-species-out cross-validation: %.4f +- %.4f"
          % (cv_scores2.mean(), cv_scores2.std()))

    # Save classifier
    print("Saving classifier...")
    dump(pipeline, "./classifier_retrained.joblib")
