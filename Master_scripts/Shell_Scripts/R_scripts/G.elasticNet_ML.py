#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def make_grouped_splits(X, y, groups, n_splits, seed, description):
    splitter = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=seed,
    )
    splits = list(splitter.split(X, y, groups))
    for fold, (train_indices, test_indices) in enumerate(splits, start=1):
        if set(groups.iloc[train_indices]) & set(groups.iloc[test_indices]):
            raise ValueError(
                f"{description} fold {fold} assigns a mother to both subsets."
            )
        if set(y.iloc[train_indices]) != {0, 1}:
            raise ValueError(
                f"{description} fold {fold} training data must contain both classes."
            )
        if set(y.iloc[test_indices]) != {0, 1}:
            raise ValueError(
                f"{description} fold {fold} held-out data must contain both classes."
            )
    return splits


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("iteration", type=int)
    args = parser.parse_args()
    if args.iteration < 1:
        parser.error("iteration must be a positive integer")

    workdir = Path(
        os.environ.get("CH2_WORKDIR", Path.cwd())
    ).expanduser().resolve()
    output_dir = workdir / "output" / "MLresults"
    ml_dir = workdir / "output" / "MLdata"
    iteration = args.iteration

    training_data_file = (
        ml_dir
        / f"hatch_trainingSamples_DMS_percMeth_df_iter{iteration}.RDS"
    )
    training_metadata_file = (
        ml_dir / f"hatch_trainingSamples_metadata_iter{iteration}.csv"
    )

    X = rdata.read_rds(training_data_file).T
    training_metadata = pd.read_csv(
        training_metadata_file,
        dtype={"Sample": str, "Mother": str},
    )
    required_metadata_columns = {"Sample", "Mother", "parasiteFactor"}
    missing_metadata_columns = required_metadata_columns - set(
        training_metadata.columns
    )
    if missing_metadata_columns:
        raise ValueError(
            "Training metadata is missing required columns: "
            f"{sorted(missing_metadata_columns)}."
        )
    training_metadata = training_metadata.set_index("Sample")
    y = training_metadata["parasiteFactor"]
    groups = training_metadata["Mother"]

    if not X.index.equals(y.index):
        raise ValueError(
            "Training sample order does not match between data and metadata."
        )
    if not X.index.is_unique or not y.index.is_unique:
        raise ValueError("Training sample IDs must be unique.")
    if y.isna().any():
        raise ValueError("Training response labels must not be missing.")
    if groups.isna().any():
        raise ValueError("Mother identifiers must not be missing.")
    y = pd.to_numeric(y)
    if set(y.unique()) != {0, 1}:
        raise ValueError("The response variable must contain both 0 and 1.")
    mothers_per_class = {
        label: groups[y.eq(label)].nunique()
        for label in (0, 1)
    }
    if min(mothers_per_class.values()) < 5:
        raise ValueError(
            "At least five distinct mothers must contribute samples to each "
            "class for 5-fold family-aware cross-validation; got "
            f"{mothers_per_class}."
        )

    base_pipeline = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    class_weight="balanced",
                    max_iter=10000,
                    random_state=42,
                ),
            ),
        ]
    )
    param_grid = {
        "model__C": np.logspace(-4, 4, 10),
        "model__l1_ratio": [0.1, 0.25, 0.5, 0.75, 0.9],
    }
    n_jobs = max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))

    outer_splits = make_grouped_splits(
        X, y, groups, n_splits=5, seed=42, description="Outer CV"
    )
    oof_probabilities = np.full(len(y), np.nan)
    for fold, (train_indices, test_indices) in enumerate(outer_splits, start=1):
        X_train, y_train = X.iloc[train_indices], y.iloc[train_indices]
        groups_train = groups.iloc[train_indices]
        inner_splits = make_grouped_splits(
            X_train,
            y_train,
            groups_train,
            n_splits=5,
            seed=42 + fold,
            description=f"Inner CV for outer fold {fold}",
        )
        fold_search = GridSearchCV(
            estimator=base_pipeline,
            param_grid=param_grid,
            scoring="roc_auc",
            cv=inner_splits,
            n_jobs=n_jobs,
            refit=True,
        )
        fold_search.fit(X_train, y_train)
        oof_probabilities[test_indices] = fold_search.predict_proba(
            X.iloc[test_indices]
        )[:, 1]

    if np.isnan(oof_probabilities).any():
        raise RuntimeError("Family-aware outer CV did not predict every sample.")

    thresholds = np.arange(0.1, 0.91, 0.01)
    threshold_scores = {
        threshold: balanced_accuracy_score(
            y, (oof_probabilities >= threshold).astype(int)
        )
        for threshold in thresholds
    }
    optimal_threshold = max(threshold_scores, key=threshold_scores.get)

    final_splits = make_grouped_splits(
        X, y, groups, n_splits=5, seed=42, description="Final model CV"
    )
    final_search = GridSearchCV(
        estimator=base_pipeline,
        param_grid=param_grid,
        scoring="roc_auc",
        cv=final_splits,
        n_jobs=n_jobs,
        refit=True,
    )
    final_search.fit(X, y)
    result = pd.DataFrame(
        [
            {
                "iteration": iteration,
                "best_C": final_search.best_params_["model__C"],
                "best_l1_ratio": final_search.best_params_["model__l1_ratio"],
                "optimal_threshold": optimal_threshold,
            }
        ]
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results_file = output_dir / f"iteration_model_results_iter{iteration}.csv"
    result.to_csv(results_file, index=False)
    print(f"Iteration results saved to: {results_file}")


if __name__ == "__main__":
    main()
