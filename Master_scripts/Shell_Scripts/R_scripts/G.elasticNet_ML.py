#!/usr/bin/env python3

import os
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import rdata
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def load_cohort(ml_dir, cohort_name):
    feature_file = ml_dir / f"hatch_{cohort_name}_DMS_percMeth_final.RDS"
    metadata_file = ml_dir / f"hatch_{cohort_name}_metadata_final.csv"

    features = rdata.read_rds(feature_file)
    if not isinstance(features, pd.DataFrame):
        raise TypeError(f"Expected an R data frame in {feature_file}.")
    X = features.T
    metadata = pd.read_csv(
        metadata_file,
        dtype={"Sample": str, "Mother": str},
    )
    required_columns = {"Sample", "Mother", "parasiteFactor"}
    missing_columns = required_columns - set(metadata.columns)
    if missing_columns:
        raise ValueError(
            f"{cohort_name} metadata is missing columns: "
            f"{sorted(missing_columns)}."
        )
    metadata = metadata.set_index("Sample")

    if not X.index.equals(metadata.index):
        raise ValueError(
            f"{cohort_name} sample order does not match feature data."
        )
    if not X.index.is_unique or not X.columns.is_unique:
        raise ValueError(
            f"{cohort_name} sample IDs and feature IDs must be unique."
        )
    if metadata["Mother"].isna().any():
        raise ValueError(f"{cohort_name} maternal group IDs must not be missing.")

    y = pd.to_numeric(metadata["parasiteFactor"], errors="raise")
    if y.isna().any() or set(y.unique()) != {0, 1}:
        raise ValueError(
            f"{cohort_name} response labels must contain both 0 and 1."
        )
    X = X.apply(pd.to_numeric, errors="raise")
    feature_values = X.to_numpy(dtype=float)
    if (
        np.isinf(feature_values).any()
        or np.isnan(feature_values).all(axis=0).any()
    ):
        raise ValueError(
            f"{cohort_name} feature data contain infinite values or an "
            "entirely missing feature."
        )
    if X.shape[1] == 0:
        raise ValueError(f"{cohort_name} feature matrix has no DMS features.")
    return X, y.astype(int), metadata


def make_grouped_splits(X, y, groups, n_splits, seed):
    splitter = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=seed,
    )
    splits = list(splitter.split(X, y, groups))
    for fold, (train_indices, test_indices) in enumerate(splits, start=1):
        if set(groups.iloc[train_indices]) & set(groups.iloc[test_indices]):
            raise ValueError(f"CV fold {fold} assigns a mother to both subsets.")
        if set(y.iloc[train_indices]) != {0, 1}:
            raise ValueError(f"CV fold {fold} training subset lacks a class.")
        if set(y.iloc[test_indices]) != {0, 1}:
            raise ValueError(f"CV fold {fold} held-out subset lacks a class.")
    return splits


def choose_threshold(y_true, probabilities):
    candidates = np.unique(np.concatenate(([0.0, 0.5, 1.0], probabilities)))
    scores = np.array(
        [
            balanced_accuracy_score(y_true, probabilities >= threshold)
            for threshold in candidates
        ]
    )
    best_score = scores.max()
    best_indices = np.flatnonzero(np.isclose(scores, best_score))
    best_index = min(
        best_indices,
        key=lambda index: (abs(candidates[index] - 0.5), candidates[index]),
    )
    return float(candidates[best_index]), float(best_score)


def main():
    workdir = Path(
        os.environ.get("CH2_WORKDIR", Path.cwd())
    ).expanduser().resolve()
    ml_dir = workdir / "output" / "MLdata"
    output_dir = workdir / "output" / "MLresults"

    X_train, y_train, train_meta = load_cohort(ml_dir, "modelTraining")
    X_calibration, y_calibration, calibration_meta = load_cohort(
        ml_dir,
        "thresholdCalibration",
    )
    X_deployment, y_deployment, deployment_meta = load_cohort(
        ml_dir,
        "deployment",
    )

    if (
        not X_train.columns.equals(X_calibration.columns)
        or not X_train.columns.equals(X_deployment.columns)
    ):
        raise ValueError("Feature IDs or order differ between final cohorts.")
    if (
        set(train_meta["Mother"]) & set(calibration_meta["Mother"])
        or set(train_meta["Mother"]) & set(deployment_meta["Mother"])
        or set(calibration_meta["Mother"]) & set(deployment_meta["Mother"])
    ):
        raise ValueError("Maternal groups must not overlap between cohorts.")

    mothers_per_class = {
        label: train_meta.loc[y_train.eq(label), "Mother"].nunique()
        for label in (0, 1)
    }
    if min(mothers_per_class.values()) < 5:
        raise ValueError(
            "At least five distinct training mothers must contribute samples "
            "to each response class for 5-fold family-aware cross-validation; "
            f"got {mothers_per_class}."
        )

    pipeline = Pipeline(
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
    cv_splits = make_grouped_splits(
        X_train,
        y_train,
        train_meta["Mother"],
        n_splits=5,
        seed=42,
    )
    search = GridSearchCV(
        estimator=pipeline,
        param_grid=param_grid,
        scoring="roc_auc",
        cv=cv_splits,
        n_jobs=max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))),
        refit=True,
    )
    search.fit(X_train, y_train)

    calibration_probabilities = search.best_estimator_.predict_proba(
        X_calibration
    )[:, 1]
    optimal_threshold, calibration_balanced_accuracy = choose_threshold(
        y_calibration.to_numpy(),
        calibration_probabilities,
    )

    deployment_probabilities = search.best_estimator_.predict_proba(
        X_deployment
    )[:, 1]
    deployment_predictions = (
        deployment_probabilities >= optimal_threshold
    ).astype(int)
    deployment_balanced_accuracy = balanced_accuracy_score(
        y_deployment,
        deployment_predictions,
    )
    deployment_auc = roc_auc_score(y_deployment, deployment_probabilities)
    deployment_accuracy = accuracy_score(
        y_deployment,
        deployment_predictions,
    )
    tn, fp, fn, tp = confusion_matrix(
        y_deployment,
        deployment_predictions,
        labels=[0, 1],
    ).ravel()

    output_dir.mkdir(parents=True, exist_ok=True)
    model_file = output_dir / "final_elastic_net_model.joblib"
    joblib.dump(
        {
            "estimator": search.best_estimator_,
            "threshold": optimal_threshold,
            "feature_ids": X_train.columns.tolist(),
        },
        model_file,
    )

    results = pd.DataFrame(
        [
            {
                "training_samples": len(y_train),
                "training_mothers": train_meta["Mother"].nunique(),
                "calibration_samples": len(y_calibration),
                "calibration_mothers": calibration_meta["Mother"].nunique(),
                "deployment_samples": len(y_deployment),
                "deployment_mothers": deployment_meta["Mother"].nunique(),
                "consensus_dms": X_train.shape[1],
                "best_C": search.best_params_["model__C"],
                "best_l1_ratio": search.best_params_["model__l1_ratio"],
                "cv_best_roc_auc": search.best_score_,
                "optimal_threshold": optimal_threshold,
                "calibration_balanced_accuracy": calibration_balanced_accuracy,
                "deployment_balanced_accuracy": deployment_balanced_accuracy,
                "deployment_accuracy": deployment_accuracy,
                "deployment_roc_auc": deployment_auc,
                "deployment_true_negative": tn,
                "deployment_false_positive": fp,
                "deployment_false_negative": fn,
                "deployment_true_positive": tp,
                "model_file": model_file.name,
            }
        ]
    )
    results_file = output_dir / "final_model_results.csv"
    results.to_csv(results_file, index=False)

    calibration_predictions = pd.DataFrame(
        {
            "Sample": calibration_meta.index,
            "Mother": calibration_meta["Mother"].to_numpy(),
            "observed_parasiteFactor": y_calibration.to_numpy(),
            "predicted_probability": calibration_probabilities,
            "predicted_parasiteFactor": (
                calibration_probabilities >= optimal_threshold
            ).astype(int),
        }
    )
    calibration_predictions.to_csv(
        output_dir / "threshold_calibration_predictions.csv",
        index=False,
    )

    deployment_predictions_df = pd.DataFrame(
        {
            "Sample": deployment_meta.index,
            "Mother": deployment_meta["Mother"].to_numpy(),
            "observed_parasiteFactor": y_deployment.to_numpy(),
            "predicted_probability": deployment_probabilities,
            "predicted_parasiteFactor": deployment_predictions,
        }
    )
    deployment_predictions_df.to_csv(
        output_dir / "deployment_predictions.csv",
        index=False,
    )

    print(f"Final model saved to: {model_file}")
    print(f"Threshold-calibration threshold: {optimal_threshold:.6g}")
    print(f"Final model results saved to: {results_file}")


if __name__ == "__main__":
    main()
