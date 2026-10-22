## Description

The pipeline uses fixed maternal-group cohorts for model development,
classification-threshold calibration, and a final proof-of-concept deployment
evaluation. All hatchlings from a mother remain in the same cohort.

## Cohort design

`A.Generate_Splits.R` expects 29 maternal groups with four hatchlings each and
creates a single, reproducible split:

- 20 mothers (80 hatchlings) for DMS identification and model development.
- 6 mothers (24 hatchlings) for classification-threshold calibration.
- 3 mothers (12 hatchlings) for final deployment evaluation.

The deployment cohort contains one randomly selected low-burden infected
mother, one high-burden infected mother, and one uninfected mother. The six
calibration mothers are selected from the remaining groups, stratified by
maternal infection status; the other 20 form the model-training cohort. The
three cohort assignments remain fixed across all iterations.

Each DMS iteration resamples 80% of the 20 model-training mothers, keeping all
of each selected mother's hatchlings together. DMS are therefore identified
only within the training cohort. Stage E selects consensus DMS present in more
than 95% of requested iterations; failed or missing iterations count as absent.

After consensus DMS selection, the pipeline constructs feature matrices for
all three fixed cohorts. The final elastic-net model is tuned by 5-fold
stratified, mother-grouped cross-validation on the 20 training mothers and is
then refit on all 20. The six calibration mothers are used to select a
classification threshold by balanced accuracy. The fitted model and threshold
are frozen before predictions are evaluated on the three deployment mothers.
The deployment cohort is not used for feature selection, hyperparameter
tuning, or threshold selection.

The deployment evaluation contains only three independent maternal groups.
It is intended as a proof-of-concept and should not be interpreted as a precise
estimate of performance in future populations. Cross-validation tunes
hyperparameters after the consensus DMS set has been selected on the training
cohort; its score is not presented as an unbiased estimate of the complete
feature-selection-and-modeling procedure.

## Running the pipeline

Set `CH2_WORKDIR` to the pipeline data directory, which must contain
`input/hatch_metadata.csv` and
`input/hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS`.
Run these commands from the `Master_scripts` directory.

For the default 1,000 DMS iterations, run:

    export CH2_WORKDIR=/path/to/pipeline-data
    Rscript Shell_Scripts/R_scripts/A.Generate_Splits.R
    sbatch submitPipeline.sh . "$CH2_WORKDIR"

To generate a selected set of DMS iterations, set `CH2_ITERATIONS` when
generating the splits and use the same value when submitting the pipeline:

    export CH2_ITERATIONS=1,2,3
    Rscript Shell_Scripts/R_scripts/A.Generate_Splits.R
    sbatch submitPipeline.sh . "$CH2_WORKDIR"

If `CH2_ITERATIONS` is not specified, the submission launcher schedules
iterations 1 through 1,000. It checks that fixed cohort metadata and
per-iteration resample metadata exist before submitting the job array. The
launcher submits the per-iteration DMS jobs, then a coordinator schedules
consensus-DMS selection, final feature-matrix preparation, model fitting,
threshold calibration, deployment evaluation, and intermediate cleanup.

Each C2 chromosome task retries its R script up to two times after the initial
attempt, waiting 30 seconds between attempts. If all three attempts fail, the
Slurm task exits with the final R script's failure status.

Final outputs are written to `output/MLresults/`:

- `stable_DMS_gt95pct.csv`: consensus DMS selected from training-only
  iterations.
- `final_elastic_net_model.joblib`: fitted scikit-learn estimator, threshold,
  and feature IDs.
- `final_model_results.csv`: selected hyperparameters, threshold, and
  calibration/deployment metrics.
- `threshold_calibration_predictions.csv`: predictions for the six calibration
  mothers.
- `deployment_predictions.csv`: predictions and observed labels for the three
  deployment mothers.

The final-model job uses Python 3 with pandas, NumPy, joblib, rdata, and
scikit-learn available in its environment. Use scikit-learn 1.1 or later for
`StratifiedGroupKFold`.
