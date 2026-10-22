## Description

The automated pipeline identifies consensus DMS using repeated, family-aware
resamples of the fixed model-training cohort. It then generates feature
matrices at those DMS for training, threshold calibration, and deployment.
Machine-learning analysis is intentionally interactive and is not submitted
automatically.

## Cohort design

`A.Generate_Splits.R` expects 29 maternal groups with four hatchlings each and
creates a single, reproducible split:

- 19 mothers (76 hatchlings) for DMS identification and model development.
- 6 mothers (24 hatchlings) for threshold calibration, with 3 infected and
  3 uninfected mothers.
- 4 mothers (16 hatchlings) for deployment evaluation.

The deployment cohort contains one randomly selected low-burden infected
mother, one high-burden infected mother, and two uninfected mothers. Calibration
mothers are selected from the remaining groups to give an equal infected and
uninfected mix; burden is not used for this selection. The other 19 mothers
form the model-training cohort. All cohort assignments remain fixed across
iterations, and all hatchlings from a mother remain together.

Each DMS iteration resamples 80% of the 19 model-training mothers. DMS are
identified only within the model-training cohort. Stage E selects consensus
DMS present in more than 95% of requested iterations; failed or missing
iterations count as absent.

## Running the automated pipeline

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
per-iteration training-resample metadata exist before submitting the job array.
After the DMS iterations finish, the coordinator selects consensus DMS and
generates the three final feature matrices. The automated workflow ends there;
it does not fit a model, select a threshold, or evaluate deployment samples.

Each C2 chromosome task retries its R script up to two times after the initial
attempt, waiting 30 seconds between attempts. If all three attempts fail, the
Slurm task exits with the final R script's failure status.

## Interactive machine learning

Open `Shell_Scripts/R_scripts/G.elasticNet_ML.ipynb` in Jupyter after the
automated pipeline completes. Set `CH2_WORKDIR` in the notebook environment,
then run the cells to load and inspect the generated cohort data, tune and fit
an example elastic-net model, select a threshold using the calibration cohort,
save the frozen model, and evaluate it on the held-out deployment cohort. The
notebook is an interactive starting point; model choice and analysis can be
edited there. Do not use deployment results to revise the model or threshold.

The notebook environment needs pandas, NumPy, joblib, and scikit-learn. Use
scikit-learn 1.1 or later for `StratifiedGroupKFold`.

## Generated data

The main automated outputs are:

- `output/MLresults/stable_DMS_gt95pct.csv`: consensus DMS selected using only
  training resamples.
- `output/MLdata/hatch_modelTraining_DMS_percMeth_final.csv` and matching
  `hatch_modelTraining_metadata_final.csv`.
- `output/MLdata/hatch_thresholdCalibration_DMS_percMeth_final.csv` and
  matching `hatch_thresholdCalibration_metadata_final.csv`.
- `output/MLdata/hatch_deployment_DMS_percMeth_final.csv` and matching
  `hatch_deployment_metadata_final.csv`.

Feature matrices are CSV files with DMS as rows and hatchling samples as
columns. Matching metadata files provide sample IDs, maternal groups, and
response labels. These files are the hand-off point from automated processing
to the interactive notebook.
