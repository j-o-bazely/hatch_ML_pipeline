## Description 

#### This pipeline generates DMS on family-aware training subsets, assesses them
#### in the corresponding validation subsets, and selects stable ML features
#### across iterations. A deployment cohort is held out from all iterations.

## Running the pipeline

Generate the training and validation split metadata manually before submitting
the remaining pipeline. Set `CH2_WORKDIR` to the pipeline data directory, which
must contain `input/hatch_metadata.csv` and
`input/hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS`.
Run the following commands from the `Master_scripts` directory.

For the default 1,000 iterations, run:

    export CH2_WORKDIR=/path/to/pipeline-data
    Rscript Shell_Scripts/R_scripts/A.Generate_Splits.R
    sbatch submitPipeline.sh . "$CH2_WORKDIR"

To generate a selected set of iterations, set `CH2_ITERATIONS` when running the
script and use the same value when submitting the pipeline:

    export CH2_ITERATIONS=1,2,3
    Rscript Shell_Scripts/R_scripts/A.Generate_Splits.R
    sbatch submitPipeline.sh . "$CH2_WORKDIR"

If `CH2_ITERATIONS` is not specified, the submission launcher schedules
iterations 1 through 1,000. It checks that the split metadata exists for every
requested iteration before submitting the job array. The launcher first
submits the per-iteration DMS jobs, then a coordinator schedules the global
stable-DMS selection and remaining pipeline stages.

Each C2 chromosome task retries its R script up to two times after the initial
attempt, waiting 30 seconds between attempts. If all three attempts fail, the
Slurm task exits with the final R script's failure status.

Once all per-iteration DMS analysis jobs finish, stage E selects DMS found in
more than 95% of the requested iterations. Failed or missing DMS iterations
count as absent in this frequency calculation. Stage F uses this common set
when writing the ML feature matrices for each iteration.

Stage G runs `Shell_Scripts/R_scripts/G.elasticNet_ML.py` once per iteration.
It tunes the elastic-net model using the training samples, selects a
classification threshold from nested out-of-fold training predictions, and
writes the selected hyperparameters and threshold to
`output/MLresults/iteration_model_results_iter<iteration>.csv`. The G job uses
Python 3 with pandas, NumPy, rdata, and scikit-learn available in its
environment. Stage H tidies intermediate files after the ML job succeeds.
Both inner hyperparameter tuning and outer out-of-fold prediction use
stratified, mother-grouped folds; the script requires at least five distinct
mothers contributing to each response class. Use scikit-learn 1.1 or later
for `StratifiedGroupKFold`.

Because the stable-DMS set is calculated across all requested iterations, a
per-iteration validation subset may be part of another iteration's DMS
selection. Treat the separate deployment cohort, which is excluded from all
iterations, as the independent held-out set.

Training and validation splits are made by mother: all samples from the same
mother are assigned to the same subset within each iteration. The 80/20 split
is therefore approximate at the sample level when family sizes differ.