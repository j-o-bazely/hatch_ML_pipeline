## Description 

#### This pipeline iteratively generates a set of DMS in a specififed training subset
#### before then assessing their stability and application in ML models 
#### in a withheld validation subset, avoiding data leakage during the feature selection stage

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
requested iteration before submitting the job array.

Each C2 chromosome task retries its R script up to two times after the initial
attempt, waiting 30 seconds between attempts. If all three attempts fail, the
Slurm task exits with the final R script's failure status.