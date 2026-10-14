#!/bin/bash
#SBATCH --job-name=masterSubmissionScript
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=240:00:00

set -euo pipefail

##  Define current iteration ##
iter="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required; submit this script as a job array}"

printf 'Script: masterSubmissionScript.sh\n'
printf 'Iteration: %s\n' "$iter"

CH2_WORKDIR="${CH2_WORKDIR:?Set CH2_WORKDIR to the pipeline data directory before submitting this script}"
if [[ ! -d "$CH2_WORKDIR" ]]; then
    echo "CH2_WORKDIR is not an existing directory: $CH2_WORKDIR" >&2
    exit 1
fi
export CH2_WORKDIR
shell_scripts_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"

##  MethylKit Objects     ##
b1_job=$(sbatch --parsable \
    --export=ALL \
    "${shell_scripts_dir}/B1.GenerateTrainingObjects.sh" "$iter")

b2_job=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${b1_job} \
    "${shell_scripts_dir}/B2.MethylKitDMA.sh" "$iter")

b3_job=$(sbatch --parsable \
    --export=ALL \
    "${shell_scripts_dir}/B3.GenerateValidationObjects.sh" "$iter")

##         PQLseq         ##
c1_job=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${b1_job} \
    "${shell_scripts_dir}/C1.PrepPQLseq.sh" "$iter")

c2_job=$(sbatch --parsable \
    --export=ALL \
    --array=1-29 \
    --dependency=afterok:${c1_job} \
    "${shell_scripts_dir}/C2.RunPQLseq.sh" "$iter")

c3_job=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${c2_job} \
    "${shell_scripts_dir}/C3.PvalAdjustSLIM.sh" "$iter")

##      Combine Data      ##
d_job=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${b2_job}:${c3_job} \
    "${shell_scripts_dir}/D.CombineDMAdata.sh" "$iter")

##       Assess DMS       ##
e_job=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${d_job}:${b3_job} \
    "${shell_scripts_dir}/E.MeasureDMSinVal.sh" "$iter")

##          Tidy up       ##
f_job=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${e_job} \
    "${shell_scripts_dir}/F.TidyEnvironment.sh" "$iter")