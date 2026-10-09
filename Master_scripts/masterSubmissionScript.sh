#!/bin/bash
#SBATCH --job-name=masterSubmissionScript
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=240:00:00

set -euo pipefail

CH2_WORKDIR="${CH2_WORKDIR:?Set CH2_WORKDIR to the pipeline data directory before submitting this script}"
if [[ ! -d "$CH2_WORKDIR" ]]; then
    echo "CH2_WORKDIR is not an existing directory: $CH2_WORKDIR" >&2
    exit 1
fi
export CH2_WORKDIR
SHELL_SCRIPTS_DIR="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"

##  Define current iteration ##
iter="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required; submit this script as a job array}"

##  MethylKit Objects     ##
B1_JOB=$(sbatch --parsable \
    --export=ALL \
    "${SHELL_SCRIPTS_DIR}/B1.GenerateTrainingObjects.sh" "$iter")

B2_JOB=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${B1_JOB} \
    "${SHELL_SCRIPTS_DIR}/B2.MethylKitDMA.sh" "$iter")

B3_JOB=$(sbatch --parsable \
    --export=ALL \
    "${SHELL_SCRIPTS_DIR}/B3.GenerateValidationObjects.sh" "$iter")

##         PQLseq         ##
C1_JOB=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${B1_JOB} \
    "${SHELL_SCRIPTS_DIR}/C1.PrepPQLseq.sh" "$iter")

C2_JOB=$(sbatch --parsable \
    --export=ALL \
    --array=1-29 \
    --dependency=afterok:${C1_JOB} \
    "${SHELL_SCRIPTS_DIR}/C2.RunPQLseq.sh" "$iter")

C3_JOB=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${C2_JOB} \
    "${SHELL_SCRIPTS_DIR}/C3.PvalAdjustSLIM.sh" "$iter")

##      Combine Data      ##
D_JOB=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${B2_JOB}:${C3_JOB} \
    "${SHELL_SCRIPTS_DIR}/D.CombineDMAdata.sh" "$iter")

##       Assess DMS       ##
E_JOB=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${D_JOB}:${B3_JOB} \
    "${SHELL_SCRIPTS_DIR}/E.MeasureDMSinVal.sh" "$iter")

##          Tidy up       ##
F_JOB=$(sbatch --parsable \
    --export=ALL \
    --dependency=afterok:${E_JOB} \
    "${SHELL_SCRIPTS_DIR}/F.TidyEnvironment.sh" "$iter")