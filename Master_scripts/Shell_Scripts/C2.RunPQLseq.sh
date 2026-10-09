#!/bin/bash
#SBATCH --job-name=C2.runPQLseq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1:29

set -euo pipefail

####            C2.runPQLseq.sh

SCRIPT_DIR="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

if command -v module >/dev/null 2>&1; then
  module load R
fi

## Pull current iteration ##
iter=$1

## Pull available cores ##
numCore="${SLURM_CPUS_PER_TASK:-1}"

## Pull chromosome ##
chr_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for chromosome array tasks}"

## Run script to generate training and validation splits ##
Rscript "${SCRIPT_DIR}/R_scripts/C2.RunPQLseq.R" "$iter" "$numCore" "$chr_id"
