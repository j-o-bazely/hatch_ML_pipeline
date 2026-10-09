#!/bin/bash
#SBATCH --job-name=D.CombineDMAdata
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00

set -euo pipefail

####            D.CombineDMAdata.sh

SCRIPT_DIR="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

if command -v module >/dev/null 2>&1; then
  module load R
fi

## Pull current iteration ##
iter=$1

## Pull available cores ##
numCore="${SLURM_CPUS_PER_TASK:-1}"

## Run script to generate training and validation splits ##
Rscript "${SCRIPT_DIR}/R_scripts/D.CombineDMAdata.R" "$iter" "$numCore"
