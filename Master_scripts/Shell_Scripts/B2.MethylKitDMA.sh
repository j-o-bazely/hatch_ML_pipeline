#!/bin/bash
#SBATCH --job-name=B2.MethylKitDMA
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

####            B2.MethylKitDMA.sh

script_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

if command -v module >/dev/null 2>&1; then
  module load R
fi

## Pull current iteration ##
iter=$1

## Pull available cores ##
num_core="${SLURM_CPUS_PER_TASK:-1}"

## Run script to generate training and validation splits ##
Rscript "${script_dir}/R_scripts/B2.MethylKitDMA.R" "$iter" "$num_core"
