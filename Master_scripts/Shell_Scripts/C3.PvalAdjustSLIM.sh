#!/bin/bash
#SBATCH --job-name=C3.pValAdjustSLIM
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00

set -euo pipefail

iter=$1

printf 'Script: C3.PvalAdjustSLIM.sh\n'
printf 'Iteration: %s\n' "$iter"

####            C3.PvalAdjustSLIM.sh

script_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

if command -v module >/dev/null 2>&1; then
  module load R
fi

## Pull available cores ##
num_core="${SLURM_CPUS_PER_TASK:-1}"

## Run script to calculate adjusted p values ##
Rscript "${script_dir}/R_scripts/C3.PvalAdjustSLIM.R" "$iter" "$num_core"
