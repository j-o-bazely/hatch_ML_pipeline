#!/bin/bash
#SBATCH --job-name=G.PrepareFinalModelData
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00

set -euo pipefail

script_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

if command -v module >/dev/null 2>&1; then
    module load R
fi

Rscript "${script_dir}/R_scripts/G.PrepareFinalModelData.R"
