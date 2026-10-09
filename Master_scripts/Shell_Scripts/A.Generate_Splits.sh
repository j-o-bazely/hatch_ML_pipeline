#!/bin/bash
#SBATCH --job-name=A.Generate_Splits
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail

####            A.Generate_Splits.sh

SCRIPT_DIR="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

if command -v module >/dev/null 2>&1; then
  module load R
fi

## Run script to generate training and validation splits ##
Rscript "${SCRIPT_DIR}/R_scripts/A.Generate_Splits.R"