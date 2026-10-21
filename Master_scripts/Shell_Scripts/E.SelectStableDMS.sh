#!/bin/bash
#SBATCH --job-name=E.SelectStableDMS
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00

set -euo pipefail

printf 'Script: E.SelectStableDMS.sh\n'

script_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

if command -v module >/dev/null 2>&1; then
    module load R
fi

Rscript "${script_dir}/R_scripts/E.SelectStableDMS.R"
