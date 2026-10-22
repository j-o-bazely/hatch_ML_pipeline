#!/bin/bash
#SBATCH --job-name=G.elasticNet_ML
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

printf 'Script: G.elasticNet_ML.sh\n'

script_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

## Activate virtual environment ##
source /data/home/bty500/.virtualenvs/hatch_ML/bin/activate

python3 "${script_dir}/R_scripts/G.elasticNet_ML.py"

deactivate