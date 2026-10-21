#!/bin/bash
#SBATCH --job-name=G.elasticNet_ML
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

if [[ "$#" -ne 1 ]]; then
    echo "Usage: $0 <iteration>" >&2
    exit 2
fi

iter="$1"
printf 'Script: G.elasticNet_ML.sh\n'
printf 'Iteration: %s\n' "$iter"

script_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

python3 "${script_dir}/R_scripts/G.elasticNet_ML.py" "$iter"
