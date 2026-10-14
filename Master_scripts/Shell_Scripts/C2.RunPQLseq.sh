#!/bin/bash
#SBATCH --job-name=C2.runPQLseq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1:29

set -euo pipefail

iter=$1

printf 'Script: %s\n' "$(basename -- "$0")"
printf 'Iteration: %s\n' "$iter"

####            C2.runPQLseq.sh

script_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"

if command -v module >/dev/null 2>&1; then
  module load R
fi

## Pull available cores ##
num_core="${SLURM_CPUS_PER_TASK:-1}"

## Pull chromosome ##
chr_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for chromosome array tasks}"

## Retry this chromosome task when R exits unsuccessfully.
max_attempts=3
retry_delay=30
attempt=1
while (( attempt <= max_attempts )); do
  if Rscript "${script_dir}/R_scripts/C2.RunPQLseq.R" "$iter" "$num_core" "$chr_id"; then
    exit 0
  else
    status=$?
  fi

  printf 'C2 failed for iteration %s, chromosome task %s (attempt %s/%s; exit status %s).\n' \
    "$iter" "$chr_id" "$attempt" "$max_attempts" "$status" >&2

  if (( attempt == max_attempts )); then
    exit "$status"
  fi

  sleep "$retry_delay"
  attempt=$((attempt + 1))
done
