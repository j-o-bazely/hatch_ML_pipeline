#!/bin/bash
#SBATCH --job-name=submitPipeline
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00

# Submit the shared split-generation job, then the per-iteration job array.

set -euo pipefail

if [[ "$#" -ne 2 ]]; then
    echo "Usage: $0 <master-scripts-directory> <pipeline-work-directory>" >&2
    exit 2
fi

MASTER_SCRIPTS_DIR="$1"
CH2_WORKDIR="$2"

if [[ ! -d "$MASTER_SCRIPTS_DIR" ]]; then
    echo "Master scripts directory does not exist: $MASTER_SCRIPTS_DIR" >&2
    exit 1
fi
if [[ ! -d "$CH2_WORKDIR" ]]; then
    echo "CH2_WORKDIR is not an existing directory: $CH2_WORKDIR" >&2
    exit 1
fi

MASTER_SCRIPTS_DIR="$(cd -- "$MASTER_SCRIPTS_DIR" && pwd -P)"
CH2_WORKDIR="$(cd -- "$CH2_WORKDIR" && pwd -P)"
CH2_SHELL_SCRIPTS_DIR="${MASTER_SCRIPTS_DIR}/Shell_Scripts"
A_SCRIPT="${CH2_SHELL_SCRIPTS_DIR}/A.Generate_Splits.sh"
MASTER_SCRIPT="${MASTER_SCRIPTS_DIR}/masterSubmissionScript.sh"

if [[ ! -r "${CH2_WORKDIR}/input/hatch_metadata.csv" ]]; then
    echo "Cannot read required metadata: ${CH2_WORKDIR}/input/hatch_metadata.csv" >&2
    exit 1
fi
if [[ ! -r "${CH2_WORKDIR}/input/hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS" ]]; then
    echo "Cannot read required methylation object: ${CH2_WORKDIR}/input/hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS" >&2
    exit 1
fi

if [[ -z "${CH2_ITERATIONS:-}" ]]; then
    CH2_ITERATIONS=""
    for iteration in {1..1000}; do
        CH2_ITERATIONS+="${CH2_ITERATIONS:+,}${iteration}"
    done
fi

if [[ ! "$CH2_ITERATIONS" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; then
    echo "CH2_ITERATIONS must be a comma-separated list of positive iteration numbers (for example: 1,2,3)." >&2
    exit 1
fi

if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch was not found; run this launcher on a Slurm system." >&2
    exit 1
fi

if [[ ! -f "$A_SCRIPT" || ! -f "$MASTER_SCRIPT" ]]; then
    echo "Could not find both pipeline submission scripts under ${MASTER_SCRIPTS_DIR}." >&2
    exit 1
fi

mkdir -p "${CH2_WORKDIR}/output/permanovaResults" "${CH2_WORKDIR}/output/MLdata"
export CH2_WORKDIR CH2_ITERATIONS CH2_SHELL_SCRIPTS_DIR

A_JOB_ID="$(sbatch --parsable --export=ALL "$A_SCRIPT")"
A_JOB_ID="${A_JOB_ID%%;*}"
if [[ ! "$A_JOB_ID" =~ ^[0-9]+$ ]]; then
    echo "Could not parse Job A's Slurm job ID: ${A_JOB_ID}" >&2
    exit 1
fi

ARRAY_JOB_ID="$(sbatch --parsable --export=ALL \
    --array="$CH2_ITERATIONS" \
    --dependency="afterok:${A_JOB_ID}" \
    "$MASTER_SCRIPT")"
 
printf 'Submitted Job A (ID %s).\n' "$A_JOB_ID"
printf 'Submitted iteration array (ID %s), dependent on Job A; iterations: %s\n' \
    "$ARRAY_JOB_ID" "$CH2_ITERATIONS"
