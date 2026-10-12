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

master_scripts_dir="$1"
CH2_WORKDIR="$2"

if [[ ! -d "$master_scripts_dir" ]]; then
    echo "Master scripts directory does not exist: $master_scripts_dir" >&2
    exit 1
fi
if [[ ! -d "$CH2_WORKDIR" ]]; then
    echo "CH2_WORKDIR is not an existing directory: $CH2_WORKDIR" >&2
    exit 1
fi

master_scripts_dir="$(cd -- "$master_scripts_dir" && pwd -P)"
CH2_WORKDIR="$(cd -- "$CH2_WORKDIR" && pwd -P)"
CH2_SHELL_SCRIPTS_DIR="${master_scripts_dir}/Shell_Scripts"
a_script="${CH2_SHELL_SCRIPTS_DIR}/A.Generate_Splits.sh"
master_script="${master_scripts_dir}/masterSubmissionScript.sh"

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

if [[ ! -f "$a_script" || ! -f "$master_script" ]]; then
    echo "Could not find both pipeline submission scripts under ${master_scripts_dir}." >&2
    exit 1
fi

mkdir -p "${CH2_WORKDIR}/output/permanovaResults" "${CH2_WORKDIR}/output/MLdata"
export CH2_WORKDIR CH2_ITERATIONS CH2_SHELL_SCRIPTS_DIR

a_job_id="$(sbatch --parsable --export=ALL "$a_script")"
a_job_id="${a_job_id%%;*}"
if [[ ! "$a_job_id" =~ ^[0-9]+$ ]]; then
    echo "Could not parse Job A's Slurm job ID: ${a_job_id}" >&2
    exit 1
fi

array_job_id="$(sbatch --parsable --export=ALL \
    --array="$CH2_ITERATIONS" \
    --dependency="afterok:${a_job_id}" \
    "$master_script")"
 
printf 'Submitted Job A (ID %s).\n' "$a_job_id"
printf 'Submitted iteration array (ID %s), dependent on Job A; iterations: %s\n' \
    "$array_job_id" "$CH2_ITERATIONS"
