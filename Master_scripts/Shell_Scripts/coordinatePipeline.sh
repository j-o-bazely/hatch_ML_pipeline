#!/bin/bash
#SBATCH --job-name=coordinatePipeline
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

set -euo pipefail

if [[ "$#" -ne 1 ]]; then
    echo "Usage: $0 <DMS-iteration-array-job-id>" >&2
    exit 2
fi

master_array_job_id="$1"
: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"
shell_scripts_dir="${CH2_SHELL_SCRIPTS_DIR:?CH2_SHELL_SCRIPTS_DIR is required; submit with submitPipeline.sh}"
if [[ -z "${CH2_ITERATIONS:-}" ]]; then
    echo "CH2_ITERATIONS is required; submit with submitPipeline.sh" >&2
    exit 1
fi

IFS=, read -r -a iterations <<< "$CH2_ITERATIONS"
manifest_dir="${CH2_WORKDIR}/tmp/pipeline_jobs_${master_array_job_id}"
d_jobs=()

for iteration in "${iterations[@]}"; do
    manifest="${manifest_dir}/iter${iteration}.jobs"
    if [[ ! -r "$manifest" ]]; then
        echo "Missing per-iteration job manifest: ${manifest}" >&2
        exit 1
    fi
    d_job="$(sed -n '1p' "$manifest")"
    if [[ -z "$d_job" ]]; then
        echo "Invalid per-iteration job manifest: ${manifest}" >&2
        exit 1
    fi
    d_jobs+=("$d_job")
done

barrier_jobs=()
for ((offset = 0; offset < ${#d_jobs[@]}; offset += 50)); do
    batch=("${d_jobs[@]:offset:50}")
    dependency_ids="$(printf '%s:' "${batch[@]}")"
    dependency_ids="${dependency_ids%:}"
    barrier_job="$(sbatch --parsable \
        --export=ALL \
        --dependency="afterany:${dependency_ids}" \
        "${shell_scripts_dir}/DMSReadyBarrier.sh")"
    barrier_jobs+=("$barrier_job")
done

dependency_ids="$(printf '%s:' "${barrier_jobs[@]}")"
dependency_ids="${dependency_ids%:}"
stable_job="$(sbatch --parsable \
    --export=ALL \
    --dependency="afterok:${dependency_ids}" \
    "${shell_scripts_dir}/E.SelectStableDMS.sh")"
printf 'Submitted stable-DMS selection (ID %s).\n' "$stable_job"

final_data_job="$(sbatch --parsable \
    --export=ALL \
    --dependency="afterok:${stable_job}" \
    "${shell_scripts_dir}/G.PrepareFinalModelData.sh")"
# The stable-DMS job depends on barriers covering every per-iteration DMS job.
# Keeping cleanup behind feature generation therefore waits for all iterations'
# upstream jobs before removing their temporary files.
cleanup_job="$(sbatch --parsable \
    --export=ALL \
    --dependency="afterok:${final_data_job}" \
    "${shell_scripts_dir}/H.TidyEnvironment.sh" "$master_array_job_id")"
printf 'Submitted three-cohort feature-matrix preparation (ID %s).\n' \
    "$final_data_job"
printf 'Submitted post-generation cleanup (ID %s).\n' "$cleanup_job"
printf 'Automation ends after the training, threshold-calibration, and deployment matrices are generated.\n'
