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
    b3_job="$(sed -n '2p' "$manifest")"
    if [[ -z "$d_job" || -z "$b3_job" ]]; then
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

for iteration in "${iterations[@]}"; do
    manifest="${manifest_dir}/iter${iteration}.jobs"
    d_job="$(sed -n '1p' "$manifest")"
    b3_job="$(sed -n '2p' "$manifest")"

    e_job="$(sbatch --parsable \
        --export=ALL \
        --dependency="afterok:${d_job}:${b3_job}:${stable_job}" \
        "${shell_scripts_dir}/F.MeasureDMSinVal.sh" "$iteration")"

    f_job="$(sbatch --parsable \
        --export=ALL \
        --dependency="afterok:${e_job}" \
        "${shell_scripts_dir}/G.elasticNet_ML.sh" "$iteration")"

    sbatch --parsable \
        --export=ALL \
        --dependency="afterok:${f_job}" \
        "${shell_scripts_dir}/H.TidyEnvironment.sh" "$iteration" >/dev/null
done

printf 'Submitted validation, ML, and tidy-up jobs for iterations: %s\n' \
    "$CH2_ITERATIONS"
