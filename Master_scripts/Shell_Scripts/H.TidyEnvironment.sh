#!/bin/bash
#SBATCH --job-name=H.TidyEnvironment
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00

printf 'Script: H.TidyEnvironment.sh\n'

#### Clean per-iteration intermediates after final feature matrices are ready ####

: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"
: "${CH2_ITERATIONS:?CH2_ITERATIONS is required; submit with submitPipeline.sh}"
if [[ "$#" -ne 1 || ! "$1" =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: $0 <DMS-iteration-array-job-id>" >&2
    exit 2
fi

array_job_id="$1"
tmpdir="${CH2_WORKDIR}/tmp"
iterdir="${tmpdir}/iterationData"
manifest_dir="${tmpdir}/pipeline_jobs_${array_job_id}"

if [[ ! -r "${CH2_WORKDIR}/output/MLresults/stable_DMS_ge50pct.csv" ]]; then
    echo "Consensus DMS output is missing; refusing to clean intermediates." >&2
    exit 1
fi
for cohort in modelTraining thresholdCalibration deployment; do
    if [[ ! -r "${CH2_WORKDIR}/output/MLdata/hatch_${cohort}_DMS_percMeth_final.csv" ||
          ! -r "${CH2_WORKDIR}/output/MLdata/hatch_${cohort}_metadata_final.csv" ]]; then
        echo "Final ${cohort} feature matrix or metadata is missing; refusing to clean intermediates." >&2
        exit 1
    fi
done

IFS=, read -r -a iterations <<< "$CH2_ITERATIONS"
for iteration in "${iterations[@]}"; do
    if [[ ! "$iteration" =~ ^[1-9][0-9]*$ ]]; then
        echo "CH2_ITERATIONS contains an invalid iteration: ${iteration}" >&2
        exit 1
    fi
done

for iteration in "${iterations[@]}"; do
    for intermediate_file in "${tmpdir}"/*iter"${iteration}"[_.]*; do
        [[ -e "$intermediate_file" ]] || continue
        filename="${intermediate_file##*/}"
        if [[ "$filename" == "hatch_trainingSamples_PQLseqOutput_DMS_q0.01_iter${iteration}.RDS" ]]; then
            continue
        fi
        rm -f -- "$intermediate_file"
    done
    rm -f -- \
        "${iterdir}/hatch_trainingSamples_metadata_iter${iteration}.csv" \
        "${iterdir}/hatch_validationSamples_metadata_iter${iteration}.csv"
done

if [[ -d "$manifest_dir" ]]; then
    rm -f -- "${manifest_dir}"/iter*.jobs
    rmdir -- "$manifest_dir"
fi

printf 'Removed per-iteration temporary files for iterations: %s\n' \
    "$CH2_ITERATIONS"
printf 'Retained q0.01 DMS files, fixed cohort metadata, consensus DMS, and final feature matrices.\n'