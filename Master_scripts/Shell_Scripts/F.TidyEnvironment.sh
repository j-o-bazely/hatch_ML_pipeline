#!/bin/bash
#SBATCH --job-name=F.TidyEnvironment
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00

####            F.TidyEnvironment.sh

: "${CH2_WORKDIR:?CH2_WORKDIR is required; submit with submitPipeline.sh}"
workdir="$CH2_WORKDIR"
tempdir="${workdir}/tmp/"

## Pull current iteration ##
iter=$1

## Check that expected outputs for this iteration have been generated   ##
if [ "$(find "${workdir}/output/permanovaResults/iter${iter}_permanovaResults.txt" -print 2>/dev/null | wc -l)" -eq 0 ]; then
    echo "Expected permanova results for iteration ${iter} not found. Exiting script."
    exit 1
fi

if [ "$(find "${workdir}/output/MLdata/" -name "*iter${iter}*" -print 2>/dev/null | wc -l)" -ne 6 ]; then
    echo "Expected ML data for iteration ${iter} not found. Exiting script."
    exit 1
fi

## Delete intermediate files generated for this iteration ##
rm -f "${tempdir}"*iter${iter}[_.]*