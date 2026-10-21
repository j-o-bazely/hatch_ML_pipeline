#!/bin/bash
#SBATCH --job-name=DMSReadyBarrier
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00

set -euo pipefail

printf 'All dependent DMS jobs have completed.\n'
