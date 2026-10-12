###          B2. Differential Methylation Analysis          ###

## Date: 06.10.2026

## Author: James B

## Description:
## This script runs methylKit DMA on training samples to generate
## % methylation difference values.



####  Environment Preparation  ####
## Set libraries ##
library(methylKit)

## Pull iteration & available cores ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
num_core <- as.numeric(args[2])

message(
  "Running Differential Methylation Analysis (DMA) for iteration ",
  iter,
  " using ",
  num_core,
  " cores."
)


## Set directories ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")

## Load methylBase object ##
hatch_training_unite_obj <- readRDS(
  file.path(tmpdir, paste0("hatch_trainingSamples_uniteObj_iter", iter, ".RDS"))
)

####  Calculate DMS ####
hatch_training_diff_meth <- calculateDiffMethDSS(
  meth = hatch_training_unite_obj,
  mc.cores = num_core
)

## Write to file ##
saveRDS(
  hatch_training_diff_meth,
  file.path(tmpdir, paste0("hatch_trainingSamples_diffMeth_iter", iter, ".RDS"))
)
