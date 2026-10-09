###          B2. Differential Methylation Analysis          ###

## Date: 06.10.2026

## Author: James B

## Description:
##    This script runs methylKit DMA on training samples to generate % methylation difference values    ##



####  Environment Preparation  ####  
## Set libraries ##
library(methylKit)

## Pull iteration & available cores ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
numCore <- as.numeric(args[2])

message("Running Differential Methylation Analysis (DMA) for iteration ", iter, " using ", numCore, " cores.")


## Set directories ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")

## Load methylBase object ##
hatch_trainingSamples_uniteObj <- readRDS(
  file.path(tmpdir, paste0("hatch_trainingSamples_uniteObj_iter", iter, ".RDS"))
)

####  Calculate DMS ####
hatch_trainingSamples_diffMeth <- calculateDiffMethDSS(
  meth = hatch_trainingSamples_uniteObj,
  mc.cores = numCore
)

saveRDS(
  hatch_trainingSamples_diffMeth,
  file.path(tmpdir, paste0("hatch_trainingSamples_diffMeth_iter", iter, ".RDS"))
)
