###          B1. Generate Training Object               ###

## Date: 06.10.2026

## Author: James B

## Description:
##    This script generates the required analysis objects from the training sample subset   ##




####  Environment Preparation  ####
## Load required packages
library(methylKit)

## Pull iteration & available cores ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
numCore <- as.numeric(args[2])

## Set directories ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")
iterdir <- file.path(workdir, "tmp", "iterationData")


####  Training Data Preparation  ####
## Load metadata for training samples of the current iteration
hatch_trainingSamples_metadata <- read.csv(
  file.path(iterdir, paste0("hatch_trainingSamples_metadata_iter", iter, ".csv"))
)

## Load full hatchling methylBase object ##
hatch_allSamples_uniteObj <- readRDS(
  file.path(workdir, "input/hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS")
)

## Use reorganize() to keep only training samples ##
hatch_trainingSamples_uniteObj <- reorganize(
  hatch_allSamples_uniteObj,
  sample.ids = hatch_trainingSamples_metadata$Sample,
  treatment = hatch_trainingSamples_metadata$parasiteFactor
)

## Save object to file ##
saveRDS(
  hatch_trainingSamples_uniteObj,
  file.path(tmpdir, paste0("hatch_trainingSamples_uniteObj_iter", iter, ".RDS"))
)

