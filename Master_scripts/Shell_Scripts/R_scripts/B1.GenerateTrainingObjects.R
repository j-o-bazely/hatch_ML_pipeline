###          B1. Generate Training Object               ###

## Date: 06.10.2026

## Author: James B

## Description:
##  This script generates the required analysis objects
##  from the training sample subset




####  Environment Preparation  ####
## Load required packages
library(methylKit)

## Pull iteration & available cores ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
num_core <- as.numeric(args[2])

## Set directories ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")
iterdir <- file.path(workdir, "tmp", "iterationData")


####  Training Data Preparation  ####
## Load metadata for training samples of the current iteration
hatch_training_meta <- read.csv(
  file.path(
    iterdir,
    paste0("hatch_trainingSamples_metadata_iter", iter, ".csv")
  )
)

## Load full hatchling methylBase object ##
hatch_all_samples_unite_obj <- readRDS(
  file.path(
    workdir,
    "input/hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS"
  )
)

## Check: Are all samples in full methylBase object ##
if (length(setdiff(
  hatch_training_meta$Sample,
  hatch_all_samples_unite_obj@sample.ids
)) > 0) {
  stop("ERROR: Training samples missing from full object.")
}

## Use reorganize() to keep only training samples ##
hatch_training_unite_obj <- reorganize(
  hatch_all_samples_unite_obj,
  sample.ids = hatch_training_meta$Sample,
  treatment = hatch_training_meta$parasiteFactor
)

## Check: Has sample <-> treatment integrity been maintained ##
stopifnot(
  identical(
    hatch_training_unite_obj@sample.ids,
    hatch_training_meta$Sample
  ),
  identical(
    hatch_training_unite_obj@treatment,
    hatch_training_meta$parasiteFactor
  )
)

## Save object to file ##
saveRDS(
  hatch_training_unite_obj,
  file.path(tmpdir, paste0("hatch_trainingSamples_uniteObj_iter", iter, ".RDS"))
)
