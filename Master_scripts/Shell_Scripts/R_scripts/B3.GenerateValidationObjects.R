###          B3. Generate Validation Objects          ###

## Date: 06.10.2026

## Author: James B

## Description:
## This script generates the required analysis objects from the
## validation sample subset.



####  Environment Preparation ####
## Set libraries ##
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
iterdir <- file.path(tmpdir, "iterationData")

## load objects ##
hatch_validation_meta <- read.csv(
  file.path(
    iterdir,
    paste0("hatch_validationSamples_metadata_iter", iter, ".csv")
  )
)
hatch_all_samples_unite_obj <- readRDS(
  file.path(
    workdir,
    "input/hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS"
  )
)


####  Generate validation object  ####
## Select only valudation samples ##
hatch_validation_unite_obj <- reorganize(
  hatch_all_samples_unite_obj,
  sample.ids = hatch_validation_meta$Sample,
  treatment = hatch_validation_meta$parasiteFactor
)

## Save object to file ##
saveRDS(
  hatch_validation_unite_obj,
  file.path(
    tmpdir,
    paste0("hatch_validationSamples_uniteObj_iter", iter, ".RDS")
  )
)
