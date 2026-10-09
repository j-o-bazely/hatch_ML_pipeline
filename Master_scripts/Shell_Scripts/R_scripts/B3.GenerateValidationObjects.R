###          B3. Generate Validation Objects          ###

## Date: 06.10.2026

## Author: James B

## Description:
##    This script generates the required analysis objects from the validation sample subset   ##



####  Environment Preparation ####
## Set libraries ##
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
iterdir <- file.path(tmpdir, "iterationData")

## load objects ##
hatch_validationSamples_metadata <- read.csv(
  file.path(iterdir, paste0("hatch_validationSamples_metadata_iter", iter, ".csv"))
)
hatch_allSamples_uniteObj <- readRDS(
  file.path(workdir, "input/hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS")
)


####  Generate validation object  ####
## Select only valudation samples ##
hatch_validationSamples_uniteObj <- reorganize(
  hatch_allSamples_uniteObj,
  sample.ids = hatch_validationSamples_metadata$Sample,
  treatment = hatch_validationSamples_metadata$parasiteFactor
)

## Save object to file ##
saveRDS(
  hatch_validationSamples_uniteObj,
  file.path(tmpdir, paste0("hatch_validationSamples_uniteObj_iter", iter, ".RDS"))
)

