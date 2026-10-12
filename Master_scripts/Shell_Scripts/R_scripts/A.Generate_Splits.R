###          A. Generate Training & Validation Splits               ###

## Date: 06.10.2026

## Author: James B

## Description:
##  This script dictates the samples
##  which will make up the training and validation subsets
##  for each iteration



####  Environment Preparation  ####
library(caret)
library(dplyr)

## Set path to working directory
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
iterdir <- file.path(workdir, "tmp", "iterationData")
dir.create(iterdir, recursive = TRUE, showWarnings = FALSE)

#### Generate training and validation splits ####
##  Load metadata ##
meta <- read.csv(paste0(workdir, "/input/hatch_metadata.csv"))

##  Pull sample name + Mother + Maternal infection status ##
meta_clean <- meta |>
  dplyr::select(
    Sample,
    Mother,
    maternal_infection_status,
    parasiteFactor
  )

## Check: Are samples unique? ##
if (anyDuplicated(meta_clean$Sample)) {
  stop("ERROR: Duplicate sample IDs detected in metadata.")
}

##  Generate a series of splits of 80%/20% ##

## Set seed for reproducbility ##
set.seed(123)

## Define number of iterations to generate ##
iterations <- as.integer(
  strsplit(
    Sys.getenv(
      "CH2_ITERATIONS",
      unset = paste(1:1000, collapse = ",")
    ),
    ","
  )[[1]]
)

## For each iteration specified...
for (iteration in iterations) {
  ##  ... Define training samples for that iteration ##
  train_index <- createDataPartition(
    meta_clean$maternal_infection_status,
    p = 0.8,
    list = FALSE,
    times = 1
  )[, 1]

  ## ... Pull metadata for training samples ##
  train <- meta_clean[train_index, , drop = FALSE]

  ## ... Pull metadata for validation samples ##
  validation <- meta_clean[-train_index, , drop = FALSE]


  ## Write metadata for each iteration to file ##
  write.csv(
    train,
    file.path(iterdir,
              paste0("hatch_trainingSamples_metadata_iter", iteration, ".csv")),
    row.names = FALSE
  )

  write.csv(
    validation,
    file.path(
      iterdir,
      paste0("hatch_validationSamples_metadata_iter", iteration, ".csv")
    ),
    row.names = FALSE
  )
}
