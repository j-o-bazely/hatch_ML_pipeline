###          A. Generate Training & Validation Splits               ###

## Date: 06.10.2026

## Author: James B

## Description:
##    This script dictates the samples which will make up the training and validation subsets for each iteration    ##



####  Environment Preparation  ####
library(caret)
library(dplyr)
Sys.getenv()

## Set path to working directory
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
iterdir <- file.path(workdir, "tmp", "iterationData")
dir.create(iterdir, recursive = TRUE, showWarnings = FALSE)

#### Generate training and validation splits ####
##  Load metadata ##
meta <- read.csv(paste0(workdir, '/input/hatch_metadata.csv'))

##  Pull sample name + Mother + Maternal infection status ##
meta_clean <- meta %>%
  dplyr::select(
    Sample,
    Mother,
    maternal_infection_status,
    parasiteFactor
  )

##  Generate a series of splits of 80%/20% ##
set.seed(123)
iterations <- as.integer(
  strsplit(Sys.getenv("CH2_ITERATIONS", unset = paste(1:100, collapse = ",")), ",")[[1]]
)

for (iteration in iterations) {
  train_index <- createDataPartition(
    meta_clean$maternal_infection_status,
    p = 0.8,
    list = FALSE,
    times = 1
  )[, 1]
  
  train <- meta_clean[train_index, , drop = FALSE]
  validation <- meta_clean[-train_index, , drop = FALSE]
  
  write.csv(
    train,
    file.path(iterdir, paste0("hatch_trainingSamples_metadata_iter", iteration, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    validation,
    file.path(iterdir, paste0("hatch_validationSamples_metadata_iter", iteration, ".csv")),
    row.names = FALSE
  )
}
