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

## Set seed for reproducbility ##
set.seed(123)


#### Identify samples to simulate model deployment ####

## Load metadata ##
meta <- read.csv(file.path(workdir, "input", "hatch_metadata.csv"))

## Select relevant metadata columns ##
meta_clean <- meta |>
  dplyr::select(
    Sample,
    Mother,
    maternal_infection_status,
    Amount_cloacal_p,
    parasiteFactor
  )

## Check sample IDs ##
if (anyDuplicated(meta_clean$Sample)) {
  stop("ERROR: Duplicate sample IDs detected in metadata.")
}

## Identify eligible mothers ##
mother_meta <- meta_clean |>
  dplyr::distinct(
    Mother,
    maternal_infection_status,
    Amount_cloacal_p
  )

## Select one infected mother from each specified group ##
infected_low <- mother_meta |>
  dplyr::filter(
    maternal_infection_status == "yes",
    Amount_cloacal_p == "1 to 5"
  ) |>
  dplyr::pull(Mother) |>
  sample(size = 1)

infected_high <- mother_meta |>
  dplyr::filter(
    maternal_infection_status == "yes",
    Amount_cloacal_p == ">26"
  ) |>
  dplyr::pull(Mother)

## Randomly select two uninfected mothers ##
uninfected_mothers <- mother_meta |>
  dplyr::filter(maternal_infection_status == "no") |>
  dplyr::pull(Mother) |>
  sample(size = 2, replace = FALSE)

## Combine selected mothers ##
deploy_mothers <- c(
  infected_low,
  infected_high,
  uninfected_mothers
)

## Hold out all associated hatchlings ##
deploy_samples <- meta_clean |>
  dplyr::filter(Mother %in% deploy_mothers)

## Check selected mothers and sample count ##
print(mother_meta |> dplyr::filter(Mother %in% deploy_mothers))
print(deploy_samples)
cat("Number of deployment samples:", nrow(deploy_samples), "\n")


## Write metadata to file ##
write.csv(
  deploy_samples,
  file.path(iterdir,
  "hatch_deploymentSamples_metadata.csv"),
  row.names = FALSE
)




####  Generate training and validation splits ####

## Remove deployment samples from metadata ##
meta_clean_no_deploy <- meta_clean %>%
  filter(! Mother %in% deploy_samples$Mother)

## Split mothers, then assign all associated samples to the same subset ##
mother_status <- meta_clean_no_deploy |>
  dplyr::distinct(Mother, maternal_infection_status)
if (anyDuplicated(mother_status$Mother)) {
  stop("ERROR: A mother has more than one maternal infection status.")
}

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
  ##  ... Partition mothers so families cannot cross subsets ##
  train_mother_index <- createDataPartition(
    mother_status$maternal_infection_status,
    p = 0.8,
    list = FALSE,
    times = 1
  )[, 1]

  training_mothers <- mother_status$Mother[train_mother_index]
  train <- meta_clean_no_deploy[
    meta_clean_no_deploy$Mother %in% training_mothers,
    ,
    drop = FALSE
  ]
  validation <- meta_clean_no_deploy[
    !(meta_clean_no_deploy$Mother %in% training_mothers),
    ,
    drop = FALSE
  ]

  if (length(intersect(train$Mother, validation$Mother)) > 0) {
    stop("ERROR: A mother was assigned to both training and validation.")
  }

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
