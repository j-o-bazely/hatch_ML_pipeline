###          C1. PQLseq object preparation          ###

## Date: 06.10.2026
## Author:
## Description:
## This script prepares training sample data for PQLseq


#### Environment Preparation ####

library(methylKit)
library(dplyr)

## Pull iteration & available cores
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
num_core <- as.numeric(args[2])

## Set paths
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")
iterdir <- file.path(tmpdir, "iterationData")


## Define functions to generate number of Cs and coverage data for PQLseq input

prep_num_cs <- function(list_of_dfs) {
  for (i in seq_along(list_of_dfs)) {
    list_of_dfs[[i]] <- list_of_dfs[[i]] |>
      dplyr::select(-chr)

    col_names <- colnames(list_of_dfs[[i]])
    cols_to_delete <- grep("^coverage", col_names)

    if (length(cols_to_delete) > 0) {
      list_of_dfs[[i]] <- list_of_dfs[[i]][
        , -cols_to_delete, drop = FALSE
      ]
    }

    list_of_dfs[[i]] <- list_of_dfs[[i]] |>
      rename_with(~gsub("^numCs", "", .x))

    colnames(list_of_dfs[[i]]) <- paste0(
      "ind",
      colnames(list_of_dfs[[i]])
    )
  }

  return(list_of_dfs) # nolint: return_linter.
}


prep_coverage <- function(list_of_dfs) {
  for (i in seq_along(list_of_dfs)) {
    list_of_dfs[[i]] <- list_of_dfs[[i]] |>
      dplyr::select(-chr)

    col_names <- colnames(list_of_dfs[[i]])
    cols_to_delete <- grep("^numCs", col_names)

    if (length(cols_to_delete) > 0) {
      list_of_dfs[[i]] <- list_of_dfs[[i]][
        , -cols_to_delete, drop = FALSE
      ]
    }

    list_of_dfs[[i]] <- list_of_dfs[[i]] |>
      rename_with(~gsub("^coverage", "", .x))

    colnames(list_of_dfs[[i]]) <- paste0(
      "ind",
      colnames(list_of_dfs[[i]])
    )
  }

  return(list_of_dfs) # nolint: return_linter.
}


#### Load Input Objects ####

## Load training sample uniteCov object
hatch_training_unite_obj <- readRDS(
  file.path(
    tmpdir,
    paste0("hatch_trainingSamples_uniteObj_iter", iter, ".RDS")
  )
)

## Load metadata for training samples
hatch_training_meta <- read.csv(
  file.path(
    iterdir,
    paste0("hatch_trainingSamples_metadata_iter", iter, ".csv")
  )
)


#### Check Metadata Integrity ####

## Check: Sample IDs are unique and non-missing
if (
  anyNA(hatch_training_meta$Sample) ||
    anyDuplicated(hatch_training_meta$Sample)
) {
  stop("ERROR: Missing or duplicated training sample IDs.")
}

message("PASS: Training sample IDs are unique and non-missing.")


## Check: Maternal group and treatment assignments are non-missing
if (
  anyNA(hatch_training_meta$Mother) ||
    anyNA(hatch_training_meta$parasiteFactor)
) {
  stop("ERROR: Missing maternal group or treatment assignments.")
}

message("PASS: Maternal group and treatment assignments are non-missing.")


## Check: Metadata matches the methylation object
if (
  !identical(
    hatch_training_unite_obj@sample.ids,
    hatch_training_meta$Sample
  ) ||
    !identical(
      hatch_training_unite_obj@treatment,
      hatch_training_meta$parasiteFactor
    )
) {
  stop(
    "ERROR: Training metadata does not match the methylation object ",
    "in sample order or treatment assignment."
  )
}

message("PASS: Training metadata matches the methylation object.")


#### Generate Kinship Matrix ####

## Build the kinship matrix by maternal group
mother_ids <- as.character(hatch_training_meta$Mother)
same_mother <- outer(mother_ids, mother_ids, FUN = "==")
same_mother[is.na(same_mother)] <- FALSE

hatch_gen_mat <- 0.5 * same_mother
diag(hatch_gen_mat) <- 1

dimnames(hatch_gen_mat) <- list(
  hatch_training_meta$Sample,
  hatch_training_meta$Sample
)

hatch_gen_df <- as.data.frame(hatch_gen_mat)


## Check: Kinship matrix integrity
if (
  !identical(rownames(hatch_gen_mat), hatch_training_meta$Sample) ||
    !identical(colnames(hatch_gen_mat), hatch_training_meta$Sample) ||
    !identical(hatch_gen_mat, t(hatch_gen_mat)) ||
    !all(diag(hatch_gen_mat) == 1)
) {
  stop(
    "ERROR: Kinship matrix is asymmetric, misaligned with sample IDs, ",
    "or has an invalid diagonal."
  )
}

message("PASS: Kinship matrix integrity checks.")


## Save kinship matrix
saveRDS(
  hatch_gen_df,
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqInput_KinshipMatrix_iter",
      iter,
      ".RDS"
    )
  )
)


#### Prepare Methylation Data Inputs for PQLseq ####

## Convert uniteCov object to data.frame
hatch_training_samples_data <- as.data.frame(
  getData(hatch_training_unite_obj)
)


## Check: Required methylation columns exist
sample_numbers <- seq_along(hatch_training_unite_obj@sample.ids)

required_cols <- c(
  paste0("numCs", sample_numbers),
  paste0("numTs", sample_numbers),
  paste0("coverage", sample_numbers)
)

if (!all(required_cols %in% colnames(hatch_training_samples_data))) {
  stop(
    "ERROR: Methylation input columns do not match the expected ",
    "sample count and naming convention."
  )
}

message("PASS: Required methylation input columns found.")


## Create genomic site names
site_names <- paste(
  hatch_training_samples_data$chr,
  hatch_training_samples_data$start,
  sep = "."
)

## Remove start, end and strand columns; keep chr
hatch_training_samples_data <- hatch_training_samples_data |>
  dplyr::select(-start, -end, -strand)

## Assign site names after column selection
rownames(hatch_training_samples_data) <- site_names


## Delete all numTs columns
col_names <- colnames(hatch_training_samples_data)
cols_to_delete <- grep("^numTs", col_names)

if (length(cols_to_delete) > 0) {
  hatch_training_samples_data <- hatch_training_samples_data[
    , -cols_to_delete, drop = FALSE
  ]
}


## Split data into chromosomes
hatch_training_samples_data$chr <- as.factor(
  hatch_training_samples_data$chr
)

hatch_training_data_list <- split(
  hatch_training_samples_data,
  hatch_training_samples_data$chr
)


## Generate number-of-Cs inputs
training_data_by_chr_num_cs <- prep_num_cs(
  hatch_training_data_list
)


## Generate coverage inputs
training_data_by_chr_coverage <- prep_coverage(
  hatch_training_data_list
)


## Check: Chromosomes, genomic sites and sample columns match
if (
  !identical(
    names(training_data_by_chr_num_cs),
    names(training_data_by_chr_coverage)
  ) ||
    !all(vapply(
      names(training_data_by_chr_num_cs),
      function(chr) {
        identical(
          rownames(training_data_by_chr_num_cs[[chr]]),
          rownames(training_data_by_chr_coverage[[chr]])
        ) &&
          identical(
            colnames(training_data_by_chr_num_cs[[chr]]),
            colnames(training_data_by_chr_coverage[[chr]])
          )
      },
      logical(1)
    ))
) {
  stop(
    "ERROR: Number-of-Cs and coverage inputs have mismatched ",
    "chromosomes, genomic sites or sample columns."
  )
}

message(
  "PASS: Chromosome, genomic site and sample-column integrity checks."
)


## Save number-of-Cs inputs per chromosome
lapply(
  names(training_data_by_chr_num_cs),
  function(chr_name) {
    saveRDS(
      training_data_by_chr_num_cs[[chr_name]],
      file.path(
        tmpdir,
        paste0(
          "hatch_trainingSamples_PQLseqInput_NumCs_iter",
          iter,
          "_",
          chr_name,
          ".RDS"
        )
      )
    )
  }
)


## Save coverage inputs per chromosome
lapply(
  names(training_data_by_chr_coverage),
  function(chr_name) {
    saveRDS(
      training_data_by_chr_coverage[[chr_name]],
      file.path(
        tmpdir,
        paste0(
          "hatch_trainingSamples_PQLseqInput_Coverage_iter",
          iter,
          "_",
          chr_name,
          ".RDS"
        )
      )
    )
  }
)


#### Generate Treatment Vector ####

treatment_no <- setNames(
  as.integer(factor(
    hatch_training_meta$parasiteFactor,
    levels = unique(hatch_training_meta$parasiteFactor)
  )) - 1L,
  hatch_training_meta$Sample
)


## Check: Treatment vector integrity
if (
  !identical(names(treatment_no), hatch_training_meta$Sample) ||
    length(treatment_no) != nrow(hatch_training_meta) ||
    anyNA(treatment_no) ||
    !all(treatment_no %in% c(0L, 1L))
) {
  stop(
    "ERROR: Treatment vector has missing values, invalid classes, ",
    "or sample IDs that do not match the training metadata."
  )
}

message("PASS: Treatment vector integrity checks.")


## Save treatment vector
saveRDS(
  treatment_no,
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqInput_TreatmentVector_iter",
      iter,
      ".RDS"
    )
  )
)

message(
  "Completed PQLseq input preparation for iteration ",
  iter,
  "."
)