###          C1. PQLseq object preparation               ###

## Date: 

## Author:

## Description:
##    This script prepares training sample data for PQLseq    ##



####  Environment Preparation ####
library(methylKit)
library(dplyr)

## Pull iteration & available cores ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
numCore <- as.numeric(args[2])

## Set paths
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")
iterdir <- file.path(workdir, "tmp", "iterationData")

## Define functions to generate number of Cs and coverage data for PQLseq input
PrepNumCsFun <- function(list_of_dfs) {
  for (i in seq_along(list_of_dfs)) {
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% dplyr::select(-chr)

    col_names <- colnames(list_of_dfs[[i]])
    cols_to_delete <- grep("^coverage", col_names)
    if (length(cols_to_delete) > 0) {
      list_of_dfs[[i]] <- list_of_dfs[[i]][, -cols_to_delete, drop = FALSE]
    }

    list_of_dfs[[i]] <- list_of_dfs[[i]] %>%
      rename_with(~gsub("^numCs", "", .x))

    colnames(list_of_dfs[[i]]) <- paste0("ind", colnames(list_of_dfs[[i]]))
  }
  return(list_of_dfs)
}

PrepCoverageFun <- function(list_of_dfs) {
  for (i in seq_along(list_of_dfs)) {
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% dplyr::select(-chr)

    col_names <- colnames(list_of_dfs[[i]])
    cols_to_delete <- grep("^numCs", col_names)
    if (length(cols_to_delete) > 0) {
      list_of_dfs[[i]] <- list_of_dfs[[i]][, -cols_to_delete, drop = FALSE]
    }

    list_of_dfs[[i]] <- list_of_dfs[[i]] %>%
      rename_with(~gsub("coverage", "", .x))

    colnames(list_of_dfs[[i]]) <- paste0("ind", colnames(list_of_dfs[[i]]))
  }
  return(list_of_dfs)
}

## Load training sample uniteCov object for object generation
hatch_trainingSamples_uniteObj <- readRDS(
  file.path(tmpdir, paste0("hatch_trainingSamples_uniteObj_iter", iter, ".RDS"))
)

## Load metadata for training samples of the current iteration
hatch_trainingSamples_metadata <- read.csv(
  file.path(iterdir, paste0("hatch_trainingSamples_metadata_iter", iter, ".csv"))
)



####  Generate Kinship Matrix   ####
## Build the kinship matrix by maternal group
mother_ids <- as.character(hatch_trainingSamples_metadata$Mother)
same_mother <- outer(mother_ids, mother_ids, FUN = "==")
same_mother[is.na(same_mother)] <- FALSE
hatch_gen_mat <- 0.5 * same_mother
diag(hatch_gen_mat) <- 1
dimnames(hatch_gen_mat) <- list(
  hatch_trainingSamples_metadata$Sample,
  hatch_trainingSamples_metadata$Sample
)
hatch_gen_df <- as.data.frame(hatch_gen_mat)

## Save kinship matrix to file
saveRDS(
  hatch_gen_df,
  file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqInput_KinshipMatrix_iter", iter, ".RDS"))
)


####  Prepare Methylation Data Inputs for PQLseq  ####
## Convert uniteCov object to data.frame
hatch_trainingSamples_data <- as.data.frame(getData(hatch_trainingSamples_uniteObj))

## Rename rows to site names
rownames(hatch_trainingSamples_data) <- paste(
  hatch_trainingSamples_data$chr,
  hatch_trainingSamples_data$start,
  sep = "."
)

## Remove site, start, end and strand columns (keep chr column)
hatch_trainingSamples_data <- hatch_trainingSamples_data %>% dplyr::select(-start, -end, -strand)

## Delete all numTs columns: not needed by either
col_names <- colnames(hatch_trainingSamples_data)
cols_to_delete <- grep("^numTs", col_names)
if (length(cols_to_delete) > 0) {
  hatch_trainingSamples_data <- hatch_trainingSamples_data[, -cols_to_delete, drop = FALSE]
}

## Split data into chromosomes
hatch_trainingSamples_data$chr <- as.factor(hatch_trainingSamples_data$chr)
hatch_trainingSamples_data_list <- split(hatch_trainingSamples_data, hatch_trainingSamples_data$chr)

## Number of Cs inputs
hatch_trainingSamples_data_list_chr_df_numCs <- PrepNumCsFun(hatch_trainingSamples_data_list)

## save per-chromosome ##
lapply(
  names(hatch_trainingSamples_data_list_chr_df_numCs),
  function(chr_name) {
    saveRDS(
      hatch_trainingSamples_data_list_chr_df_numCs[[chr_name]],
      file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqInput_NumCs_iter", iter, "_", chr_name, ".RDS"))
    )
  }
)

## Coverage inputs
hatch_trainingSamples_data_list_chr_df_coverage <- PrepCoverageFun(hatch_trainingSamples_data_list)

## save per-chromosome ##
lapply(
  names(hatch_trainingSamples_data_list_chr_df_coverage),
  function(chr_name) {
    saveRDS(
      hatch_trainingSamples_data_list_chr_df_coverage[[chr_name]],
      file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqInput_Coverage_iter", iter, "_", chr_name, ".RDS"))
    )
  }
)

## Treatment vector
treatment.no <- setNames(
  as.integer(factor(
    hatch_trainingSamples_metadata$parasiteFactor,
    levels = unique(hatch_trainingSamples_metadata$parasiteFactor)
  )) - 1L,
  hatch_trainingSamples_metadata$Sample
)

## Save treatment to file ##
saveRDS(
  treatment.no,
  file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqInput_TreatmentVector_iter", iter, ".RDS"))
)

