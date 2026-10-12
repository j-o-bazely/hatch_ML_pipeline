###          C2. Run PQLseq               ###

## Date: 06.10.2026

## Author:

## Description:
##    This script runs PQLseq on the training sample data   ##


####  Environment Preparation ####

## Load custom PQLseq function ##
custom_pqlseq_script <- Sys.getenv(
  "CH2_CUSTOM_PQLSEQ_SCRIPT",
  unset = file.path(
    Sys.getenv("CH2_WORKDIR", unset = getwd()),
    "customPQLseq.R"
  )
)

if (file.exists(custom_pqlseq_script)) {
  source(custom_pqlseq_script)
} else {
  if (!requireNamespace("PQLseq", quietly = TRUE)) {
    stop(
      "PQLseq package is required, and no local customPQLseq.R file was found."
    )
  }
  custom_pqlseq <- get("pqlseq", asNamespace("PQLseq"))
}

## Set up environment for custom function ##
environment(custom_pqlseq) <- asNamespace("PQLseq")
assignInNamespace("pqlseq", custom_pqlseq, ns = "PQLseq")

## Pull iteration, available cores and chromosome ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
num_core <- as.numeric(args[2])
chr <- as.numeric(args[3]) - 1
# The array starts at 1, but the chromosome IDs range from 0 to 28.

## Set path to working directory ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")



####    Load Data for PQLseq    ####

chr_label <- paste0("SLK063_ragtag_chr", chr)

## Load methylation count inputs ##
num_cs_df <- readRDS(
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqInput_NumCs_iter",
      iter,
      "_",
      chr_label,
      ".RDS"
    )
  )
)

## Load coverage inputs ##
coverage_df <- readRDS(
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqInput_Coverage_iter",
      iter,
      "_",
      chr_label,
      ".RDS"
    )
  )
)

## Load genetic kinship matrix ##
relatedness_mat <- readRDS(
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqInput_KinshipMatrix_iter",
      iter,
      ".RDS"
    )
  )
)

## Load treatment vector ##
treatment_no <- readRDS(
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqInput_TreatmentVector_iter",
      iter,
      ".RDS"
    )
  )
)

#### Check" Samples alignment is consistent across objects ####
sample_cols <- paste0("ind", seq_along(treatment_no))

if (
!identical(colnames(num_cs_df), sample_cols) ||
!identical(colnames(coverage_df), sample_cols) ||
!identical(names(treatment_no), sample_mapping$Sample) ||
!identical(rownames(relatedness_mat), sample_mapping$Sample) ||
!identical(colnames(relatedness_mat), sample_mapping$Sample) ||
ncol(num_cs_df) != length(treatment_no) ||
ncol(coverage_df) != length(treatment_no)
) {
stop(
"ERROR: Sample alignment failed. Phenotype, kinship matrix, ",
"methylation counts and coverage must contain the same samples ",
"in the same order."
)
}
message(
"PASS: Phenotype, kinship matrix, methylation counts and coverage ",
"are aligned by sample."
)


####      Run PQLseq        ####
fit <- custom_pqlseq(
  RawCountDataSet = num_cs_df,
  Phenotypes = treatment_no,
  RelatednessMatrix = relatedness_mat,
  LibSize = coverage_df,
  fit.model = "BMM",
  verbose = FALSE,
  numCore = num_core
)


## Save fit object ##
saveRDS(
  fit,
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqOutput_Fit_iter",
      iter,
      "_",
      chr_label,
      ".RDS"
    )
  )
)
