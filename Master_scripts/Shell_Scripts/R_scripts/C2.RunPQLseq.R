###          C2. Run PQLseq               ###

## Date: 

## Author:

## Description:
##    This script runs PQLseq on the training sample data   ##


####  Environment Preparation ####

## Load custom PQLseq function if a project-local override is available.
## Otherwise, fall back to the installed PQLseq package.
custom_pqlseq_script <- Sys.getenv(
  "CH2_CUSTOM_PQLSEQ_SCRIPT",
  unset = file.path(Sys.getenv("CH2_WORKDIR", unset = getwd()), "customPQLseq.R")
)

if (file.exists(custom_pqlseq_script)) {
  source(custom_pqlseq_script)
} else {
  if (!requireNamespace("PQLseq", quietly = TRUE)) {
    stop("PQLseq package is required, and no local customPQLseq.R file was found.")
  }
  custom_pqlseq <- get("pqlseq", asNamespace("PQLseq"))
}

## Set up environment for custom function ##
environment(custom_pqlseq) <- asNamespace('PQLseq')
assignInNamespace("pqlseq", custom_pqlseq, ns = "PQLseq")

## Pull iteration, available cores and chromosome ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
numCore <- as.numeric(args[2])
chr <- as.numeric(args[3])-1 # *-1 IS *IMPORTANT. We want to run chr0 to 28 but array has to start with 1, so we run from 1-29

## Set path to working directory ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")



####    Load Data for PQLseq    ####    

chr_label <- paste0("SLK063_ragtag_chr", chr)

## Load methylation count inputs ##
numCs_df <- readRDS(file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqInput_NumCs_iter", iter, "_", chr_label, ".RDS")))
  
## Load coverage inputs ##
coverage_df <- readRDS(file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqInput_Coverage_iter", iter, "_", chr_label, ".RDS")))

## Load genetic kinship matrix ##
relatedness_mat <- readRDS(file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqInput_KinshipMatrix_iter", iter, ".RDS")))

## Load treatment vector ##
treatment.no <- readRDS(file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqInput_TreatmentVector_iter", iter, ".RDS")))


####      Run PQLseq        ####
fit = custom_pqlseq(RawCountDataSet = numCs_df , Phenotypes = treatment.no ,
                    RelatednessMatrix = relatedness_mat , LibSize = coverage_df ,
                    fit.model = "BMM", verbose=FALSE, numCore = numCore)


## Save fit object ##
saveRDS(fit, file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqOutput_Fit_iter", iter, "_", chr_label, ".RDS")))
