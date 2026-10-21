###          F. Measure DMS in Validation Samples                ##

## Date: 06.10.2026

## Author: James B

## Description:
## This code analyzes validation sample clustering at training-identified DMS,
## then creates objects for downstream machine-learning tests.


####  Environment Preparation ####
## Load libraries ##
library(methylKit)
library(GenomicRanges)
library(dplyr)
library(vegan)

## Pull iteration & available cores ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
num_core <- as.numeric(args[2])

## Set paths ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")
iterdir <- file.path(workdir, "tmp", "iterationData")
outdir <- file.path(workdir, "output")
ml_dir <- file.path(outdir, "MLdata")

## Set seed for adonis calculations ##
set.seed(123)

## Load objects ##
hatch_training_unite_obj <- readRDS(
  file.path(tmpdir, paste0("hatch_trainingSamples_uniteObj_iter", iter, ".RDS"))
)

iteration_dms <- readRDS(
  file.path(
    tmpdir,
    paste0("hatch_trainingSamples_PQLseqOutput_DMS_q0.05_iter", iter, ".RDS")
  )
)

## Load metadata ##
hatch_training_meta <- read.csv(
  file.path(
    iterdir,
    paste0("hatch_trainingSamples_metadata_iter", iter, ".csv")
  )
) |>
  dplyr::select(
    Sample,
    parasiteFactor
  )

hatch_validation_meta <- read.csv(
  file.path(
    iterdir,
    paste0("hatch_validationSamples_metadata_iter", iter, ".csv")
  )
) |>
  dplyr::select(
    Sample,
    parasiteFactor
  )




####    Matrix of Iteration DMS in Training Samples     ####

## methylKit DMS results use chr/start/end genomic coordinates.
iteration_dms_df <- as.data.frame(getData(iteration_dms))
required_columns <- c("chr", "start", "end")
if (!all(required_columns %in% names(iteration_dms_df))) {
  stop("DMS object must contain chr, start, and end columns.")
}

granges_iteration_dms <- makeGRangesFromDataFrame(
  iteration_dms_df,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = FALSE
)

## Selects records that match list of DMS locations ##
hatch_training_dms_unite_obj <- selectByOverlap(
  hatch_training_unite_obj,
  granges_iteration_dms
)

## Generate percent methylation matrix ##
hatch_training_dms_perc_meth <- percMethylation(
  hatch_training_dms_unite_obj,
  rowids = TRUE
)
training_dms_matrix <- as.data.frame(
  hatch_training_dms_perc_meth
)



####    Matrix of Iteration DMS in Validation Samples     ####
## Load object ##
hatch_validation_unite_obj <- readRDS(
  file.path(
    tmpdir,
    paste0("hatch_validationSamples_uniteObj_iter", iter, ".RDS")
  )
)

### Use methylKit select by overlap ###
# Selects records that match list of DMS locations
hatch_validation_dms_unite_obj <- selectByOverlap(
  hatch_validation_unite_obj,
  granges_iteration_dms
)

## Create percent methylation matrix ##
hatch_validation_dms_perc_meth <- percMethylation(
  hatch_validation_dms_unite_obj,
  rowids = TRUE
)
validation_dms_matrix <- as.data.frame(
  hatch_validation_dms_perc_meth
)

## Build ML features from the cross-iteration stable DMS set ##
stable_dms <- read.csv(
  file.path(workdir, "output", "MLresults", "stable_DMS_gt95pct.csv"),
  colClasses = c(chr = "character")
)
if (nrow(stable_dms) == 0) {
  stop("The stable DMS feature set is empty.")
}
required_columns <- c("chr", "start", "end")
if (!all(required_columns %in% names(stable_dms))) {
  stop("Stable DMS file is missing genomic coordinate columns.")
}

stable_dms_ranges <- makeGRangesFromDataFrame(
  stable_dms,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = FALSE
)
training_stable_dms_unite_obj <- selectByOverlap(
  hatch_training_unite_obj,
  stable_dms_ranges
)
validation_stable_dms_unite_obj <- selectByOverlap(
  hatch_validation_unite_obj,
  stable_dms_ranges
)
training_stable_dms_matrix <- as.data.frame(
  percMethylation(training_stable_dms_unite_obj, rowids = TRUE)
)
validation_stable_dms_matrix <- as.data.frame(
  percMethylation(validation_stable_dms_unite_obj, rowids = TRUE)
)

if (
  nrow(training_stable_dms_matrix) != nrow(stable_dms) ||
  nrow(validation_stable_dms_matrix) != nrow(stable_dms)
) {
  stop("Not all stable DMS features are present in both sample sets.")
}
if (
  !identical(
    rownames(training_stable_dms_matrix),
    rownames(validation_stable_dms_matrix)
  )
) {
  stop("Training and validation stable DMS features do not match.")
}
if (
  !identical(
    colnames(training_stable_dms_matrix),
    hatch_training_meta$Sample
  ) ||
  !identical(
    colnames(validation_stable_dms_matrix),
    hatch_validation_meta$Sample
  )
) {
  stop("Stable DMS matrix columns do not match their sample metadata.")
}



####      PERMANOVA of DMS        ####
###   Training Samples    ###
## Generate distance matrix for training samples ##
training_dms_complete <-
  hatch_training_dms_perc_meth[
    rowSums(is.na(hatch_training_dms_perc_meth)) == 0,
  ]

hatch_training_dms_dist <- dist(
  t(training_dms_complete),
  method = "euclidean"
)

## Check that metadata order matches data ##
identical(
  rownames(as.matrix(hatch_training_dms_dist)),
  hatch_training_meta$Sample
)

## Run PERMANOVA test on this distance matrix ##
hatch_training_dms_adonis <- vegan::adonis2(
  hatch_training_dms_dist ~ parasiteFactor,
  data = hatch_training_meta,
  permutations = 999
)

###   Validation Samples    ###
## Generate distance matrix for validation samples ##
validation_dms_complete <-
  hatch_validation_dms_perc_meth[
    rowSums(is.na(hatch_validation_dms_perc_meth)) == 0,
  ]

hatch_validation_dms_dist <- dist(
  t(validation_dms_complete),
  method = "euclidean"
)

## Check that metadata order matches data ##
identical(
  rownames(as.matrix(hatch_validation_dms_dist)),
  hatch_validation_meta$Sample
)


## Run PERMANOVA test on this distance matrix ##
hatch_validation_dms_adonis <- vegan::adonis2(
  hatch_validation_dms_dist ~ parasiteFactor,
  data = hatch_validation_meta,
  permutations = 999
)

## Save adonis results to file ##
permanova_results <- data.frame(
  iteration = iter,
  trainingPval = hatch_training_dms_adonis$`Pr(>F)`[1],
  trainingR2 = hatch_training_dms_adonis$R2[1],
  validationPval = hatch_validation_dms_adonis$`Pr(>F)`[1],
  validationR2 = hatch_validation_dms_adonis$R2[1]
)

## Save table to file ##
write.table(
  permanova_results,
  paste0(outdir, "/permanovaResults/iter", iter, "_permanovaResults.txt")
)



#### Prepare Data for ML ####
## Check that orders match beetween matix and metadata ##
if (!identical(
  hatch_training_meta$Sample,
  colnames(training_dms_matrix)
)) {
  stop(
    paste0(
      "Training metadata order does not match the ",
      "training methylation matrix columns."
    )
  )
}
if (!identical(
  hatch_validation_meta$Sample,
  colnames(validation_dms_matrix)
)) {
  stop(
    paste0(
      "Validation metadata order does not match the ",
      "validation methylation matrix columns."
    )
  )
}

## Save objects for use in ML ##
saveRDS(
  validation_stable_dms_matrix,
  file.path(
    ml_dir,
    paste0("hatch_validationSamples_DMS_percMeth_df_iter", iter, ".RDS")
  )
)
saveRDS(
  validation_stable_dms_unite_obj,
  file = file.path(
    ml_dir,
    paste0("hatch_validationSamples_DMS_uniteObj_iter", iter, ".RDS")
  )
)

saveRDS(
  training_stable_dms_matrix,
  file.path(
    ml_dir,
    paste0("hatch_trainingSamples_DMS_percMeth_df_iter", iter, ".RDS")
  )
)
saveRDS(
  training_stable_dms_unite_obj,
  file = file.path(
    ml_dir,
    paste0("hatch_trainingSamples_DMS_uniteObj_iter", iter, ".RDS")
  )
)

write.csv(
  hatch_training_meta,
  file.path(
    ml_dir,
    paste0("hatch_trainingSamples_metadata_iter", iter, ".csv")
  ),
  row.names = FALSE
)

write.csv(
  hatch_validation_meta,
  file.path(
    ml_dir,
    paste0("hatch_validationSamples_metadata_iter", iter, ".csv")
  ),
  row.names = FALSE
)
