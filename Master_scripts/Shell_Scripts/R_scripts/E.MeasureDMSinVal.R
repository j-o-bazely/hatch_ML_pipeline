###          E. Measure DMS in Validation Samples                ##

## Date: 06.10.2026

## Author: James B

## Description:
##  This code analyses the clustering of validation samples at DMS identified in training samples ##
##  before then generating the required objects for downstream machine learning testing           ##


####  Environment Preparation ####
## Load libraries ##
library(methylKit)
library(GenomicRanges)
library(dplyr)
library(vegan)

## Pull iteration & available cores ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
numCore <- as.numeric(args[2])

## Set paths ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")
iterdir <- file.path(workdir, "tmp", "iterationData")
outdir <- file.path(workdir, "output")
MLdir <- file.path(outdir, "MLdata")

## Set seed for adonis calculations ##
set.seed(123)

## Load objects ##
hatch_trainingSamples_uniteObj <- readRDS(
  file.path(tmpdir, paste0("hatch_trainingSamples_uniteObj_iter", iter, ".RDS"))
)

iterationDMS <- readRDS(file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqOutput_DMS_q0.05_iter", iter, ".RDS")))

## Load metadata ##
hatch_training_meta <- read.csv(file.path(iterdir, paste0("hatch_trainingSamples_metadata_iter", iter, ".csv"))) %>%
  dplyr::select(
    Sample,
    parasiteFactor
  )

hatch_validation_meta <- read.csv(file.path(iterdir, paste0("hatch_validationSamples_metadata_iter", iter, ".csv"))) %>%
  dplyr::select(
    Sample,
    parasiteFactor
  )



####    Matrix of Iteration DMS in Training Samples     ####

## methylKit DMS results use chr/start/end genomic coordinates.
iterationDMS_df <- as.data.frame(getData(iterationDMS))
required_columns <- c("chr", "start", "end")
if (!all(required_columns %in% names(iterationDMS_df))) {
  stop("DMS object must contain chr, start, and end columns.")
}

GRanges_iteration_DMS <- makeGRangesFromDataFrame(
  iterationDMS_df,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = FALSE
)

## Selects records that match list of DMS locations ##
hatch_trainingSamples_DMS_uniteObj <- selectByOverlap(hatch_trainingSamples_uniteObj, GRanges_iteration_DMS)

## Generate percent methylation matrix ##
hatch_trainingSamples_DMS_percMeth <- percMethylation(hatch_trainingSamples_DMS_uniteObj, rowids = T)
hatch_trainingSamples_DMS_percMeth_df <- as.data.frame(hatch_trainingSamples_DMS_percMeth)



####    Matrix of Iteration DMS in Validation Samples     ####
## Load object ##
hatch_validationSamples_uniteObj <- readRDS(
  file.path(tmpdir, paste0("hatch_validationSamples_uniteObj_iter", iter, ".RDS"))
)

### Use methylKit select by overlap ###
# Selects records that match list of DMS locations
hatch_validationSamples_DMS_uniteObj <- selectByOverlap(hatch_validationSamples_uniteObj, GRanges_iteration_DMS)

## Create percent methylation matrix ##
hatch_validationSamples_DMS_percMeth <- percMethylation(hatch_validationSamples_DMS_uniteObj, rowids = T)
hatch_validationSamples_DMS_percMeth_df <- as.data.frame(hatch_validationSamples_DMS_percMeth)


####      PERMANOVA of DMS        ####
###   Training Samples    ###
## Generate distance matrix for training samples ##
hatch_trainingSamples_DMS_percMeth_complete <- 
  hatch_trainingSamples_DMS_percMeth[rowSums(is.na(hatch_trainingSamples_DMS_percMeth)) == 0,]

hatch_trainingSamples_DMS_dist <- dist(
  t(hatch_trainingSamples_DMS_percMeth_complete),
  method = "euclidean"
)

## Check that metadata order matches data ##
identical(rownames(as.matrix(hatch_trainingSamples_DMS_dist)), hatch_training_meta$Sample)

## Run PERMANOVA test on this distance matrix ##
hatch_trainingSamples_DMS_adonis2 <- vegan::adonis2(
  hatch_trainingSamples_DMS_dist ~ parasiteFactor,
  data = hatch_training_meta,
  permutations = 999
)

###   Validation Samples    ###
## Generate distance matrix for validation samples ##
hatch_validationSamples_DMS_percMeth_complete <- 
  hatch_validationSamples_DMS_percMeth[rowSums(is.na(hatch_validationSamples_DMS_percMeth)) == 0,]

hatch_validationSamples_DMS_dist <- dist(t(hatch_validationSamples_DMS_percMeth_complete), method = 'euclidean') # Check best method for genetic data #

## Check that metadata order matches data ##
identical(rownames(as.matrix(hatch_validationSamples_DMS_dist)), hatch_validation_meta$Sample)


## Run PERMANOVA test on this distance matrix ##
hatch_validationSamples_DMS_adonis2 <- vegan::adonis2(
  hatch_validationSamples_DMS_dist ~ parasiteFactor,
  data = hatch_validation_meta,
  permutations = 999
)

## Save adonis results to file ##
permanovaResults <- data.frame(
  iteration = iter,
  trainingPval = hatch_trainingSamples_DMS_adonis2$`Pr(>F)`[1],
  trainingR2 = hatch_trainingSamples_DMS_adonis2$R2[1],
  validationPval = hatch_validationSamples_DMS_adonis2$`Pr(>F)`[1],
  validationR2 = hatch_validationSamples_DMS_adonis2$R2[1]
)

## Save table to 
write.table(permanovaResults, paste0(outdir, '/permanovaResults/iter', iter, '_permanovaResults.txt'))



#### Prepare Data for ML ####
## Check that orders match beetween matix and metadata ##
if (!identical(hatch_training_meta$Sample, colnames(hatch_trainingSamples_DMS_percMeth_df))) {
  stop("Training metadata order does not match the training methylation matrix columns.")
}
if (!identical(hatch_validation_meta$Sample, colnames(hatch_validationSamples_DMS_percMeth_df))) {
  stop("Validation metadata order does not match the validation methylation matrix columns.")
}

## Save objects for use in ML ##
saveRDS(hatch_validationSamples_DMS_percMeth_df, file.path(MLdir, paste0("hatch_validationSamples_DMS_percMeth_df_iter", iter, ".RDS")))
saveRDS(hatch_validationSamples_DMS_uniteObj, file = file.path(MLdir, paste0("hatch_validationSamples_DMS_uniteObj_iter", iter, ".RDS")))

saveRDS(hatch_trainingSamples_DMS_percMeth_df, file.path(MLdir, paste0("hatch_trainingSamples_DMS_percMeth_df_iter", iter, ".RDS")))
saveRDS(hatch_trainingSamples_DMS_uniteObj, file = file.path(MLdir, paste0("hatch_trainingSamples_DMS_uniteObj_iter", iter, ".RDS")))

write.csv(
  hatch_training_meta,
  file.path(MLdir, paste0("hatch_trainingSamples_metadata_iter", iter, ".csv")),
  row.names = FALSE
)

write.csv(
  hatch_validation_meta,
  file.path(MLdir, paste0("hatch_validationSamples_metadata_iter", iter, ".csv")),
  row.names = FALSE
)

