###          D. Combine PQLseq and methylKit results; select DMS          ###

## Date: 06.10.2026

## Author: James B

## Description:
## This script combines significance values from PQLseq with methylation
## difference values from methylKit, then selects significant DMS.


####  Environment Preparation ####
## Load required libraries ##
library(methylKit)

## Pull iteration & available cores ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
num_core <- as.numeric(args[2])

## Set directories ##
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")

## Load slim output ##
slim_output <- readRDS(
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqFit_iter",
      iter,
      "_wholeGenome_corrected.RDS"
    )
  )
)

rownames(slim_output) <- gsub("chr[0-9]+\\.S", "S", rownames(slim_output))

head(rownames(slim_output))

## Load methylKit methdiff object ##
methdiff <- readRDS(
  file.path(tmpdir, paste0("hatch_trainingSamples_diffMeth_iter", iter, ".RDS"))
)




#### Prepare Data for Merging ####
methdiff_df <- as.data.frame(getData(methdiff))
methdiff_df$pvalue <- NULL
methdiff_df$qvalue <- NULL
methdiff_df$siteloc <- paste(methdiff_df$chr, methdiff_df$start, sep = ".")
rownames(methdiff_df) <- methdiff_df$siteloc
slim_output$siteloc <- rownames(slim_output)

## Check genomic identifiers ##

if (
  anyNA(methdiff_df$siteloc) ||
  anyNA(slim_output$siteloc) ||
  any(methdiff_df$siteloc == "") ||
  any(slim_output$siteloc == "")
) {
  stop("ERROR: Missing genomic site identifiers.")
}

if (
  anyDuplicated(methdiff_df$siteloc) ||
  anyDuplicated(slim_output$siteloc)
) {
  stop("ERROR: Duplicate genomic locations prevent matching results.")
}

reorder_idx <- match(methdiff_df$siteloc, slim_output$siteloc)
converged_match <- !is.na(reorder_idx)
if (any(!converged_match)) {
  message(
    sum(!converged_match),
    " methylKit sites were omitted because their PQLseq fits did not converge."
  )
}
methdiff_df <- methdiff_df[converged_match, , drop = FALSE]
slim_output_reorder <- slim_output[reorder_idx[converged_match], , drop = FALSE]

merged_df <- tibble::add_column(
  methdiff_df,
  pvalue = slim_output_reorder$pvalue,
  .before = "meth.diff"
)
merged_df <- tibble::add_column(
  merged_df,
  qvalue = slim_output_reorder$qvalue,
  .before = "meth.diff"
)

####  Generate Updated methylDiff Object  ####
merged_df$converged <- slim_output_reorder$converged
merged_df <- merged_df[merged_df$converged == "TRUE", , drop = FALSE]
merged_df$siteloc <- NULL
merged_df$converged <- NULL

new_methdiff <- new(
  Class = "methylDiff",
  merged_df,
  sample.ids = methdiff@sample.ids,
  assembly = methdiff@assembly,
  context = methdiff@context,
  treatment = methdiff@treatment,
  destranded = methdiff@destranded,
  resolution = methdiff@resolution,
  names = methdiff@names,
  row.names = rownames(merged_df)
)

saveRDS(
  new_methdiff,
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqOutput_DiffMeth_Corrected_iter",
      iter,
      ".RDS"
    )
  )
)

pql_dms_10pc_q0_05_no0q <- getMethylDiff(
  new_methdiff,
  difference = 10,
  qvalue = 0.05
)

pql_dms_10pc_q0_01_no0q <- getMethylDiff(
  new_methdiff,
  difference = 10,
  qvalue = 0.01
)

saveRDS(
  pql_dms_10pc_q0_05_no0q,
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqOutput_DMS_q0.05_iter",
      iter,
      ".RDS"
    )
  )
)

saveRDS(
  pql_dms_10pc_q0_01_no0q,
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqOutput_DMS_q0.01_iter",
      iter,
      ".RDS"
    )
  )
)
