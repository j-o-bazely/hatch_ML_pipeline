### Prepare full-cohort feature matrices for the final ML model ###

library(methylKit)
library(GenomicRanges)

workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
iterdir <- file.path(workdir, "tmp", "iterationData")
ml_dir <- file.path(workdir, "output", "MLdata")
dir.create(ml_dir, recursive = TRUE, showWarnings = FALSE)

read_cohort_metadata <- function(filename) {
  metadata <- read.csv(
    file.path(iterdir, filename),
    colClasses = c(Sample = "character", Mother = "character")
  )
  required <- c("Sample", "Mother", "parasiteFactor")
  if (!all(required %in% names(metadata))) {
    stop("Cohort metadata is missing required columns: ", filename)
  }
  if (
    anyNA(metadata$Sample) ||
    anyNA(metadata$Mother) ||
    anyNA(metadata$parasiteFactor) ||
    anyDuplicated(metadata$Sample)
  ) {
    stop("Cohort metadata contains invalid sample IDs or labels: ", filename)
  }
  metadata
}

training_meta <- read_cohort_metadata("hatch_modelTrainingCohort_metadata.csv")
calibration_meta <- read_cohort_metadata(
  "hatch_thresholdCalibration_metadata.csv"
)
deployment_meta <- read_cohort_metadata("hatch_deploymentSamples_metadata.csv")

cohort_mothers <- list(
  training = unique(training_meta$Mother),
  calibration = unique(calibration_meta$Mother),
  deployment = unique(deployment_meta$Mother)
)
if (
  length(intersect(cohort_mothers$training, cohort_mothers$calibration)) > 0 ||
  length(intersect(cohort_mothers$training, cohort_mothers$deployment)) > 0 ||
  length(intersect(cohort_mothers$calibration, cohort_mothers$deployment)) > 0
) {
  stop("A maternal group occurs in more than one final modeling cohort.")
}

stable_dms <- read.csv(
  file.path(workdir, "output", "MLresults", "stable_DMS_ge50pct.csv"),
  colClasses = c(chr = "character")
)

required_dms_columns <- c("chr", "start", "end")

if (!all(required_dms_columns %in% names(stable_dms)) || nrow(stable_dms) == 0) {
  stop("Stable DMS file is empty or missing genomic coordinate columns.")
}
if (
  anyNA(stable_dms[, required_dms_columns]) ||
  anyDuplicated(paste(stable_dms$chr, stable_dms$start, stable_dms$end))
) {
  stop("Stable DMS coordinates must be complete and unique.")
}

stable_dms_ranges <- makeGRangesFromDataFrame(
  stable_dms,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = FALSE
)

expected_feature_ids <- paste(stable_dms$chr, stable_dms$start, sep = ".")

full_methylation <- readRDS(
  file.path(
    workdir,
    "input",
    "hatchUnite_100minSmallestGroup_minCov5_filtered_normal.RDS"
  )
)

make_feature_matrix <- function(metadata, cohort_name) {
  missing_samples <- setdiff(metadata$Sample, full_methylation@sample.ids)
  if (length(missing_samples) > 0) {
    stop(
      cohort_name,
      " samples are missing from the methylation object: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  ## Reorganise object to keep only cohort-specific samples ##
  cohort_object <- reorganize(
    full_methylation,
    sample.ids = metadata$Sample,
    treatment = as.numeric(as.character(metadata$parasiteFactor))
  )

  ## Subset for consensus DMS ##
  feature_object <- selectByOverlap(cohort_object, stable_dms_ranges)

  feature_matrix <- as.data.frame(
    percMethylation(feature_object, rowids = TRUE)
  )

  if (
    nrow(feature_matrix) != length(expected_feature_ids) ||
    !setequal(
      sub("\\.([0-9]+)$", "", rownames(feature_matrix)),
      expected_feature_ids
    )
  ) {
    stop(
      cohort_name,
      " feature matrix does not contain exactly the consensus DMS."
    )
  }

  feature_matrix <- feature_matrix[expected_feature_ids, , drop = FALSE]
  if (!identical(colnames(feature_matrix), metadata$Sample)) {
    stop(cohort_name, " feature matrix sample order does not match metadata.")
  }
  feature_matrix
}

cohorts <- list(
  modelTraining = training_meta,
  thresholdCalibration = calibration_meta,
  deployment = deployment_meta
)

for (cohort_name in names(cohorts)) {
  metadata <- cohorts[[cohort_name]]
  feature_matrix <- make_feature_matrix(metadata, cohort_name)
  write.csv(
    feature_matrix,
    file.path(
      ml_dir,
      paste0("hatch_", cohort_name, "_DMS_percMeth_final.csv")
    ),
    row.names = TRUE
  )
  write.csv(
    metadata,
    file.path(ml_dir, paste0("hatch_", cohort_name, "_metadata_final.csv")),
    row.names = FALSE
  )
}

message(
  "Prepared final model matrices for ",
  nrow(training_meta), " training, ",
  nrow(calibration_meta), " threshold-calibration, and ",
  nrow(deployment_meta), " deployment hatchlings across ",
  nrow(stable_dms), " consensus DMS."
)
