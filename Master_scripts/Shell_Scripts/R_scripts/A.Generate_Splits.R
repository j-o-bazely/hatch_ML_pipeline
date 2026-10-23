###          A. Generate fixed training, calibration, and deployment cohorts ###

library(dplyr)

workdir <- Sys.getenv("CH2_WORKDIR")
if (!nzchar(workdir)) {
  stop("CH2_WORKDIR must be set to the pipeline data directory.")
}
workdir <- normalizePath(workdir, winslash = "/", mustWork = TRUE)
iterdir <- file.path(workdir, "tmp", "iterationData")
dir.create(iterdir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

meta <- read.csv(
  file.path(workdir, "input", "hatch_metadata.csv"),
  colClasses = c(Sample = "character", Mother = "character")
) |>
  dplyr::select(
    Sample,
    Mother,
    maternal_infection_status,
    Amount_cloacal_p,
    parasiteFactor
  )

required_columns <- c(
  "Sample",
  "Mother",
  "maternal_infection_status",
  "Amount_cloacal_p",
  "parasiteFactor"
)
if (!all(required_columns %in% names(meta))) {
  stop("Metadata is missing required sample, mother, or response columns.")
}
if (
  anyNA(meta$Sample) ||
  anyNA(meta$Mother) ||
  anyNA(meta$maternal_infection_status) ||
  anyNA(meta$Amount_cloacal_p) ||
  anyNA(meta$parasiteFactor)
) {
  stop("Sample metadata must not contain missing cohort or response values.")
}
if (anyDuplicated(meta$Sample)) {
  stop("Duplicate sample IDs detected in metadata.")
}
if (!setequal(unique(as.character(meta$parasiteFactor)), c("0", "1"))) {
  stop("parasiteFactor must contain exactly the binary labels 0 and 1.")
}

mother_meta <- meta |>
  dplyr::group_by(Mother) |>
  dplyr::summarise(
    maternal_infection_status = dplyr::first(maternal_infection_status),
    Amount_cloacal_p = dplyr::first(Amount_cloacal_p),
    .groups = "drop"
  )
mother_consistency <- meta |>
  dplyr::group_by(Mother) |>
  dplyr::summarise(
    statuses = dplyr::n_distinct(maternal_infection_status),
    burdens = dplyr::n_distinct(Amount_cloacal_p),
    .groups = "drop"
  )
if (
  any(mother_consistency$statuses != 1) ||
  any(mother_consistency$burdens != 1)
) {
  stop("Each mother must have one maternal infection status and burden stratum.")
}

samples_per_mother <- table(meta$Mother)
if (any(samples_per_mother != 4)) {
  stop("Expected exactly four hatchling samples per mother.")
}
if (nrow(mother_meta) != 29) {
  stop(
    "The fixed cohort design requires exactly 29 maternal groups: ",
    "19 training, 6 threshold calibration, and 4 deployment."
  )
}

select_mothers <- function(candidates, n, description) {
  if (nrow(candidates) < n) {
    stop(
      "Need ", n, " eligible maternal groups for ", description,
      ", but found ", nrow(candidates), "."
    )
  }
  sample(as.character(candidates$Mother), size = n, replace = FALSE)
}

deploy_mothers <- c(
  select_mothers(
    dplyr::filter(
      mother_meta,
      maternal_infection_status == "yes",
      Amount_cloacal_p == "1 to 5"
    ),
    n = 1,
    description = "the low-burden infected deployment cohort"
  ),
  select_mothers(
    dplyr::filter(
      mother_meta,
      maternal_infection_status == "yes",
      Amount_cloacal_p == ">26"
    ),
    n = 1,
    description = "the high-burden infected deployment cohort"
  ),
  select_mothers(
    dplyr::filter(mother_meta, maternal_infection_status == "no"),
    n = 2,
    description = "the two uninfected deployment mothers"
  )
)
if (anyDuplicated(deploy_mothers)) {
  stop("Deployment strata selected the same mother more than once.")
}

remaining_mothers <- dplyr::filter(
  mother_meta,
  !(Mother %in% deploy_mothers)
)

sample_stratified_mothers <- function(group_meta, n) {
  group_counts <- table(group_meta$maternal_infection_status)
  exact_allocation <- n * group_counts / sum(group_counts)
  allocation <- floor(exact_allocation)
  remainder <- n - sum(allocation)
  tie_break <- runif(length(group_counts))
  allocation_order <- order(
    -(exact_allocation - allocation),
    tie_break
  )

  while (remainder > 0) {
    eligible <- allocation_order[
      allocation[allocation_order] < group_counts[allocation_order]
    ]
    if (length(eligible) == 0) {
      stop("Could not allocate the requested number of maternal groups.")
    }
    allocation[eligible[[1]]] <- allocation[eligible[[1]]] + 1
    remainder <- remainder - 1
  }

  selected <- unlist(
    lapply(names(group_counts), function(status) {
      status_mothers <- as.character(
        group_meta$Mother[
          group_meta$maternal_infection_status == status
        ]
      )
      sample(
        status_mothers,
        size = allocation[[status]],
        replace = FALSE
      )
    }),
    use.names = FALSE
  )
  selected
}

remaining_infected <- dplyr::filter(
  remaining_mothers,
  maternal_infection_status == "yes"
)
remaining_uninfected <- dplyr::filter(
  remaining_mothers,
  maternal_infection_status == "no"
)
if (nrow(remaining_infected) < 3 || nrow(remaining_uninfected) < 3) {
  stop(
    "At least three infected and three uninfected mothers must remain after ",
    "deployment selection to create a balanced calibration cohort."
  )
}
calibration_mothers <- c(
  sample(as.character(remaining_infected$Mother), size = 3),
  sample(as.character(remaining_uninfected$Mother), size = 3)
)
training_mothers <- setdiff(
  as.character(remaining_mothers$Mother),
  calibration_mothers
)

training_cohort <- meta[meta$Mother %in% training_mothers, , drop = FALSE]
calibration_cohort <- meta[meta$Mother %in% calibration_mothers, , drop = FALSE]
deployment_cohort <- meta[meta$Mother %in% deploy_mothers, , drop = FALSE]
deployment_mother_meta <- dplyr::filter(
  mother_meta,
  Mother %in% deploy_mothers
)
if (
  sum(deployment_mother_meta$maternal_infection_status == "yes") != 2 ||
  sum(deployment_mother_meta$maternal_infection_status == "no") != 2
) {
  stop("Deployment cohort must contain two infected and two uninfected mothers.")
}

cohort_mothers <- list(
  training = unique(as.character(training_cohort$Mother)),
  calibration = unique(as.character(calibration_cohort$Mother)),
  deployment = unique(as.character(deployment_cohort$Mother))
)
if (
  length(cohort_mothers$training) != 19 ||
  length(cohort_mothers$calibration) != 6 ||
  length(cohort_mothers$deployment) != 4 ||
  length(intersect(cohort_mothers$training, cohort_mothers$calibration)) > 0 ||
  length(intersect(cohort_mothers$training, cohort_mothers$deployment)) > 0 ||
  length(intersect(cohort_mothers$calibration, cohort_mothers$deployment)) > 0
) {
  stop("The fixed training, calibration, and deployment cohorts are invalid.")
}

calibration_mother_meta <- dplyr::filter(
  mother_meta,
  Mother %in% calibration_mothers
)
if (
  sum(calibration_mother_meta$maternal_infection_status == "yes") != 3 ||
  sum(calibration_mother_meta$maternal_infection_status == "no") != 3
) {
  stop("Threshold calibration cohort must contain three mothers per class.")
}

for (cohort in list(training_cohort, calibration_cohort, deployment_cohort)) {
  if (!setequal(unique(as.character(cohort$parasiteFactor)), c("0", "1"))) {
    stop("Each fixed cohort must include samples from both response classes.")
  }
}

write.csv(
  training_cohort,
  file.path(iterdir, "hatch_modelTrainingCohort_metadata.csv"),
  row.names = FALSE
)
write.csv(
  calibration_cohort,
  file.path(iterdir, "hatch_thresholdCalibration_metadata.csv"),
  row.names = FALSE
)
write.csv(
  deployment_cohort,
  file.path(iterdir, "hatch_deploymentSamples_metadata.csv"),
  row.names = FALSE
)

message(
  "Fixed cohorts: ",
  length(cohort_mothers$training), " training mothers (",
  nrow(training_cohort), " hatchlings), ",
  length(cohort_mothers$calibration), " calibration mothers (",
  nrow(calibration_cohort), " hatchlings), and ",
  length(cohort_mothers$deployment), " deployment mothers (",
  nrow(deployment_cohort), " hatchlings)."
)

iteration_text <- Sys.getenv(
  "CH2_ITERATIONS",
  unset = paste(1:1000, collapse = ",")
)
iterations <- as.integer(strsplit(iteration_text, ",", fixed = TRUE)[[1]])
if (
  length(iterations) == 0 ||
  anyNA(iterations) ||
  any(iterations < 1) ||
  anyDuplicated(iterations)
) {
  stop("CH2_ITERATIONS must contain unique positive iteration numbers.")
}

training_mother_meta <- dplyr::filter(
  mother_meta,
  Mother %in% training_mothers
)
for (iteration in iterations) {
  ## Resample only within the fixed training cohort for DMS stability selection.
  iteration_mothers <- sample_stratified_mothers(
    training_mother_meta,
    n = round(0.8 * nrow(training_mother_meta))
  )
  train <- training_cohort[
    training_cohort$Mother %in% iteration_mothers,
    ,
    drop = FALSE
  ]
  if (!setequal(unique(as.character(train$parasiteFactor)), c("0", "1"))) {
    stop(
      "DMS resample ", iteration,
      " does not include both response classes; revise the cohort allocation."
    )
  }

  write.csv(
    train,
    file.path(
      iterdir,
      paste0("hatch_trainingSamples_metadata_iter", iteration, ".csv")
    ),
    row.names = FALSE
  )
  write.csv(
    calibration_cohort,
    file.path(
      iterdir,
      paste0("hatch_validationSamples_metadata_iter", iteration, ".csv")
    ),
    row.names = FALSE
  )
}
