#### E. Select DMS found in more than 95% of requested iterations ####

library(methylKit)

workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)

iteration_text <- Sys.getenv("CH2_ITERATIONS")
if (!nzchar(iteration_text)) {
  stop("CH2_ITERATIONS is required.")
}
iterations <- unique(as.integer(strsplit(iteration_text, ",", fixed = TRUE)[[1]]))
if (anyNA(iterations) || any(iterations < 1)) {
  stop("CH2_ITERATIONS must contain positive iteration numbers.")
}

site_tables <- vector("list", length(iterations))
for (index in seq_along(iterations)) {
  iteration <- iterations[[index]]
  marker_file <- file.path(
    workdir,
    "tmp",
    paste0("dms_ready_iter", iteration, ".success")
  )
  dms_file <- file.path(
    workdir,
    "tmp",
    paste0("hatch_trainingSamples_PQLseqOutput_DMS_q0.05_iter", iteration, ".RDS")
  )

  if (!file.exists(marker_file)) {
    message("Iteration ", iteration, " did not produce DMS results; counting as absent.")
    next
  }
  if (file.access(dms_file, mode = 4) != 0) {
    stop("DMS results are marked complete but cannot be read: ", dms_file)
  }

  dms <- readRDS(dms_file)
  dms_data <- as.data.frame(getData(dms))
  required_columns <- c("chr", "start", "end")
  if (!all(required_columns %in% names(dms_data))) {
    stop("DMS results are missing genomic coordinate columns: ", dms_file)
  }
  sites <- unique(dms_data[, required_columns, drop = FALSE])
  sites$chr <- as.character(sites$chr)
  sites$start <- as.integer(sites$start)
  sites$end <- as.integer(sites$end)
  if (anyNA(sites)) {
    stop("DMS results contain missing or invalid coordinates: ", dms_file)
  }
  site_tables[[index]] <- sites
}

site_tables <- Filter(Negate(is.null), site_tables)
if (length(site_tables) == 0) {
  stop("No requested iterations produced readable DMS results.")
}

all_sites <- do.call(rbind, site_tables)
site_keys <- paste(all_sites$chr, all_sites$start, all_sites$end, sep = "\t")
site_frequency <- as.data.frame(table(site_keys), stringsAsFactors = FALSE)
names(site_frequency) <- c("key", "iterations_present")
site_frequency$iterations_present <- as.integer(site_frequency$iterations_present)
stable_keys <- site_frequency$key[
  site_frequency$iterations_present / length(iterations) > 0.95
]
if (length(stable_keys) == 0) {
  stop("No DMS were present in more than 95% of all requested iterations.")
}

first_occurrence <- !duplicated(site_keys)
unique_sites <- all_sites[first_occurrence, , drop = FALSE]
unique_keys <- site_keys[first_occurrence]
stable_sites <- unique_sites[match(stable_keys, unique_keys), , drop = FALSE]
stable_sites$iterations_present <- site_frequency$iterations_present[
  match(stable_keys, site_frequency$key)
]
stable_sites$frequency <- stable_sites$iterations_present / length(iterations)

results_dir <- file.path(workdir, "output", "MLresults")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
results_file <- file.path(results_dir, "stable_DMS_gt95pct.csv")
write.csv(stable_sites, results_file, row.names = FALSE)

message(
  "Selected ", nrow(stable_sites), " stable DMS from ",
  length(iterations), " requested iterations. Saved: ", results_file
)
