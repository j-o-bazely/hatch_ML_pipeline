###          C3. PQLseq P Value adjustment with SLIM                ###

## Date: 06.10.2026

##  Author: James B, adapted from Charley Y

## Description:
## This code runs multiple testing correction on PQLseq-generated
## P values using SLIM.


## Pull iteration, available cores and chromosome ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
num_core <- as.numeric(args[2])

## Set paths
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")


####    Merge all chromosomes   ####

## Create chromosome vector to loop through ##
chroms <- strsplit(
  Sys.getenv(
    "CH2_CHROMOSOMES",
    unset = paste(c("chr0", paste0("chr", 1:28)), collapse = ",")
  ),
  ","
)[[1]]

## Create emopty lists to store results ##
list_fit <- list()
list_dfs_fit <- list()

## Load in data for each chromosome ##
for (chr in chroms) {
  print(chr)
  list_fit[[chr]] <- readRDS(
    file.path(
      tmpdir,
      paste0(
        "hatch_trainingSamples_PQLseqOutput_Fit_iter",
        iter,
        "_SLK063_ragtag_",
        chr,
        ".RDS"
      )
    )
  )
  list_dfs_fit[[chr]] <- list_fit[[chr]]
}

##  Bind chromosomes into a single dataframe ##
df_all_chrms <- do.call(rbind, list_dfs_fit)

##  Check the % of NaNs for h2  ##
length(df_all_chrms$h2[grep("NaN", df_all_chrms$h2)]) / nrow(df_all_chrms) * 100

##  Save object  ##
saveRDS(
  df_all_chrms,
  file.path(
    tmpdir,
    paste0("hatch_trainingSamples_PQLseqFit_iter", iter, "_wholeGenome.RDS")
  )
)


####    Run SLIM correction   ####


# SLIMcodes.R downloaded from: https://www.aspendb.org/downloads

#### Define SLIM functions ####

# SLIMcodes.R uses the same defaults as methylKit.
slim_func <- function(
  raw_p,
  sta = .1,
  divi = 10,
  p_z = 0.05,
  b = 100,
  b_plot = FALSE
) {
  ####################
  m <- length(raw_p) # nolint # nolint

  ########################
  alpha_mtx <- NULL #
  pi0s_est_com <- NULL
  sz_cluster_mtx <- NULL # nolint
  p_pi1_mtx <- NULL
  pi1_act_mtx <- NULL
  dist_group <- NULL # nolint
  num_mtx <- NULL
  gamma_values <- NULL
  pi0_matrix <- NULL

  #############
  ## observed points
  lambda_ga <- seq(0, 1, 0.001)
  gamma_ga <- sapply(lambda_ga, f1, raw_p = raw_p)
  gamma_values <- c(gamma_values, gamma_ga)
  alpha_mtx <- c(alpha_mtx, gamma_ga[which(lambda_ga == 0.05)])


  ### esimation
  pi0_mtx <- NULL
  x_axis <- NULL # nolint
  y_axis <- NULL # nolint
  interval <- (1 - sta) / divi
  for (i in 1:divi) ## 10 for uniform data
  { # nolint
    cutoff <- sta + (i / divi) * (1 - sta) ## 10 for uniform data
    lambda <- seq(cutoff - interval, cutoff, interval / 10)
    gamma_mtx <- sapply(lambda, f1, raw_p = raw_p) # nolint
    lm_model <- lm(gamma_mtx ~ lambda)
    pi0_mtx <- c(pi0_mtx, coefficients(lm_model)[2])
  }


  ##################################
  ######## searching
  n_com <- NULL
  n_raw_p <- NULL # nolint
  max_fdr_mtx <- NULL
  qua_point_mtx <- NULL
  if (b <= 1) b <- 100
  qua_point_mtx <- seq(0.01, 0.99, 1 / b)
  for (k in 1:length(qua_point_mtx)) # nolint
  { # nolint
    qua_point <- qua_point_mtx[k]

    pi0_comb_lr <- min(quantile(pi0_mtx, qua_point), 1)
    ## 0.4 for uniform or normal distribution;
    pi0_est <- pi0_comb_lr

    ########### Calculate independent index of raw p vlaues
    pi0_matrix <- rbind(pi0_matrix, pi0_mtx)

    pi0s_est_com <- c(pi0s_est_com, pi0_est)
    ## Condition1
    p_pi1 <- sort(raw_p)[max(length(raw_p) * (1 - pi0_est), 1)]
    p_pi1_mtx <- c(p_pi1_mtx, p_pi1)

    pi0 <- pi0_est
    if (is.null(p_z)) p_z <- 0.05
    max_fdr <- p_z * pi0 / (1 - (1 - p_z) * pi0)
    max_fdr_mtx <- c(max_fdr_mtx, max_fdr)

    q_values_comb_lr <- q_values_fun(raw_p, pi0)
    q_value_cutoff <- max_fdr
    selected <- which(q_values_comb_lr < q_value_cutoff)
    selected_q_values_comb_lr <- selected

    pi1_act_mtx <- c(
      pi1_act_mtx,
      length(selected_q_values_comb_lr) / length(raw_p)
    )
    n_com <- c(n_com, list(selected_q_values_comb_lr))
    num_mtx <- c(num_mtx, length(selected_q_values_comb_lr))
  }
  length(n_com)
  length(qua_point_mtx)

  #### doing judging
  ## by max FDR
  pi1s_est_com <- 1 - pi0s_est_com # nolint
  difference <- sum(raw_p <= p_z) / length(raw_p) - pi1_act_mtx

  ###
  location <- which.min(abs(difference))
  difference_location <- difference[location]
  selected_quantile <- qua_point_mtx[location]
  pi0_est <- min(1, pi0s_est_com[location])
  max_fdr_p_z <- p_z * pi0_est / (1 - (1 - p_z) * pi0_est)

  if (b_plot) {
    windows() # nolint
    par(mfrow = c(1, 2))
    hist(raw_p, main = "Histogram of p-value")
    gamma_ga <- sapply(lambda_ga, f1, raw_p = raw_p)
    plot(
      lambda_ga,
      gamma_ga,
      type = "l",
      main = "Relationship of p- and q value",
      xlab = expression(lambda),
      ylab = expression(gamma),
      cex.lab = 1.45,
      cex.axis = 1.42
    )
    # par(xaxp=c(0,1,10)); # nolint
    # axis(1); # nolint
    ## qvalues
    q_values <- q_values_fun(raw_p, pi0 = pi0_est)
    gamma_q_ga <- sapply(lambda_ga, f1, raw_p = q_values)
    lines(lambda_ga, gamma_q_ga, col = "blue", lwd = 2, lty = "dashed")
    abline(v = p_z, col = "black", lwd = 2, lty = "dotdash")
    abline(v = max_fdr_p_z, col = "blue", lwd = 2, lty = "dotdash")
    text(
      0.75,
      0.6,
      labels = paste("L=", round(abs(difference_location), 4), sep = "")
    )
    legend_labels <- list(
      bquote("CPD of p-value"),
      bquote("CPD of q-value"),
      bquote("Pmax" == .(p_z)),
      bquote("FDRmax" == .(round(max_fdr_p_z, 2)))
    )
    legend(
      "bottomright",
      legend = as.expression(legend_labels),
      lwd = 2,
      lty = c("solid", "dashed", "dotdash", "dotdash"),
      col = c("black", "blue", "black", "blue")
    )
  }

  return(list(pi0_est = pi0_est, selected_quantile = selected_quantile)) # nolint
}
f1 <- function(cutoff, raw_p) {
  sum(raw_p < cutoff) / length(raw_p) * 100
}
q_values_fun <- function(raw_p, pi0) {
  ordered_raw_p <- sort(raw_p)
  q_values <- pi0 * length(ordered_raw_p) * ordered_raw_p /
    seq_along(ordered_raw_p)
  temporary_q_values <- cummin(q_values[seq(length(q_values), 1, -1)])
  q_values <- temporary_q_values[seq(length(temporary_q_values), 1, -1)]
  q_values[order(order(raw_p))]
}

#### Run SLIM on PQLseq output ####

## Load PQLseq whole genome object generated previously ##
df_all_chrms <- readRDS(
  file.path(
    tmpdir,
    paste0("hatch_trainingSamples_PQLseqFit_iter", iter, "_wholeGenome.RDS")
  )
)

## Pull p values ##
pvals <- df_all_chrms$pvalue

# Run SLIM on these p values ##
qvals <- {
  q_values_fun(
    pvals,
    slim_func(
      pvals,
      sta = .1,
      divi = 10,
      p_z = 0.05,
      b = 100,
      b_plot = FALSE
    )$pi0_est
  )
}

## Add calculated Q values to PQLseq output ##
slim_output <- df_all_chrms
slim_output$qvalue <- qvals

####    Checks    ####
## Check number of sites which failed to converge ##
algorithm_failure <- nrow(
  slim_output[
    slim_output$qvalue == 0 & slim_output$converged == "FALSE",
  ]
)

message(
  paste0(algorithm_failure, " sites failed to converge and will be removed")
)

## Filter out failed convergences
slim_output_filtered <- slim_output[slim_output$converged == "TRUE", ]

## Save filtered object
saveRDS(
  slim_output_filtered,
  file.path(
    tmpdir,
    paste0(
      "hatch_trainingSamples_PQLseqFit_iter",
      iter,
      "_wholeGenome_corrected.RDS"
    )
  )
)
