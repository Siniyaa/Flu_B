#!/usr/bin/env Rscript
# Firth-penalized logistic-regression sensitivity analysis for separated binary data.
#
# Usage:
#   Rscript scripts/firth_logistic_sensitivity.R input.csv output.txt
#
# Expected columns in input.csv:
#   focal_NA_any        Binary outcome: 1 if at least one focal NA substitution is present.
#   HA_epitope_burden   Numeric HA epitope-substitution count.
#   clade               Clade label.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript firth_logistic_sensitivity.R input.csv output.txt")
}

input_file <- args[1]
output_file <- args[2]

dat <- read.csv(input_file, stringsAsFactors = TRUE)
required <- c("focal_NA_any", "HA_epitope_burden", "clade")
missing_cols <- setdiff(required, names(dat))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

if (!requireNamespace("logistf", quietly = TRUE)) {
  stop("Package 'logistf' is required. Install it with install.packages('logistf').")
}

dat$high_HA_epitope_burden <- as.integer(dat$HA_epitope_burden > 5)

fit <- logistf::logistf(
  focal_NA_any ~ high_HA_epitope_burden + clade,
  data = dat
)

sink(output_file)
cat("Firth-penalized logistic-regression sensitivity analysis\n")
cat("=======================================================\n\n")
print(summary(fit))
sink()
