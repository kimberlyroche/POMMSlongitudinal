# This script fits the (fido::pibble) model.

library(driver)
library(fido)
library(uuid)

source("00_functions.R")

# ------------------------------------------------------------------------------
#   Global vars
# ------------------------------------------------------------------------------

# See README.txt for a list of possible scenarios here.

# Use baseline samples only.
baseline_only <- FALSE

# Fit model on permutation of data.
permutation_test <- TRUE

# Use only high Akkermansia subjects.
filter_akkermansia_subjects <- FALSE

# If TRUE, includes P95 BMI as a covariate in the model, theoretically to
# suppress the effect of variation driven by differences in subject BMI.
# As of 2021-03-21, I've not evaluated how effective this is.
include_BMI <- TRUE

# If TRUE, subsets to 5 subjects so the model fits quickly. This is just for
# testing the effects of various parameterizations.
testing <- FALSE

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

fit_model <- function(subjects_to_use, Y, X, metadata,
                      permute = FALSE, MAP = FALSE) {
  # Set priors a la
  # https://jsilve24.github.io/fido/articles/introduction-to-fido.html
  upsilon <- nrow(Y) + 3
  Omega <- diag(nrow(Y))
  G <- cbind(diag(nrow(Y)-1), -1)
  Xi <- (upsilon-nrow(Y))*G%*%Omega%*%t(G)
  alr_Y <- alr_array(Y + 1, parts = 1)
  Theta <- matrix(rowMeans(alr_Y), nrow(Y)-1, nrow(X))
  Gamma <- diag(nrow(X))
  
  if(permute) {
    # Shuffle taxa randomly across observations.
    Y <- t(apply(Y, 1, function(x) sample(x)))
  }
  
  # Fit model and convert to CLR.
  n_samples <- 1000
  if(MAP) {
    n_samples <- 0
  }
  fit <- pibble(as.matrix(Y),
                X,
                upsilon,
                Theta,
                Gamma,
                Xi,
                n_samples = 0,
                ret_mean = TRUE)
  fit.clr <- to_clr(fit)
  return(list(fit = fit, fit.clr = fit.clr, subjects = subject_labels))
}

# ------------------------------------------------------------------------------
#   Parse data
# ------------------------------------------------------------------------------

data_obj <- readRDS(file.path("data", "processed_data.rds"))
counts <- data_obj$counts
tax <- data_obj$tax
metadata <- data_obj$metadata
rm(data_obj)

subject_tallies <- table(metadata$ind_id)
if(baseline_only) {
  subjects <- as.numeric(names(subject_tallies))
} else {
  # Filter to subjects with at least 4 visits; we want as many repeated
  # observations as we can get since these series are short.
  subjects <- as.numeric(names(subject_tallies)[which(subject_tallies >= 4)])
}

akkermansia_list <- file.path("data", "akkermansia_subjects.rds")
if(filter_akkermansia_subjects) {
  if(file.exists(akkermansia_list)) {
    akkermansia_subjects <- readRDS(akkermansia_list)
    subjects <- intersect(subjects, akkermansia_subjects)
  } else {
    stop(paste0("List of high Akkermansia subjects not found at: ",
                akkermansia_list,
                "\n"))
  }
}

if(baseline_only) {
  # Pull first visit only.
  retain_samples <- metadata[(metadata$visit == 1 & metadata$ind_id %in% subjects),]$sample_id
  counts <- counts[,colnames(counts) %in% retain_samples]
  metadata <- metadata[metadata$sample_id %in% retain_samples,]
}

# ------------------------------------------------------------------------------
#   Pull BMI data for inclusion as a covariate.
# ------------------------------------------------------------------------------

if(include_BMI) {
  met <- read.table("data/MET-S-nonCook_7-23-20pedsobesity_r24_phenotype_wide_20200723.tsv",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)
  BMI_P95 <- sapply(subjects,
                    function(x) unname(met[met$IND == x,]$BMI_P95))
  names(BMI_P95) <- subjects
  
  # Restrict to subjects with non-NA BMI95 values. Otherwise, it's not clear how
  # best to impute the missing ones here.
  selected_subjects <- subjects[which(!is.na(BMI_P95))]
  
  # Center and scale these parameters. This bungles interpretation but I don't
  # think we want to directly interpret the regression coefficients for these
  # anyway.
  BMI_P95 <- BMI_P95[names(BMI_P95) %in% selected_subjects]
  saved_names <- names(BMI_P95)
  BMI_P95 <- as.vector(scale(BMI_P95, center = FALSE, scale = TRUE))
  names(BMI_P95) <- saved_names
} else {
  selected_subjects <- subjects
}

# ------------------------------------------------------------------------------
#   Fit model
# ------------------------------------------------------------------------------

if(testing) {
  selected_subjects <- selected_subjects[1:5]
}

# Build observation matrix Y with samples from the selected subjects only.
subject_labels <- c()
if(include_BMI) {
  bmi_covariate <- c()
}
Y <- NULL
for(subject in selected_subjects) {
  subject_counts <- counts[,colnames(counts) %in% metadata[metadata$ind_id == subject,]$sample_id,drop=FALSE]
  if(is.null(Y)) {
    Y <- subject_counts
  } else {
    Y <- cbind(Y, subject_counts)
  }
  subject_labels <- c(subject_labels, rep(subject, ncol(subject_counts)))
  if(include_BMI) {
    baseline_scaled_bmi <- unname(BMI_P95[which(names(BMI_P95) == subject)])
    bmi_covariate <- c(bmi_covariate,
                       rep(baseline_scaled_bmi, ncol(subject_counts)))
  }
}
# Filter out taxa totally absent in this subset of subject samples. (It happens
# if we're subsetting -- e.g. for testing.)
Y <- Y[rowSums(Y) > 0,]
Y <- as.matrix(Y)

# Build covariate matrix X with subject labels, etc.
subject_labels <- as.factor(subject_labels)
if(include_BMI) {
  if(baseline_only) {
    X <- t(model.matrix(~ bmi_covariate))
  } else {
    X <- t(model.matrix(~ subject_labels + bmi_covariate))
  }
} else {
  if(baseline_only) {
    # Intercept-only model
    X <- matrix(1, 1, ncol(Y))
  } else {
    X <- t(model.matrix(~ subject_labels))
  }
}

# Fit the model.
fit_obj <- fit_model(subjects_to_use,
                     Y,
                     X,
                     metadata,
                     permute = permutation_test,
                     MAP = TRUE)
if(permutation_test) {
  saveRDS(fit_obj,
          file = file.path("output",
                           "fitted_models",
                           paste0("model_",
                                  UUIDgenerate(),
                                  ".rds")))
} else {
  saveRDS(fit_obj, file = file.path("output",
                                    "fitted_models",
                                    "model_canonical.rds"))
}
