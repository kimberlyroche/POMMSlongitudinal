# This script fits the (fido::pibble) model.

library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(driver)
library(fido)

source("00_functions.R")

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

filter_CIs <- function(Sigma, correlators, threshold) {
  df <- data.frame(tag1 = c(),
                   tag2 = c(),
                   pair_name = c(),
                   left = c(),
                   middle = c(),
                   right = c())
  for(idx in correlators) {
    pair <- combos[,idx]
    x <- Sigma[pair[1],pair[2],]
    # Filter to a large mean
    mu <- mean(x)
    if((sign(threshold) == 1 & mu > threshold) |
       (sign(threshold) == -1 & mu < threshold)) {
      pair_name <- paste0(get_tax_label(pair[1], tax),
                          " x ",
                          get_tax_label(pair[2], tax))
      df <- rbind(df, data.frame(tag1 = pair[1],
                                 tag2 = pair[2],
                                 pair_name = pair_name,
                                 left = quantile(x, 0.025)[[1]],
                                 middle = mu,
                                 right = quantile(x, 0.975)[[1]]))
    }
  }
  return(df)
}

# ------------------------------------------------------------------------------
#   Parse data
# ------------------------------------------------------------------------------

filter_akkermansia_subjects <- FALSE
filter_akkermansia_results <- TRUE
baseline_only <- FALSE

data_obj <- readRDS(file.path("data", "processed_data.rds"))
counts <- data_obj$counts
tax <- data_obj$tax
metadata <- data_obj$metadata
rm(data_obj)

subject_tallies <- table(metadata$ind_id)
subjects <- as.numeric(names(subject_tallies)[which(subject_tallies >= 4)])

akkermansia_list <- file.path("data", "akkermansia_subjects.rds")
if(filter_akkermansia_subjects & file.exists(akkermansia_list)) {
  akkermansia_subjects <- readRDS(akkermansia_list)
  subjects <- intersect(subjects, akkermansia_subjects)
}

if(baseline_only) {
  # Pull first visit only
  retain_samples <- metadata[(metadata$visit == 1 & metadata$ind_id %in% subjects),]$sample_id
  counts <- counts[,colnames(counts) %in% retain_samples]
  metadata <- metadata[metadata$sample_id %in% retain_samples,]
}

# ------------------------------------------------------------------------------
#   Pull existing metabolic metadata
# ------------------------------------------------------------------------------

met <- read.table("data/MET-S-nonCook_7-23-20pedsobesity_r24_phenotype_wide_20200723.tsv",
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE)

# Baseline BMI95
BMI_P95 <- sapply(subjects,
                  function(x) unname(met[met$IND == x,]$BMI_P95))
names(BMI_P95) <- subjects

# Change in BMI95 and HOMA data
# delta_BMI <- sapply(subjects,
#                     function(x) unname(met[met$IND == x,]$Delta_BMI_P95))
# names(delta_BMI) <- subjects
# PC_HOMA_IR <- sapply(subjects,
#                      function(x) unname(met[met$IND == x,]$PC_HOMA_IR))
# names(PC_HOMA_IR) <- subjects

# Restrict to subjects with extant BMI95 values. It's not clear how best to
# impute these.
selected_subjects <- subjects[which(!is.na(BMI_P95))]

# Center and scale these parameters. This bungles interpretation but I don't
# think we want to directly interpret the regression coefficients for these
# anyway.

BMI_P95 <- BMI_P95[names(BMI_P95) %in% selected_subjects]
saved_names <- names(BMI_P95)
BMI_P95 <- as.vector(scale(BMI_P95, center = FALSE, scale = TRUE))
names(BMI_P95) <- saved_names

# delta_BMI <- delta_BMI[names(delta_BMI) %in% selected_subjects]
# saved_names <- names(delta_BMI)
# delta_BMI <- as.vector(scale(delta_BMI, center = FALSE, scale = TRUE))
# names(delta_BMI) <- saved_names
# 
# PC_HOMA_IR <- PC_HOMA_IR[names(PC_HOMA_IR) %in% selected_subjects]
# saved_names <- names(PC_HOMA_IR)
# PC_HOMA_IR <- as.vector(scale(PC_HOMA_IR, center = FALSE, scale = TRUE))
# names(PC_HOMA_IR) <- saved_names

ggplot(data.frame(x = BMI_P95), aes(x = x)) +
  geom_histogram(bins = 10, fill = "gray", color = "black") +
  xlab("Baseline BMI95 distribution")

# ------------------------------------------------------------------------------
#   Fit model
# ------------------------------------------------------------------------------

testing <- FALSE

if(!file.exists(file.path("output", "fitted_model.rds"))) {

  subjects_to_use <- selected_subjects
  if(testing) {
    subjects_to_use <- selected_subjects[1:10]
  }
  
  subject_labels <- c()
  bmi_covariate <- c()
  Y <- NULL
  for(subject in subjects_to_use) {
    subject_counts <- counts[,colnames(counts) %in% metadata[metadata$ind_id == subject,]$sample_id,drop=FALSE]
    if(is.null(Y)) {
      Y <- subject_counts
    } else {
      Y <- cbind(Y, subject_counts)
    }
    subject_labels <- c(subject_labels, rep(subject, ncol(subject_counts)))
    baseline_scaled_bmi <- unname(BMI_P95[which(names(BMI_P95) == subject)])
    bmi_covariate <- c(bmi_covariate,
                       rep(baseline_scaled_bmi, ncol(subject_counts)))
  }
  
  # Build design matrix
  subject_labels <- as.factor(subject_labels)
  if(baseline_only) {
    X <- t(model.matrix(~bmi_covariate))
  } else {
    X <- t(model.matrix(~subject_labels + bmi_covariate))
    # X <- matrix(1, 1, ncol(Y)) # intercept-only model
  }
  
  # Filter out taxa totally absent in this subset of subject samples
  Y <- Y[rowSums(Y) > 0,]
  Y <- as.matrix(Y)
  
  # Set priors a la
  # https://jsilve24.github.io/fido/articles/introduction-to-fido.html
  upsilon <- nrow(Y) + 3
  Omega <- diag(nrow(Y))
  G <- cbind(diag(nrow(Y)-1), -1)
  Xi <- (upsilon-nrow(Y))*G%*%Omega%*%t(G)
  alr_Y <- alr_array(Y + 1, parts = 1)
  Theta <- matrix(rowMeans(alr_Y), nrow(Y)-1, nrow(X))
  Gamma <- diag(nrow(X))
  
  # Fit model and convert to CLR
  cat("Fitting model...\n")
  fit <- pibble(as.matrix(Y),
                X,
                upsilon,
                Theta,
                Gamma,
                Xi,
                n_samples = 1000,
                ret_mean = TRUE)
  
  cat("Model fit!\n")
  fit.clr <- to_clr(fit)
  
  saveRDS(list(fit = fit, fit.clr = fit.clr, subjects = subject_labels),
          file = file.path("output", "fitted_model.rds"))
  quit()
} else {
  fit_obj <- readRDS(file.path("output", "fitted_model.rds"))
  fit <- fit_obj$fit
  fit.clr <- fit_obj$fit.clr
  subjects_to_use <- fit_obj$subjects
}

# ------------------------------------------------------------------------------
#   Visualize taxon-taxon correlation matrix
# ------------------------------------------------------------------------------

# Scale Sigma to correlation
Sigma_corr <- fit.clr$Sigma
for(i in 1:fit$iter) {
  Sigma_corr[,,i] <- cov2cor(Sigma_corr[,,i])
}

# Get mean posterior predictions for a few parameters
mean_Sigma <- apply(Sigma_corr, c(1,2), mean)

# Sigma visualization
show_diag <- FALSE
temp <- mean_Sigma
if(show_diag) {
  temp[upper.tri(temp, diag = TRUE)] <- NA
}
temp <- as.data.frame(cbind(1:nrow(temp), temp))
colnames(temp) <- c("taxon1", 1:nrow(temp))
temp <- pivot_longer(temp,
                     !taxon1,
                     names_to = "taxon2",
                     values_to = "correlation")
temp$taxon2 <- as.numeric(temp$taxon2)
if(show_diag) {
  temp <- temp[complete.cases(temp),]
}

ggplot(temp, aes(x = taxon1, y = taxon2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", high = "darkred")
ggsave(file.path("output", "images", "Sigma.png"), units = "in", dpi = 100, height = 6, width = 7.5)

# ------------------------------------------------------------------------------
#   Filter to non-zero-spanning posterior intervals (if you have a posterior!!!)
# ------------------------------------------------------------------------------

combos <- combn(1:nrow(mean_Sigma), m = 2)

if(filter_akkermansia_results) {
  akkermansia_idx <- which(tax[,6] == "Akkermansia")
  ak_cols <- which(combos[1,] == akkermansia_idx | combos[2,] == akkermansia_idx)
  net_sign <- sapply(1:ncol(combos), function(x) {
    if(x %in% ak_cols) {
      sum(sign(quantile(Sigma_corr[combos[1,x],combos[2,x],], probs = c(0.025, 0.975))))
    } else {
      0
    }
  })
} else {
  # Threshold to correlators with non-zero posterior 95% CIs
  net_sign <- sapply(1:ncol(combos), function(x) {
    pair <- combos[,x]
    sum(sign(quantile(Sigma_corr[pair[1],pair[2],], probs = c(0.025, 0.975))))
  })
}

negative_correlators <- which(net_sign < 0)
positive_correlators <- which(net_sign > 0)

# Plot distribution of negative correlators in order to eyeball threshold...
if(filter_akkermansia_results) {
  par(mfrow = c(1,2))
  test <- sapply(negative_correlators, function(x) {
    Sigma_corr[combos[1,x],combos[2,x],]
  })
  hist(c(test))
  test <- sapply(positive_correlators, function(x) {
    Sigma_corr[combos[1,x],combos[2,x],]
  })
  hist(c(test))
  par(mfrow = c(1,1))
}

threshold <- 0.5
if(filter_akkermansia_results) {
  threshold <- 0.15
}
if(baseline_only) {
  threshold <- 0.5
}

df_neg <- filter_CIs(Sigma_corr, negative_correlators, threshold = -threshold)
df_pos <- filter_CIs(Sigma_corr, positive_correlators, threshold = threshold)

if(nrow(df_neg) > 0) {
  write.table(df_neg,
              file = file.path("output", "negative_correlators.tsv"),
              sep = "\t")
}

if(nrow(df_pos) > 0) {
  write.table(df_pos,
              file = file.path("output", "positive_correlators.tsv"),
              sep = "\t")
}

# Diagnostic plotting of "hits"
df_plot <- df_neg # choose pairs to plot
for(i in 1:nrow(df_plot)) {
  x <- as.vector(fit.clr$Eta[df_plot[i,]$tag1,,1] - fit.clr$Lambda[df_plot[i,]$tag1,,1]%*%fit.clr$X)
  y <- as.vector(fit.clr$Eta[df_plot[i,]$tag2,,1] - fit.clr$Lambda[df_plot[i,]$tag2,,1]%*%fit.clr$X)
  tax1 <- get_tax_label(df_plot[i,]$tag1, tax)
  tax2 <- get_tax_label(df_plot[i,]$tag2, tax)
  if(exists("subject_labels") && length(subject_labels) == length(x)) {
    p <- ggplot(data.frame(x = x, y = y, subject = subject_labels), aes(x = x, y = y, color = subject)) +
      scale_color_manual(values = palette)
  } else {
    p <- ggplot(data.frame(x = x, y = y), aes(x = x, y = y))
  }
  p <- p +
    geom_point(size = 2) +
    xlab(paste0("CLR(",tax1,")")) +
    ylab(paste0("CLR(",tax2,")")) +
    theme_bw() +
    theme(legend.position = "none")
  show(p)
  ggsave(file.path("output", "images", paste0("check_",i,".png")),
         p,
         units = "in",
         dpi = 100,
         height = 4,
         width = 6)
}
