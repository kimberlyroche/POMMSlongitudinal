library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(driver)
library(fido)

# -------------------------------------------------------------------------------------------------
#   Functions
# -------------------------------------------------------------------------------------------------

# Pull per-group counts and metadata
get_cohort_samples <- function(cohort = "OB") {
  sample_ids <- metadata[metadata$group == cohort,]$sample_id
  ind_ids <- metadata[metadata$group == cohort,]$ind_id
  visits <- metadata[metadata$group == cohort,]$visit
  tags <- paste0(ind_ids, "_", visits)
  reorder_idx <- order(tags)
  ind_ids <- ind_ids[reorder_idx]
  visits <- visits[reorder_idx]
  cohort_counts <- counts[rownames(counts) %in% sample_ids,]
  cohort_counts <- cohort_counts[reorder_idx,]
  return(list(cohort = cohort, counts = cohort_counts, subjects = ind_ids, visits = visits))
}

generate_highcontrast_palette <- function(S) {
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  sample(getPalette(S))
}

get_tax_label <- function(taxon_idx, tax_map) {
  tax_pieces <- tax_map[taxon_idx,]
  level_idx <- max(which(!is.na(tax_pieces)))
  paste0(names(tax_pieces[level_idx]), " ", as.character(tax_pieces[[level_idx]]))
}

# `data` is presumed to be taxa x samples
# The last row of tax is assumed to be an "other" category (all <NA>)
collapse_below_minimum <- function(data, tax_map, threshold) {
  if(threshold < 1 & max(colSums(data)) > 1) {
    # Convert to proportions if that seems appropriate
    data <- apply(data, 2, function(x) x/sum(x))
  }
  retain_taxa <- apply(data, 1, max) >= threshold
  retain_taxa[length(retain_taxa)] <- FALSE
  data <- rbind(data[retain_taxa,], colSums(data[!retain_taxa,]))
  retain_taxa[length(retain_taxa)] <- TRUE
  tax_map <- tax_map[retain_taxa,]
  return(list(data = data, tax = tax_map))
}

# -------------------------------------------------------------------------------------------------
#   Parse data
# -------------------------------------------------------------------------------------------------

data_obj <- readRDS("processed_data.rds")
counts <- data_obj$counts
tax <- data_obj$tax
metadata <- data_obj$metadata
rm(data_obj)

# -------------------------------------------------------------------------------------------------
#   Fit model to 10 subjects with 5 samples (test)
# -------------------------------------------------------------------------------------------------

subject_tallies <- table(metadata$ind_id)
subjects <- as.numeric(names(subject_tallies)[which(subject_tallies >= 4)])
subject_labels <- c()
Y <- NULL
for(subject in subjects) {
# for(subject in subjects[1:5]) {
  subject_counts <- counts[,colnames(counts) %in% metadata[metadata$ind_id == subject,]$sample_id]
  if(is.null(Y)) {
    Y <- subject_counts
  } else {
    Y <- cbind(Y, subject_counts)
  }
  subject_labels <- c(subject_labels, rep(subject, ncol(subject_counts)))
}

# Build design matrix
subject_labels <- as.factor(subject_labels)
X <- t(model.matrix(~subject_labels))
# X <- matrix(1, 1, ncol(Y)) # intercept-only model

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
# Theta <- matrix(0, nrow(Y)-1, nrow(X))
Gamma <- diag(nrow(X))*2

# Fit model and convert to CLR
cat("Fitting model...\n")
# fit <- pibble(Y, X, upsilon, Theta, Gamma, Xi, n_samples = 500) # takes about 5 sec. at 116 taxa x 15 samples
fit <- pibble(Y, X, upsilon, Theta, Gamma, Xi, n_samples = 500,
              b2 = 0.98, step_size = 0.002, eps_f = 1e-11, eps_g = 1e-05,
              max_iter = 10000L, optim_method = "adam")

cat("Model fit!\n")
fit.clr <- to_clr(fit)

# 4 cores -- all failed
#  fitted_model: pibble w/ LBFGS
#  fitted_model_2: pibble w/ Adam (b = 0.99, step_size = 0.004)
#  fitted_model_3: pibble w/ Adam (as above) and larger covariance Gamma (to accommodate heterogeneity across users?)
# 1 core
#  fitted_model_4: pibble w/ Adam (b = 0.98, step_size = 0.002)
#   job 24319681

saveRDS(list(fit = fit, fit.clr = fit.clr), file = "fitted_model_4.rds")
quit()

# -------------------------------------------------------------------------------------------------
#   RESULTS: TAKE A LOOK AT STRONG CORRELATORS
# -------------------------------------------------------------------------------------------------

fit_obj <- readRDS("fitted_model.rds")
fit <- fit_obj$fit
fit.clr <- fit_obj$fit.clr
subject_labels <- fit_obj$subjects
rm(fit_obj)

# Scale Sigma to correlation
Sigma_corr <- fit.clr$Sigma
for(i in 1:fit$iter) {
  Sigma_corr[,,i] <- cov2cor(Sigma_corr[,,i])
}

# Get mean posterior predictions for a few parameters
mean_Sigma <- apply(Sigma_corr, c(1,2), mean)
mean_Eta <- apply(fit.clr$Eta, c(1,2), mean)

image(mean_Sigma) # TBD: replace with decent ggplot
image(Sigma_corr[,,sample(1:500, size = 1)])

# -------------------------------------------------------------------------------------------------
#   FILTER TO NON-ZERO-SPANNING POSTERIOR INTERVALS
# -------------------------------------------------------------------------------------------------

# Test with a submatrix (this is slow)
sub_mat <- Sigma_corr[1:20,1:20,]
vectorized_mat <- NULL
for(i in 1:dim(sub_mat)[3]) {
  temp <- sub_mat[,,i]
  temp[upper.tri(temp, diag = TRUE)] <- NA
  temp <- c(temp)
  temp <- temp[!is.na(temp)]
  if(is.null(vectorized_mat)) {
    vectorized_mat <- temp
  } else {
    vectorized_mat <- rbind(vectorized_mat, temp)
  }
}
rownames(vectorized_mat) <- NULL

# Set column names to the form "2_1" or "23_17" denoting taxon pairs
pairs <- t(combn(1:nrow(sub_mat), m = 2))[,c(2,1)]
pair_labels <- paste0(pairs[,1],"_",pairs[,2])
colnames(vectorized_mat) <- pair_labels

pairs_long <- pivot_longer(as.data.frame(vectorized_mat), everything(), names_to = "pair", values_to = "correlation")

# Filter to 95% CIs that don't span zero
results <- as.data.frame(pairs_long %>%
                           group_by(pair) %>%
                           summarize(q1 = quantile(correlation, probs = c(0.025)),
                                     q2 = quantile(correlation, probs = c(0.975)),
                                     matched_sign = sign(q1) == sign(q2), .groups = 'drop') %>%
                           filter(matched_sign)) %>%
                           select(pair)
head(results)

