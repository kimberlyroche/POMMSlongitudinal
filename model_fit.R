library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(driver)
library(fido)

setwd("C:/Users/kim/Documents/POMMSlongitudinal/")

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

# -------------------------------------------------------------------------------------------------
#   Parse data
# -------------------------------------------------------------------------------------------------

data <- readRDS("data/phyloseq_r24_complete_16S_metadata_corrected.rds")
metadata <- sample_data(data)
counts <- otu_table(data)@.Data # 960 samples x ~151K taxa
tax <- tax_table(data)@.Data
rownames(tax) <- NULL # kill the huge sequence labels

# -------------------------------------------------------------------------------------------------
#   Filter
# -------------------------------------------------------------------------------------------------

# Filter to cohorts we care about: OB and HWC

retain_samples <- metadata[metadata$group %in% c("OB", "HWC"),]$sample_id
counts <- counts[rownames(counts) %in% retain_samples,]
metadata <- metadata[metadata$sample_id %in% retain_samples,]

cat("Starting with",nrow(counts),"samples x",ncol(counts),"taxa\n")

# Collapse taxa below a minimum relative abundance to an "other" category

props <- t(apply(counts, 1, function(x) x/sum(x)))
max_relative_abundance <- apply(props, 2, max)
retain_taxa <- max_relative_abundance >= 0.05
cat(paste0("Retaining ", round((sum(retain_taxa) / length(max_relative_abundance))*100, 1),
           " percent of taxa (",sum(retain_taxa),")\n"))

retained_total <- sum(counts[,retain_taxa])
total <- sum(counts)
cat("Retained taxa account for",round((retained_total / total)*100, 1),"percent of total counts\n")

counts_collapsed <- cbind(counts[,retain_taxa], rowSums(counts[,!retain_taxa]))
counts <- counts_collapsed
cat("Collapsed to",nrow(counts_collapsed),"samples x",ncol(counts_collapsed),"taxa\n")
tax <- tax[retain_taxa,]
tax <- rbind(tax, NA)

# -------------------------------------------------------------------------------------------------
#   Fit model to 10 subjects with 5 samples (test)
# -------------------------------------------------------------------------------------------------

data_OB <- get_cohort_samples(cohort = "OB")

subject_tallies <- table(data_OB$subjects)
subjects <- as.numeric(names(subject_tallies)[which(subject_tallies == 5)])
subject_labels <- c()
Y <- NULL
for(subject in subjects[1:10]) {
  subject_counts <- data_OB$counts[which(data_OB$subjects == subject),]
  if(is.null(subject_counts)) {
    Y <- subject_counts
  } else {
    Y <- rbind(Y, subject_counts)
  }
  subject_labels <- c(subject_labels, rep(subject, nrow(subject_counts)))
}
Y <- t(Y) # D x N
rownames(Y) <- NULL

# Build design matrix
subject_labels <- as.factor(subject_labels)
X <- t(model.matrix(~subject_labels))
# X <- matrix(1, 1, ncol(Y)) # intercept-only model

# Filter out taxa totally absent in this subset of subject samples
Y <- Y[rowSums(Y) > 0,]

# Set priors a la
# https://jsilve24.github.io/fido/articles/introduction-to-fido.html
upsilon <- nrow(Y) + 3
Omega <- diag(nrow(Y))
G <- cbind(diag(nrow(Y)-1), -1)
Xi <- (upsilon-nrow(Y))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y)-1, nrow(X))
Gamma <- diag(nrow(X))

# Fit model and convert to CLR
fit <- pibble(Y, X, upsilon, Theta, Gamma, Xi) # takes about 5 sec. at 116 taxa x 15 samples
fit.clr <- to_clr(fit)

# -------------------------------------------------------------------------------------------------
#   SANITY CHECK #1 - Do strong correlators seem plausible?
# -------------------------------------------------------------------------------------------------

# Scale Sigma to correlation
Sigma_corr <- fit.clr$Sigma
for(i in 1:fit$iter) {
  Sigma_corr[,,i] <- cov2cor(Sigma_corr[,,i])
}

# Get mean posterior predictions for a few parameters
mean_Sigma <- apply(Sigma_corr, c(1,2), mean)
mean_Eta <- apply(fit.clr$Eta, c(1,2), mean)
mean_Lambda <- apply(fit.clr$Lambda, c(1,2), mean)

# Find (and visualize) strong correlators
mean_Sigma[upper.tri(mean_Sigma, diag = TRUE)] <- NA
pairs <- which(mean_Sigma > 0.75, arr.ind = TRUE)
image(mean_Sigma)

# Sample a pair of correlated CLR taxa
sample_row <- sample(1:nrow(pairs), size = 1)
idx1 <- pairs[sample_row,1]
idx2 <- pairs[sample_row,2]

# Visualize the (apparently correlated) residuals
set1 <- mean_Eta[idx1,] - mean_Lambda[idx1,]
set2 <- mean_Eta[idx2,] - mean_Lambda[idx2,]
plot(set1, set2)
cor(set1, set2)
cat("Pair:",idx1,"x",idx2,"\n")

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
                                     matched_sign = sign(q1) == sign(q2)) %>%
                           filter(matched_sign)) %>%
                           select(pair)
head(results)

# -------------------------------------------------------------------------------------------------
#   TO DO
# -------------------------------------------------------------------------------------------------

# (1) Filter based on OB cohort only; we're not interested in HWC samples here
# (2) Prior choices?
# (3) With a few subjects the "rare" group has high correlation because of total absense in some
#       subjects and relatively high abundance in others. If this doens't disappear with more
#       subjects, this probably needs a re-think.
# (4) Fit across many subjects and build "rug" plot



