library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(driver)
library(fido)

# setwd("C:/Users/kim/Documents/POMMSlongitudinal/")

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
collapse_below_minimum <- function(data, tax_map, threshold) {
  if(threshold < 1 & max(colSums(data)) > 1) {
    # Convert to proportions if that seems appropriate
    data <- apply(data, 2, function(x) x/sum(x))
  }
  retain_taxa <- apply(data, 1, max) >= threshold
  data <- rbind(data[retain_taxa,], colSums(data[!retain_taxa,]))
  tax_map <- tax_map[retain_taxa,]
  tax_map <- rbind(tax_map, NA)
  return(list(data = data, tax = tax_map))
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

# Filter to cohort we care about: OB

retain_samples <- metadata[metadata$group == "OB",]$sample_id
counts <- counts[rownames(counts) %in% retain_samples,]
metadata <- metadata[metadata$sample_id %in% retain_samples,]
colnames(counts) <- NULL

cat("Starting with",nrow(counts),"samples x",ncol(counts),"taxa\n")

# (1) Remove all-zero taxa (almost 1/4)
present_taxa <- which(colSums(counts) > 0)
counts <- counts[,present_taxa]
tax <- tax[present_taxa,]

# Show percent taxa missing each taxonomic label
# na_tax <- is.na(tax)
# names(na_tax) <- names(tax)
# apply(na_tax, 2, function(x) sum(x) / nrow(na_tax))

# (2) Collapse to highest common taxonomic level
agglomerated_counts <- NULL
for(tax_level in 1:6) {
  cat("Agglomerating at tax level:",colnames(tax)[tax_level],"\n")
  counts_subset <- as.data.frame(cbind(tax, t(counts)))
  counts_subset <- counts_subset[which(!is.na(counts_subset[,tax_level]) & is.na(counts_subset[,tax_level+1])),]
  copy_long <- pivot_longer(counts_subset,
                            !c("domain", "phylum", "class", "order", "family", "genus"),
                            names_to = "sample",
                            values_to = "count")
  if(nrow(copy_long) > 0) {
    copy_long$count <- as.numeric(copy_long$count)
    result <- copy_long %>%
      group_by(domain, phylum, class, order, family, genus, sample) %>%
      summarize(total_counts = sum(count), .groups = 'drop')
    counts_subset <- as.data.frame(pivot_wider(result, names_from = "sample", values_from = "total_counts"))
    if(is.null(agglomerated_counts)) {
      agglomerated_counts <- counts_subset
    } else {
      agglomerated_counts <- rbind(agglomerated_counts, counts_subset)
    }
  }
}

# The pivots will have disordered the samples; get them back in the original (subject x visit) order
orig_sample_order <- rownames(counts)
new_tax <- agglomerated_counts[,1:6]
new_counts <- agglomerated_counts[,7:ncol(agglomerated_counts)]
new_counts <- new_counts[,orig_sample_order]

counts <- new_counts
tax <- new_tax
# Note: `counts` and `tax` now have taxa as rows

# (3) Count features with max abundance < 1%
props <- apply(counts, 2, function(x) x/sum(x))
max_relative_abundance <- apply(props, 1, max)
retain_taxa <- max_relative_abundance >= 0.01
cat(paste0("Retaining ", round((sum(retain_taxa) / length(max_relative_abundance))*100, 1),
           " percent of taxa (",sum(retain_taxa),")\n"))

retained_total <- sum(counts[retain_taxa,])
total <- sum(counts)
cat("Retained taxa account for",round((retained_total / total)*100, 1),"percent of total counts\n")

counts_collapsed <- rbind(counts[retain_taxa,], rowSums(counts[!retain_taxa,]))
counts <- counts_collapsed
rownames(counts) <- NULL
cat("Collapsed to",nrow(counts),"taxa x",ncol(counts),"samples\n")
tax <- tax[retain_taxa,]
tax <- rbind(tax, NA)
rownames(tax) <- NULL

# Visualize these compositions quickly
# Pick a subset of the data that corresponds to visits 1-5 from subject 192
#   metadata[metadata$sample_id %in% colnames(counts)[15:19],c("ind_id","visit")]
# counts_subset <- counts[,15:19]

# props <- collapse_below_minimum(counts_subset, tax, threshold = 0.01)

# # Strip down to taxa
# tax_labels <- unname(sapply(1:(nrow(props$data)-1), function(x) get_tax_label(x, props$tax)))

# rownames(props$data) <- c(tax_labels, "assorted low abundance")
# props$data <- rbind(sample_index = 1:ncol(props$data), props$data)
# plot_data <- pivot_longer(as.data.frame(t(props$data)), !sample_index, names_to = "taxon", values_to = "relative_abundance")
# plot_data$taxon <- as.factor(plot_data$taxon)
# palette <- generate_highcontrast_palette(nrow(props$data))
# p <- ggplot(plot_data, aes(fill = taxon, y = relative_abundance, x = sample_index)) +
#   geom_bar(position = "stack", stat = "identity") +
#   scale_fill_manual(values = palette)
# p

# saveRDS(list(counts = counts, tax = tax, metadata = metadata), file = "processed_data.rds")

# -------------------------------------------------------------------------------------------------
#   Fit model to 10 subjects with 5 samples (test)
# -------------------------------------------------------------------------------------------------

# processed_data <- readRDS("processed_data.rds")
# counts <- processed_data$counts
# tax <- processed_data$tax
# metadata <- processed_data$metadata

subject_tallies <- table(metadata$ind_id)
subjects <- as.numeric(names(subject_tallies)[which(subject_tallies == 5)])
subject_labels <- c()
Y <- NULL
for(subject in subjects[1:5]) {
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

# Set priors a la
# https://jsilve24.github.io/fido/articles/introduction-to-fido.html
upsilon <- nrow(Y) + 3
Omega <- diag(nrow(Y))
G <- cbind(diag(nrow(Y)-1), -1)
Xi <- (upsilon-nrow(Y))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y)-1, nrow(X))
Gamma <- diag(nrow(X))

# Fit model and convert to CLR
fit <- pibble(as.matrix(Y), X, upsilon, Theta, Gamma, Xi) # takes about 5 sec. at 116 taxa x 15 samples
fit.clr <- to_clr(fit)

saveRDS(list(fit = fit, fit.clr = fit.clr), file = "fitted_model.rds")

# -------------------------------------------------------------------------------------------------
#   SANITY CHECK #1 - Do strong correlators seem plausible?
# -------------------------------------------------------------------------------------------------

# Scale Sigma to correlation
Sigma_corr <- fit.clr$Sigma
for(i in 1:fit$iter) {
  Sigma_corr[,,i] <- cov2cor(Sigma_corr[,,i])
}

# # Get mean posterior predictions for a few parameters
# mean_Sigma <- apply(Sigma_corr, c(1,2), mean)
# mean_Eta <- apply(fit.clr$Eta, c(1,2), mean)
# mean_Lambda <- apply(fit.clr$Lambda, c(1,2), mean)

# # Find (and visualize) strong correlators
# mean_Sigma[upper.tri(mean_Sigma, diag = TRUE)] <- NA
# pairs <- which(mean_Sigma > 0.75, arr.ind = TRUE)
# image(mean_Sigma)

# # Sample a pair of correlated CLR taxa
# sample_row <- sample(1:nrow(pairs), size = 1)
# idx1 <- pairs[sample_row,1]
# idx2 <- pairs[sample_row,2]

# # Visualize the (apparently correlated) residuals
# set1 <- mean_Eta[idx1,] - mean_Lambda[idx1,]
# set2 <- mean_Eta[idx2,] - mean_Lambda[idx2,]
# plot(set1, set2)
# cor(set1, set2)
# cat("Pair:",idx1,"x",idx2,"\n")

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

# (0) Fix issue with double NA's in taxonomy (after re-closing)
# (1) Presence-absence across subjects seems to drive correlation; scale to lots of subjects
# (2) Prior choices OK?
# (3) Fit across many subjects and build "rug" plot



