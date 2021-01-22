library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(driver)
library(fido)

setwd("C:/Users/kim/Documents/POMMSlongitudinal/")

# Pull per-group sample identifiers
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

# Note: combined_tax is assumed global
get_tax_label <- function(taxon_idx) {
  tax_pieces <- combined_tax[taxon_idx,]
  level_idx <- max(which(!is.na(tax_pieces)))
  paste0(names(tax_pieces[level_idx]), " ", tax_pieces[level_idx])
}

plot_sampled_bars <- function(A, B = NULL, palette = NULL) {
  if(is.null(B)) {
    combined_counts <- A
  } else {
    sample_sz <- 10
  # Subsample for combined counts (otherwise the plot is visually awful)
    combined_counts <- rbind(A[sample(1:nrow(A), size = sample_sz),],
                             B[sample(1:nrow(B), size = sample_sz),])
  }
  combined_props <- t(apply(combined_counts, 1, function(x) x/sum(x)))

  tax_labels <- unname(sapply(1:(ncol(A)-1), get_tax_label))
  # Some of these are liable not to be unique; fluff up the strings to they are
  redundant_labels <- names(which(table(tax_labels) > 1))
  counters <- numeric(length(redundant_labels))
  names(counters) <- redundant_labels
  for(i in 1:length(tax_labels)) {
    if(tax_labels[i] %in% names(counters)) {
      counters[[tax_labels[i]]] <- counters[[tax_labels[i]]] + 1
      tax_labels[i] <- paste0("ASV #",counters[[tax_labels[i]]]," in ",tax_labels[i])
    } else {
      tax_labels[i] <- paste0("ASV #1 in ",tax_labels[i])
    }
  }

  colnames(combined_props) <- c(tax_labels, "assorted low abundance")
  combined_props <- cbind(sample_index = 1:nrow(combined_props), combined_props)
  plot_data <- pivot_longer(as.data.frame(combined_props), !sample_index, names_to = "taxon", values_to = "relative_abundance")
  plot_data$taxon <- as.factor(plot_data$taxon)

  if(is.null(palette)) {
    palette <- generate_highcontrast_palette(ncol(combined_props))
  }
  p <- ggplot(plot_data, aes(fill = taxon, y = relative_abundance, x = sample_index)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette)
  if(length(tax_labels) > 20) {
    p <- p + theme(legend.position = "none")
  }
  p
}

data <- readRDS("data/phyloseq_r24_complete_16S_metadata_corrected.rds")
metadata <- sample_data(data)
counts <- otu_table(data)@.Data # 960 samples x ~151K taxa
tax <- tax_table(data)@.Data
rownames(tax) <- NULL # kill the huge sequence labels

# First filter to cohorts we care about: OB and HWC

retain_samples <- metadata[metadata$group %in% c("OB", "HWC"),]$sample_id
counts <- counts[rownames(counts) %in% retain_samples,]
metadata <- metadata[metadata$sample_id %in% retain_samples,]

cat("Starting with",nrow(counts),"samples x",ncol(counts),"taxa\n")

# Collapse taxa below a minimum relative abundance to an "other" category
props <- t(apply(counts, 1, function(x) x/sum(x)))
max_relative_abundance <- apply(props, 2, max)
plot(max_relative_abundance[1:200])
# retain_taxa <- max_relative_abundance >= 0.01
retain_taxa <- max_relative_abundance >= 0.05

cat("Retaining",round((sum(retain_taxa) / length(max_relative_abundance))*100, 1),"percent of taxa (",sum(retain_taxa),")\n")

# What proportion of total counts do the retained taxa represent?
retained_total <- sum(counts[,retain_taxa])
total <- sum(counts)

cat("Retained taxa account for",round((retained_total / total)*100, 1),"percent of total counts\n")

# Collapse counts
counts_collapsed <- cbind(counts[,retain_taxa], rowSums(counts[,!retain_taxa]))
counts <- counts_collapsed
cat("Collapsed to",nrow(counts_collapsed),"samples x",ncol(counts_collapsed),"taxa\n")
tax <- tax[retain_taxa,]
tax <- rbind(tax, NA)

# How many ASVs are defined at each taxonomic level?
tax_defined_count <- sapply(1:7, function(x) sum(!is.na(tax[,x])))
names(tax_defined_count) <- colnames(tax)
tax_defined_count

# We probably want to collapse to the genus level (TBD)

# Note: Total abundances are broadly comparable but there's a definite trend
# associated with extraction batch. (TBD)

# Order by batch
batch_labels <- metadata$dna_extraction_batch
batch_order <- order(batch_labels)

total_abundances <- apply(counts, 1, sum)
data <- data.frame(sample_index = 1:nrow(counts), total = total_abundances[batch_order], group = as.factor(batch_labels[batch_order]))
ggplot(data, aes(x = sample_index, y = total, color = group)) +
  geom_point()
ggsave("batch_effect.png", units = "in", dpi = 100, height = 6, width = 8)

# Sampling occurred every 6 weeks for 6 months giving at most 5 samples
# per subject
# 1         2         3         4         5
#  0 mo      1.5 mo    3 mo      4.5 mo    6 mo
#
# There are so few samples it might be worthwhile just to use pibble here
# That would allow us to (attempt to) model out differences on the basis of
#  sex, etc.

# Note: The HWC group only has a baseline sample
# table(metadata[metadata$group == "HWC",]$visit)

data_OB <- get_cohort_samples(cohort = "OB")
data_HWC <- get_cohort_samples(cohort = "HWC")

# Further strip down taxa to a visualizable subset of taxa
# This will have all taxa below a pretty large relative abundance collapsed into other
# We're shooting for ~20-50 ASVs

viz_counts_OB <- data_OB$counts
n_OB <- nrow(viz_counts_OB)

viz_counts_HWC <- data_HWC$counts
n_HWC <- nrow(viz_counts_HWC)

combined_counts <- rbind(viz_counts_OB, viz_counts_HWC)

props <- t(apply(combined_counts, 1, function(x) x/sum(x)))
mean_relative_abundance <- apply(props, 2, mean)
retain_taxa <- mean_relative_abundance >= 0.001

cat("Retaining",sum(retain_taxa),"taxa for visualization\n")

retain_taxa[length(retain_taxa)] <- FALSE # keep other in other

combined_counts_collapsed <- cbind(combined_counts[,retain_taxa], rowSums(combined_counts[,!retain_taxa]))
combined_counts <- combined_counts_collapsed
cat("Collapsed to",nrow(combined_counts),"samples x",ncol(combined_counts),"taxa\n")
combined_tax <- tax[retain_taxa,]
combined_tax <- rbind(combined_tax, NA)

# Pull apart
viz_counts_OB <- combined_counts[1:n_OB,]
viz_counts_HWC <- combined_counts[(n_OB + 1):nrow(combined_counts),]

# Quick visualization of variation

# (1) HWC baseline x OBS baseline

palette <- generate_highcontrast_palette(ncol(viz_counts_OB))

A <- viz_counts_HWC
B <- viz_counts_OB[data_OB$visits == 1,]
plot_sampled_bars(A, B, palette = palette)

# (2) OBS baseline x OBS final time point

A <- viz_counts_OB[data_OB$visits == 1,]
B <- viz_counts_OB[data_OB$visits == 5,]
plot_sampled_bars(A, B, palette = palette)

# (3) OBS series for a few subjects

subject <- sample(names(which(table(data_OB$subjects) == 5)), size = 1)
A <- viz_counts_OB[data_OB$subject == subject,]
plot_sampled_bars(A, B = NULL, palette = palette)

# PCA w/ samples
# OB individuals seem to have a lot more subject-specific stuff -- true?

# Compute Aitchison distance on 6 categories of samples:
# HWC (baseline)
# OB (observations 1-5)

combined_counts <- rbind(data_OB$counts[data_OB$visits == 1,], data_HWC$counts)

# Label by cohort x time point
label_cohort <- c(rep("HWC", n_HWC), rep("OB", n_OB))
# label_time <- c(data_HWC$visits, data_OB$visits)
# label <- as.factor(paste0(label_cohort, label_time))
label <- c(rep("HWC", n_HWC), rep("OB", sum(data_OB$visits == 1)))

# Label by subject
# label <- as.factor(c(data_HWC$subjects, data_OB$subjects))

clr_counts <- clr_array(combined_counts + 0.5, parts = 2)
d <- dist(clr_counts)
coords <- cmdscale(d, k = 2)
data <- data.frame(x = coords[,1], y = coords[,2], cohort = label)

ggplot(data, aes(x = x, y = y, color = label)) +
  geom_point(size = 2) +
  scale_color_manual(values = generate_highcontrast_palette(length(unique(label)))) #+
  # theme(legend.position = "none")

# Interesting questions
# - What's different between extremes of PC1, PC2?
# - Pull other patient metadata
# - Start pibble model for these

# Start combining two subjects
Y <- NULL
subjects <- c(11,124)
subject_labels <- c()
for(subject in subjects) {
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
subject_labels <- as.factor(subject_labels)
X <- t(model.matrix(~subject_labels))

upsilon <- nrow(Y) + 3
Omega <- diag(nrow(Y))
G <- cbind(diag(nrow(Y)-1), -1)
Xi <- (upsilon-nrow(Y))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y)-1, nrow(X))
Gamma <- diag(nrow(X))

fit <- pibble(Y, X, upsilon, Theta, Gamma, Xi)
fit.clr <- to_clr(fit)

A <- cov2cor(fit.clr$Sigma[,,sample(1:fit$iter, size = 1)])
B <- cov2cor(fit.clr$Sigma[,,sample(1:fit$iter, size = 1)])

par(mfrow = c(1, 2))
image(A)
image(B)

dim(A)
A[1:3,1:3]
