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
  if(taxon_idx == nrow(tax_map)) {
    "assorted low abundance"
  } else if(taxon_idx > nrow(tax_map)) {
    "N/A"
  } else {
    tax_pieces <- tax_map[taxon_idx,]
    level_idx <- max(which(!is.na(tax_pieces)))
    paste0(names(tax_pieces[level_idx]), " ", as.character(tax_pieces[[level_idx]]))
  }
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
Gamma <- diag(nrow(X))

# Fit model and convert to CLR
cat("Fitting model...\n")
fit <- pibble(as.matrix(Y), X, upsilon, Theta, Gamma, Xi, n_samples = 500) # takes about 5 sec. at 116 taxa x 15 samples
cat("Model fit!\n")
fit.clr <- to_clr(fit)

saveRDS(list(fit = fit, fit.clr = fit.clr), file = "fitted_model.rds")
quit()

# -------------------------------------------------------------------------------------------------
#   RESULTS: TAKE A LOOK AT STRONG CORRELATORS
# -------------------------------------------------------------------------------------------------

MAP_estimate <- TRUE
if(MAP_estimate) {
  fit_obj <- readRDS("fitted_model_2.rds")
} else {
  # Note: this output is erroneous; wrong input data!
  fit_obj <- readRDS("fitted_model.rds")
}
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

temp <- mean_Sigma
# temp[upper.tri(temp, diag = TRUE)] <- NA
temp <- as.data.frame(cbind(1:nrow(temp), temp))
colnames(temp) <- c("taxon1", 1:nrow(temp))
temp <- pivot_longer(temp, !taxon1, names_to = "taxon2", values_to = "correlation")
temp$taxon2 <- as.numeric(temp$taxon2)
# temp <- temp[complete.cases(temp),]

ggplot(temp, aes(x = taxon1, y = taxon2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", high = "darkred")
ggsave("Sigma.png", units = "in", dpi = 100, height = 6, width = 8)

# -------------------------------------------------------------------------------------------------
#   FILTER TO NON-ZERO-SPANNING POSTERIOR INTERVALS
# -------------------------------------------------------------------------------------------------

filter_CIs <- function(Sigma, correlators, threshold) {
  df <- data.frame(tag1 = c(), tag2 = c(), pair_name = c(), left = c(), middle = c(), right = c())
  for(idx in correlators) {
    pair <- combos[,idx]
    x <- Sigma[pair[1],pair[2],]
    # Filter to a large mean
    mu <- mean(x)
    if((sign(threshold) == 1 & mu > threshold) | (sign(threshold) == -1 & mu < threshold)) {
      pair_name <- paste0(get_tax_label(pair[1], tax), " x ", get_tax_label(pair[2], tax))
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

combos <- combn(1:nrow(mean_Sigma), m = 2)

# Threshold to correlators with non-zero posterior 95% CIs
non_zero <- sapply(1:ncol(combos), function(x) {
  pair <- combos[,x]
  sum(sign(quantile(Sigma_corr[pair[1],pair[2],], probs = c(0.025, 0.975))))
})

negative_correlators <- which(non_zero < 0)
positive_correlators <- which(non_zero > 0)

df_neg <- filter_CIs(Sigma_corr, negative_correlators, threshold = -0.5)
df_pos <- filter_CIs(Sigma_corr, positive_correlators, threshold = 0.5)

if(nrow(df_neg) > 0) {
  write.table(df_neg, file = "negative_correlators.tsv", sep = "\t")
}

if(nrow(df_pos) > 0) {
  write.table(df_pos, file = "positive_correlators.tsv", sep = "\t")
}

for(i in 1:nrow(df_pos)) {
  x <- as.vector(fit.clr$Eta[df_pos[i,]$tag1,,1] - fit.clr$Lambda[df_pos[i,]$tag1,,1]%*%fit.clr$X)
  y <- as.vector(fit.clr$Eta[df_pos[i,]$tag2,,1] - fit.clr$Lambda[df_pos[i,]$tag2,,1]%*%fit.clr$X)
  tax1 <- get_tax_label(df_pos[i,]$tag1, tax)
  tax2 <- get_tax_label(df_pos[i,]$tag2, tax)
  p <- ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(size = 2) +
    theme_bw() +
    xlab(paste0("CLR(",tax1,")")) +
    ylab(paste0("CLR(",tax2,")"))
  show(p)
  ggsave(paste0("check_",i,".png"), p, units = "in", dpi = 100, height = 4, width = 6)
}

# Bar plots -- need intervals for these to be decent looking
# df <- rbind(cbind(df_neg, type = "negative"), cbind(df_pos, type = "positive"))
# df$pair_name <- as.factor(df$pair_name)
# df$type <- as.factor(df$type)
# 
# palette2 <- c("#3d32a1", "#c23a25")
# p <- ggplot(df[!sapply(as.character(df$pair_name), function(x) { str_detect(x, "N/A") }),],
#             aes(x = reorder(pair_name, middle))) +
#   geom_boxplot(aes(ymin = left, lower = left, middle = middle, upper = right, ymax = right, fill = type, color = type),
#                 stat = "identity") +
#   # coord_flip() +
#   scale_fill_manual(values = palette2) +
#   scale_color_manual(values = palette2) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
#   xlab("bacterial pair") +
#   ylab("correlation")
# show(p)
# ggsave("strong_correlators.png", p, units = "in", dpi = 100, height = 10, width = 14)

