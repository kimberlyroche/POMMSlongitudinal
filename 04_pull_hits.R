# This script parses model output to find correlated pairs of (CLR) taxa.

library(stringr)

source("00_functions.R")

# ------------------------------------------------------------------------------
#   Global vars
# ------------------------------------------------------------------------------

# See INSTRUCTIONS.txt for a list of possible scenarios here.

# Filter to only Akkermansia-related correlations.
filter_akkermansia_results <- FALSE

# Plot correlated pairs for (eyeball) validation.
plot_results <- TRUE

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

# Test to see if an association is outside the 90% interval of permutations.
test_association <- function(S1, S_arr, idx1, idx2, alpha = 0.05, verbose = FALSE) {
  obs <- S_arr[idx1,idx2,]
  bounds <- quantile(obs, probs = c(alpha/2, 1-alpha/2))
  result <- S1[idx1,idx2] < bounds[[1]] | S1[idx1,idx2] > bounds[[2]]
  if(verbose) {
    return(list(result = result, obs_corr = S1[idx1,idx2], possible_corr = obs))
  } else {
    return(list(result = result))
  }
}

# Plot the residual of the denoised relative abundances (i.e. the covariation
# theoretically not driven by the covariates).
plot_pair <- function(model_fit, idx1, idx2, label1, label2) {
  Eta <- model_fit$fit.clr$Eta[,,1]
  Lambda <- model_fit$fit.clr$Lambda[,,1]
  x <- as.vector(Eta[idx1,,drop=F] - Lambda[idx1,,drop=F]%*%model_fit$fit.clr$X)
  y <- as.vector(Eta[idx2,,drop=F] - Lambda[idx2,,drop=F]%*%model_fit$fit.clr$X)
  tax1 <- get_tax_label(idx1, tax)
  tax2 <- get_tax_label(idx2, tax)
  p <- ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(size = 2) +
    xlab(paste0("CLR(",tax1,")")) +
    ylab(paste0("CLR(",tax2,")")) +
    theme_bw() +
    theme(legend.position = "none")
  show(p)
  ggsave(file.path("output",
                   "images",
                   paste0("pair_",
                          gsub("[[:space:]]", "", tax1),
                          "_",
                          gsub("[[:space:]]", "", tax2),
                          ".png")),
         p,
         units = "in",
         dpi = 100,
         height = 4,
         width = 6)
}

# ------------------------------------------------------------------------------
#   Pull "hits"
# ------------------------------------------------------------------------------

data_obj <- readRDS(file.path("data", "processed_data.rds"))
tax <- data_obj$tax

canonical_fit <- readRDS(file.path("output", "fitted_models", "model_canonical.rds"))

# Pull all permuted output
permutation_files <- list.files(path = file.path("output", "fitted_models"),
                                pattern = "model_\\w{8}-\\w{4}-\\w{4}-\\w{4}-\\w{12}\\.rds",
                                full.names = TRUE)
canonical_Sigma <- canonical_fit$fit.clr$Sigma[,,1]
canonical_Sigma <- cov2cor(canonical_Sigma)
permuted_Sigmas <- array(NA, dim = c(nrow(canonical_Sigma),
                                     ncol(canonical_Sigma),
                                     length(permutation_files)))
for(f in 1:length(permutation_files)) {
  fit_obj <- readRDS(permutation_files[f])
  permuted_Sigmas[,,f] <- cov2cor(fit_obj$fit.clr$Sigma[,,1])
}

# Find "hits"
pairs <- combn(1:nrow(canonical_Sigma), m = 2)
if(filter_akkermansia_results) {
  akkermansia_idx <- which(tax[,6] == "Akkermansia")
  ak_cols <- which(pairs[1,] == akkermansia_idx | pairs[2,] == akkermansia_idx)
  
  results <- t(sapply(ak_cols, function(c_idx) {
    res <- test_association(canonical_Sigma, permuted_Sigmas, pairs[1,c_idx], pairs[2,c_idx], verbose = TRUE)
    res[1:2]
  }))
  results <- data.frame(result = unlist(results[,1]), corr = unlist(results[,2]))
  results$tax1_idx <- pairs[1,ak_cols]
  results$tax2_idx <- pairs[2,ak_cols]
  results$tax1_label <- sapply(results$tax1_idx, function(x) get_tax_label(x, tax))
  results$tax2_label <- sapply(results$tax2_idx, function(x) get_tax_label(x, tax))
  results <- results %>%
    select(!result)
} else {
  results <- t(sapply(1:ncol(pairs), function(c_idx) {
    res <- test_association(canonical_Sigma, permuted_Sigmas, pairs[1,c_idx], pairs[2,c_idx], verbose = TRUE)
    res[1:2]
  }))
  results <- data.frame(result = unlist(results[,1]), corr = unlist(results[,2]))
  results$tax1_idx <- pairs[1,]
  results$tax2_idx <- pairs[2,]
  results$tax1_label <- sapply(results$tax1_idx, function(x) get_tax_label(x, tax))
  results$tax2_label <- sapply(results$tax2_idx, function(x) get_tax_label(x, tax))
  results <- results %>%
    select(!result)
}

# Plot top 10 results (in terms of magnitude of correlations)
plot_results <- results %>%
  arrange(desc(abs(corr))) %>%
  slice(1:10)

for(i in 1:nrow(plot_results)) {
  plot_pair(canonical_fit,
            plot_results$tax1_idx[i],
            plot_results$tax2_idx[i],
            plot_results$tax1_label[i],
            plot_results$tax2_label[i])
}

# Write out results.
simple_results <- results %>%
  select(tax1_label, tax2_label, corr)
colnames(simple_results) <- c("taxon1", "taxon2", "correlation")

write.table(simple_results, file.path("output", "fitted_models", "correlators.tsv"), row.names = FALSE, sep = "\t")
