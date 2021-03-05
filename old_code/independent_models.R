# This script fits the (fido::pibble) model.

library(fido)
library(driver)
library(tidyverse)

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

fit_model <- function(subject, counts, metadata) {
  subject_samples <- metadata[metadata$ind_id == subject,]$sample_id
  # Treat samples as independent for now
  
  Y <- counts[,colnames(counts) %in% subject_samples]
  Y <- as.matrix(Y)
  
  # Build design matrix
  X <- matrix(1, 1, ncol(Y)) # intercept-only model
  
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
  fit <- pibble(Y, X, upsilon, Theta, Gamma, Xi, ret_mean = TRUE, n_samples = 0)
  cat("Model fit!\n")
  fit.clr <- fido::to_clr(fit)
  
  # Scale Sigma to correlation
  Sigma_corr <- fit.clr$Sigma
  for(i in 1:fit.clr$iter) {
    Sigma_corr[,,i] <- cov2cor(Sigma_corr[,,i])
  }
  
  # Get mean posterior predictions for a few parameters
  mean_Sigma <- apply(Sigma_corr, c(1,2), mean)

  # Note: Optimization over eta constantly fails to converge. This isn't obvs.
  #   because of terrible choice of priors but I still haven't fixed it. Best
  #   just to return MAP estimates for now.
  
  return(list(fit = fit.clr,
              mean_Sigma = mean_Sigma,
              sample_Sigma = Sigma_corr[,,1]))
}

# Vectorize unique elements of the covariance matrix. If `labels` is TRUE, the
# function returns returns a column-wise lookup of taxon-taxon pairs.
vectorize_Sigma <- function(Sigma, labels = FALSE) {
  combos <- combn(nrow(Sigma), m = 2)
  vec <- sapply(1:ncol(combos), function(i) {
    Sigma[combos[1,i], combos[2,i]]
  })
  if(labels) {
    return(list(vec = vec, pairs = combos))
  } else {
    return(vec)
  }
}

# ------------------------------------------------------------------------------
#   Parse data and fit model
# ------------------------------------------------------------------------------

testing <- TRUE

data_obj <- readRDS("processed_data.rds")
counts <- data_obj$counts
tax <- data_obj$tax
metadata <- data_obj$metadata
rm(data_obj)

sample_ids <- colnames(counts)
subjects <- metadata[metadata$sample_id %in% sample_ids,]$ind_id

full_sample_subj <- which(table(subjects) > 4)

pairs <- NULL
rug_consensus <- NULL
rug_sampled <- NULL
fits <- list()
for(s_idx in 1:length(full_sample_subj)) {
# for(s_idx in 1:5) {
  cat("Evaluating subject",s_idx,"\n")
  subject <- names(full_sample_subj)[s_idx]
  res <- fit_model(subject, counts, metadata)
  fits[[s_idx]] <- res$fit
  mean_Sigma <- res$mean_Sigma
  sample_Sigma <- res$sample_Sigma
  if(is.null(pairs)) {
    res <- vectorize_Sigma(mean_Sigma, labels = TRUE)
    pairs <- res$pairs
    rug_consensus <- res$vec
    rug_sampled <- vectorize_Sigma(sample_Sigma)
  } else {
    rug_consensus <- rbind(rug_consensus, vectorize_Sigma(mean_Sigma))
    rug_sampled <- rbind(rug_sampled, vectorize_Sigma(sample_Sigma))
  }
}

# -------------------------------------------------------------------------------------------------
#   Visualize consensus rug
# -------------------------------------------------------------------------------------------------

visualize_rug <- function(rug, row_ordering = NULL, column_ordering = NULL) {
  if(is.null(row_ordering) | is.null(column_ordering)) {
    d <- dist(rug)
    clustering.hosts <- hclust(d)
    d <- dist(t(rug))
    clustering.interactions <- hclust(d)
    row_ordering <- clustering.hosts$order
    column_ordering <- clustering.interactions$order
  }
  # Reorder all
  # row_ordering <- 1:nrow(rug)
  # column_ordering <- 1:ncol(rug)
  
  rug.reordered <- rug[row_ordering,]
  rug.reordered <- rug.reordered[,column_ordering]
  
  rug.reordered <- as.data.frame(cbind(1:nrow(rug.reordered), rug.reordered))
  colnames(rug.reordered) <- c("host", 1:(ncol(rug.reordered)-1))
  
  df <- pivot_longer(rug.reordered, !host, names_to = "pair", values_to = "correlation")

  # Plot
  p <- ggplot(df, aes(pair, host)) +
    geom_tile(aes(fill = correlation)) +
    scale_fill_gradient2(low = "darkblue", high = "darkred")
  
  return(list(plot = p, row_ordering = row_ordering, column_ordering = column_ordering))
}

# -------------------------------------------------------------------------------------------------
#   Apply the consensus clustering to the random sample to see how much change there is
#   This is one way (?) to plainly visualize posterior uncertainty: disagreement in clustering
# -------------------------------------------------------------------------------------------------

rug_obj <- visualize_rug(rug_consensus, row_ordering = NULL, column_ordering = NULL)
show(rug_obj$plot)
ggsave("rug_v1.png", plot = rug_obj$plot, units = "in", dpi = 100, height = 8, width = 12)

rug_obj_sampled <- visualize_rug(rug_sampled, row_ordering = rug_obj$row_ordering, column_ordering = rug_obj$column_ordering)
show(rug_obj_sampled$plot)
ggsave("rug_v2.png", plot = rug_obj_sampled$plot, units = "in", dpi = 100, height = 8, width = 12)

