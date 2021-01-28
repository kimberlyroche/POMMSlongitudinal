library(fido)
library(tidyverse)

# -------------------------------------------------------------------------------------------------
#   Functions
# -------------------------------------------------------------------------------------------------

fit_model <- function(subject, counts, metadata) {
  subject_samples <- metadata[metadata$ind_id == subject,]$sample_id
  # Treat samples as independent for now
  
  Y <- counts[,colnames(counts) %in% subject_samples]
  
  # Build design matrix
  X <- matrix(1, 1, ncol(Y)) # intercept-only model
  
  # -------------------------------------------------------------------------------------------------
  #   What to do about taxa absent in this individual?
  #   Spike-in a random 1-count?
  # -------------------------------------------------------------------------------------------------
  
  # Set priors a la
  # https://jsilve24.github.io/fido/articles/introduction-to-fido.html
  upsilon <- nrow(Y) + 3
  Omega <- diag(nrow(Y))
  G <- cbind(diag(nrow(Y)-1), -1)
  Xi <- (upsilon-nrow(Y))*G%*%Omega%*%t(G)
  Theta <- matrix(0, nrow(Y)-1, nrow(X))
  Gamma <- diag(nrow(X))
  
  # Fit model and convert to CLR
  cat("Fitting model...\n")
  fit <- pibble(as.matrix(Y), X, upsilon, Theta, Gamma, Xi) # takes about 5 sec. at 116 taxa x 15 samples
  cat("Model fit!\n")
  fit.clr <- to_clr(fit)
  
  # Scale Sigma to correlation
  Sigma_corr <- fit.clr$Sigma
  for(i in 1:fit$iter) {
    Sigma_corr[,,i] <- cov2cor(Sigma_corr[,,i])
  }
  
  # Get mean posterior predictions for a few parameters
  mean_Sigma <- apply(Sigma_corr, c(1,2), mean)
  
  # -------------------------------------------------------------------------------------------------
  #   What to do about tremendous posterior uncertainty?
  # -------------------------------------------------------------------------------------------------
  
  return(list(mean_Sigma = mean_Sigma, sample_Sigma = Sigma_corr[,,1]))
}

# Vectorize unique elements of covariance matrix
# If TRUE, labels returns a column-wise lookup of pairs
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

# -------------------------------------------------------------------------------------------------
#   Parse data and fit model to the first 10 subjects with 5 visits
# -------------------------------------------------------------------------------------------------

data_obj <- readRDS("processed_data.rds")
counts <- data_obj$counts
tax <- data_obj$tax
metadata <- data_obj$metadata
rm(data_obj)

sample_ids <- colnames(counts)
subjects <- metadata[metadata$sample_id %in% sample_ids,]$ind_id

full_sample_subj <- which(table(subjects) == 5)

pairs <- NULL
rug_consensus <- NULL
rug_sampled <- NULL
for(s_idx in 1:10) {
  cat("Evaluating subject",s_idx,"\n")
  subject <- names(full_sample_subj)[s_idx]
  Sigmas <- fit_model(subject, counts, metadata)
  mean_Sigma <- Sigmas$mean_Sigma
  sample_Sigma <- Sigmas$sample_Sigma
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
  rug.reordered <- rug[row_ordering,]
  rug.reordered <- rug.reordered[,column_ordering]
  pairs.reordered <- pairs[column_ordering]
  
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

rug_obj <- visualize_rug(rug_sampled, row_ordering = rug_obj$row_ordering, column_ordering = rug_obj$column_ordering)
show(rug_obj$plot)







