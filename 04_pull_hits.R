
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

# Test to see if an association is outside the 90% interval of permutations
test_assoc <- function(S1, S_arr, idx1, idx2, alpha = 0.05, verbose = FALSE) {
  obs <- S_arr[idx1,idx2,]
  bounds <- quantile(obs, probs = c(alpha/2, 1-alpha/2))
  result <- S1[idx1,idx2] < bounds[[1]] | S1[idx1,idx2] > bounds[[2]]
  if(verbose) {
    return(list(result = result, obs_corr = S1[idx1,idx2], possible_corr = obs))
  } else {
    return(list(result = result))
  }
}


# ------------------------------------------------------------------------------
#   Find "hits" by permutation
# ------------------------------------------------------------------------------


canonical_Sigma <- fit_obj$fit.clr$Sigma[,,1]
canonical_Sigma <- cov2cor(canonical_Sigma)
if(permutation_test) {
  n_permute <- 100
  permuted_Sigmas <- array(NA, dim = c(nrow(tax), nrow(tax), n_permute))
  for(pp in 1:n_permute) {
    fit_obj <- fit_model(subjects_to_use, Y, X, metadata, permute = TRUE, MAP = TRUE)
    permuted_Sigmas[,,pp] <- cov2cor(fit_obj$fit.clr$Sigma[,,1])
  }
}






# Find "hits"
pairs <- combn(1:nrow(tax), m = 2)
results <- sapply(1:100, function(c_idx) {
  res <- test_assoc(canonical_Sigma, permuted_Sigmas, pairs[1,c_idx], pairs[2,c_idx], verbose = TRUE)
  if(res$result) {
    plot(density(res$possible_corr), main = paste0("Pair ",c_idx))
    lines(x = rep(res$obs_corr, 2), y = c(0, 100), col = "red")
  }
})

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
