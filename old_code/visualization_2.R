# Note: This model is kind of weird.
#   ~ subject_labels + bmi_covariate1 + bmi_covariate2 + ir_covariate
# We're trying to associate average logratio abundance with total-change covariates.
# Not sure if that's interesting/relevant.

# We're assuming `get_tax_label` is in the workspace (loaded from 02_joint_model.R)

data <- readRDS("fitted_model_MAP2.rds")

Lambda <- data$fit.clr$Lambda[,,1]
Eta <- data$fit.clr$Eta[,,1]
X <- data$fit.clr$X

input_data_obj <- readRDS("processed_data.rds")
tax <- input_data_obj$tax

params <- list(colname = "bmi_covariate1",
               label = "baseline BMI (scaled)",
               coef_pos_threshold = 1.5,
               coef_neg_threshold = -0.5)

params <- list(colname = "bmi_covariate2",
               label = "change in BMI (scaled)",
               coef_pos_threshold = 0.5,
               coef_neg_threshold = -0.5)

params <- list(colname = "ir_covariate",
               label = "change in HOMA IR (scaled)",
               coef_pos_threshold = 0.75,
               coef_neg_threshold = -0.4)

cov_idx <- which(rownames(data$fit.clr$X) == params$colname)
linear_predictor <- Lambda[,cov_idx,drop=F]%*%X[cov_idx,,drop=F]

coef_cov <- Lambda[,cov_idx]
coef_order <- order(coef_cov)
plot_data <- data.frame(x = 1:length(coef_cov), y = coef_cov[coef_order])
ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(size = 2) +
  xlab("taxon index (sorted)") +
  ylab(params$label)
ggsave("output/coefficients.png", units = "in", dpi = 100, height = 5, width = 6)

tax_idx <- sample(which(coef_cov > params$coef_pos_threshold), size = 2, replace = FALSE)
for(i in 1:2) {
  plot_data <- data.frame(x = X[cov_idx,], y = Eta[tax_idx[i],])
  ggplot(plot_data, aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab(params$label) +
    ylab(paste0("CLR(", get_tax_label(tax_idx[i], tax), ")"))
  ggsave(paste0("output/pos_Lambda_0",i,".png"), units = "in", dpi = 100, height = 5, width = 6)
}

tax_idx <- sample(which(coef_cov < params$coef_neg_threshold), size = 2, replace = FALSE)
for(i in 1:2) {
  plot_data <- data.frame(x = X[cov_idx,], y = Eta[tax_idx[i],])
  ggplot(plot_data, aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab(params$label) +
    ylab(paste0("CLR(", get_tax_label(tax_idx[i], tax), ")"))
  ggsave(paste0("output/neg_Lambda_0",i,".png"), units = "in", dpi = 100, height = 5, width = 6)
}

