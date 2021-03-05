library(driver)
library(lme4)
library(tidyverse)

get_tax_label <- function(taxon_idx, tax_map) {
  tax_pieces <- tax_map[taxon_idx,]
  level_idx <- max(which(!is.na(tax_pieces)))
  paste0(names(tax_pieces[level_idx]), " ", as.character(tax_pieces[[level_idx]]))
}

input_data_obj <- readRDS("processed_data.rds")
counts <- input_data_obj$counts
metadata <- input_data_obj$metadata
tax <- input_data_obj$tax

# The code below was taken from 02_joint_model.R

subject_tallies <- table(metadata$ind_id)
subjects <- as.numeric(names(subject_tallies)[which(subject_tallies >= 4)])

# Pull BMI and HOMA data
met <- read.table("data/MET-S-nonCook_7-23-20pedsobesity_r24_phenotype_wide_20200723.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
delta_BMI <- sapply(subjects, function(x) unname(met[met$IND == x,]$Delta_BMI_P95))
names(delta_BMI) <- subjects
BMI_P95 <- sapply(subjects, function(x) unname(met[met$IND == x,]$BMI_P95))
names(BMI_P95) <- subjects
PC_HOMA_IR <- sapply(subjects, function(x) unname(met[met$IND == x,]$PC_HOMA_IR))
names(PC_HOMA_IR) <- subjects

selected_subjects <- subjects[which(!is.na(delta_BMI) & !is.na(BMI_P95) & !is.na(PC_HOMA_IR))]

# For now, omit subjects with missing covariates
# I'm not sure how we would otherwise impute them (etc.)

# Center and scale these parameters
# This bungles interpretation but I don't think we want to directly interpret the regression coefs

delta_BMI <- delta_BMI[names(delta_BMI) %in% selected_subjects]
saved_names <- names(delta_BMI)
delta_BMI <- as.vector(scale(delta_BMI, center = FALSE, scale = TRUE))
names(delta_BMI) <- saved_names
BMI_P95 <- BMI_P95[names(BMI_P95) %in% selected_subjects]
saved_names <- names(BMI_P95)
BMI_P95 <- as.vector(scale(BMI_P95, center = FALSE, scale = TRUE))
names(BMI_P95) <- saved_names
PC_HOMA_IR <- PC_HOMA_IR[names(PC_HOMA_IR) %in% selected_subjects]
saved_names <- names(PC_HOMA_IR)
PC_HOMA_IR <- as.vector(scale(PC_HOMA_IR, center = FALSE, scale = TRUE))
names(PC_HOMA_IR) <- saved_names

# Resample counts; wildly inefficient
# props <- apply(counts, 2, function(x) x/sum(x))
# resampled_counts <- array(NA, dim = c(nrow(counts), ncol(counts), 100))
# resampled_clr <- array(NA, dim = c(nrow(counts), ncol(counts), 100))
# for(i in 1:100) {
#   for(j in 1:ncol(counts)) {
#     x <- counts[,j]
#     y <- props[,j]
#     resampled_counts[,j,i] <- rmultinom(1, size = sum(x), prob = y)
#     resampled_clr[,j,i] <- clr(resampled_counts[,j,i] + 0.1)
#   }
# }

tax_idx <- 1

# # Summarize the relevant individuals' series.
# subject_labels <- c()
# Y <- NULL
# for(subject in selected_subjects) {
#   col_selector <- colnames(clr.counts) %in% metadata[metadata$ind_id == subject,]$sample_id
#   subject_logratios <- resampled_clr[,col_selector,]
#   final_timepoint <- subject_logratios[,dim(subject_logratios)[2],]
#   initial_timepoint <- subject_logratios[,1,]
#   delta_logratios <- final_timepoint - initial_timepoint
#   tax_delta <- delta_logratios[tax_idx,]
#   if(is.null(Y)) {
#     Y <- tax_delta
#   } else {
#     Y <- cbind(Y, tax_delta)
#   }
#   subject_labels <- c(subject_labels, subject)
# }
# subject_labels <- as.factor(subject_labels)

clr.counts <- clr_array(as.matrix(counts) + 0.1, parts = 1)
subject_labels <- c()
Y <- NULL
for(subject in selected_subjects) {
  col_selector <- colnames(clr.counts) %in% metadata[metadata$ind_id == subject,]$sample_id
  subject_logratios <- clr.counts[,col_selector]
  delta_logratios <- subject_logratios[,ncol(subject_logratios)] - subject_logratios[,1]
  if(is.null(Y)) {
    Y <- delta_logratios
  } else {
    Y <- cbind(Y, delta_logratios)
  }
  subject_labels <- c(subject_labels, subject)
}
subject_labels <- as.factor(subject_labels)

# colnames(Y) <- names(delta_BMI)
# Y <- rbind(Y, delta_BMI)
# Y <- t(Y)
# colnames(Y) <- c(1:(ncol(Y)-1), "delta_BMI")
# lm_data <- pivot_longer(as.data.frame(Y), !delta_BMI, names_to = "taxon", values_to = "delta_logratio")

colnames(Y) <- names(PC_HOMA_IR)
Y <- rbind(Y, PC_HOMA_IR)
Y <- t(Y)
colnames(Y) <- c(1:(ncol(Y)-1), "PC_HOMA_IR")
lm_data <- pivot_longer(as.data.frame(Y), !PC_HOMA_IR, names_to = "taxon", values_to = "delta_logratio")

sig_tax <- c()
for(tax_idx in unique(lm_data$taxon)) {
  # fit <- lm(delta_logratio ~ delta_BMI, data = lm_data[lm_data$taxon == tax_idx,])
  fit <- lm(delta_logratio ~ PC_HOMA_IR, data = lm_data[lm_data$taxon == tax_idx,])
  pval <- coef(summary(fit))[2,4]
  if(pval < 0.05) {
    sig_tax <- c(sig_tax, tax_idx)
  }
}

plot_counter <- 1
for(tax_idx in sig_tax) {
  # fit <- lm(delta_logratio ~ delta_BMI, data = lm_data[lm_data$taxon == tax_idx,])
  fit <- lm(delta_logratio ~ PC_HOMA_IR, data = lm_data[lm_data$taxon == tax_idx,])
  # plot_data <- data.frame(x = lm_data[lm_data$taxon == tax_idx,]$delta_BMI,
  #                         y = lm_data[lm_data$taxon == tax_idx,]$delta_logratio)
  plot_data <- data.frame(x = lm_data[lm_data$taxon == tax_idx,]$PC_HOMA_IR,
                          y = lm_data[lm_data$taxon == tax_idx,]$delta_logratio)
  ggplot(plot_data, aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    # xlab("change in BMI P95 (scaled)") +
    xlab("change in HOMA IR (scaled)") +
    ylab(paste0("change in CLR(", get_tax_label(tax_idx, tax), ")"))
  ggsave(paste0("output/",plot_counter,".png"), units = "in", dpi = 100, height = 5, width = 6)
  plot_counter <- plot_counter + 1
}

# Look at overall changes in abundance of CLR features
Y <- cbind(c(sapply(1:(nrow(Y)-1), function(i) get_tax_label(i, tax)), "low abundance group"), Y)
rownames(Y) <- NULL
colnames(Y) <- c("taxon", subject_labels)
plot_data <- pivot_longer(as.data.frame(Y), !taxon, names_to = "subject", values_to = "delta_CLR")
plot_data$delta_CLR <- as.numeric(plot_data$delta_CLR)

plot_data <- plot_data %>%
  group_by(taxon) %>%
  mutate(mid = mean(delta_CLR)) %>%
  mutate(lower = mean(delta_CLR) - sd(delta_CLR)) %>%
  mutate(upper = mean(delta_CLR) + sd(delta_CLR)) %>%
  # mutate(qlower = quantile(delta_CLR, c(0.05))) %>%
  # mutate(qmid = quantile(delta_CLR, c(0.5))) %>%
  # mutate(qupper = quantile(delta_CLR, c(0.95))) %>%
  select(taxon, lower, mid, upper) %>%
  distinct(taxon, lower, mid, upper)

# palette2 <- c("#3d32a1", "#c23a25")
p <- ggplot(plot_data, aes(x = reorder(taxon, mid))) +
  geom_boxplot(aes(ymin = lower, lower = lower, middle = mid, upper = upper, ymax = upper),
               stat = "identity") +
  # scale_fill_manual(values = palette2) +
  # scale_color_manual(values = palette2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
  xlab("taxon") +
  ylab("change in CLR abundance")
show(p)
ggsave("delta_clr_boxplot.png", units = "in", dpi = 100, height = 6, width = 12)

# Means only (misleading)
# p <- ggplot(plot_data, aes(x = reorder(taxon, qmid))) +
#   geom_point(aes(y = qmid)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
#   xlab("taxon") +
#   ylab("change in CLR abundance")
# show(p)
# ggsave("delta_clr_means.png", units = "in", dpi = 100, height = 6, width = 12)

