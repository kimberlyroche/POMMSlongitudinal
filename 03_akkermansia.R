# This code identifies individuals with high* Akkermansia.

library(tidyverse)

source("00_functions.R")

# ------------------------------------------------------------------------------
#   Parse data
# ------------------------------------------------------------------------------

data_obj <- readRDS(file.path("data", "processed_data.rds"))
counts <- data_obj$counts
tax <- data_obj$tax
metadata <- data_obj$metadata
rm(data_obj)

subjects <- unique(metadata$ind_id)
akkermansia_idx <- which(tax[,6] == "Akkermansia")

# Filter to subjects with minimum Akkermansia relative abundance >= 0.1%

threshold <- 0.1/100
keep_subject_mean <- logical(length(subjects))
keep_subject_min <- logical(length(subjects))
akkermansia_range <- data.frame(value = c(), type = c(), subject = c())
for(s in 1:length(subjects)) {
  subject <- subjects[s]
  subj_counts <- counts[, which(metadata$ind_id == subject), drop = F]
  subj_props <- apply(subj_counts, 2, function(x) x/sum(x))
  if(mean(subj_props[akkermansia_idx,]) >= threshold) {
    keep_subject_mean[s] <- TRUE
  }
  if(min(subj_props[akkermansia_idx,]) >= threshold) {
    keep_subject_min[s] <- TRUE
  }
  akkermansia_range <- rbind(akkermansia_range,
                             data.frame(value = c(min(subj_props[akkermansia_idx,]),
                                                  mean(subj_props[akkermansia_idx,]),
                                                  max(subj_props[akkermansia_idx,])),
                                        type = c("min", "mean", "max"),
                                        subject = subject))
}

# Show super high scores
# akkermansia_range[akkermansia_range$type == "max" & akkermansia_range$value > 0.2,]

cat(paste0("Subject meeting MEAN threshold: ",
           sum(keep_subject_mean),
           " / ",
           length(keep_subject_mean),
           " (",
           round(sum(keep_subject_mean)/length(keep_subject_mean),2),
           ")\n"))

cat(paste0("Subject meeting MIN threshold: ",
           sum(keep_subject_min),
           " / ",
           length(keep_subject_min),
           " (",
           round(sum(keep_subject_min)/length(keep_subject_min),2),
           ")\n"))

# Show range of Akkermansia relative abundances across subjects

akkermansia_range$idx <- rep(1:length(subjects), 3)
ggplot(akkermansia_range, aes(x = idx,
                              y = value,
                              color = type)) +
  geom_point() +
  geom_segment(aes(x = 1, xend = length(subjects),
                   y = threshold, yend = threshold), color = "black") +
  facet_grid(. ~ type) +
  xlab("subject index") +
  ylab("relative abundance (Akkermansia)")
ggsave(file.path("output", "images", "Akkermansia_thresholds.png"),
       units = "in",
       dpi = 100,
       height = 5,
       width = 10)

# Alternatively, filter to subjects with top 10 mean Akk. relative abundance
top10 <- akkermansia_range %>%
  filter(type == "mean") %>%
  arrange(desc(value)) %>%
  slice(1:10) %>%
  select(subject)

# Threshold 1
# saveRDS(as.numeric(subjects[keep_subject_mean]), file = file.path("data", "akkermansia_subjects.rds"))

# Threshold 2
saveRDS(as.numeric(top10$subject), file = file.path("data", "akkermansia_subjects.rds"))

