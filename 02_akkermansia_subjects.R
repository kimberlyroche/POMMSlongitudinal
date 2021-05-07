# This code identifies individuals with high Akkermansia.

source("00_functions.R")

# ------------------------------------------------------------------------------
#   Global vars
# ------------------------------------------------------------------------------

filter_method <- 2
# 1: Filter to subjects with mean Akkermansia relative abundance > 0.1%.
# 2: Filter to 10 subjects with highest Akkermansia relative abundance.

plot_subject_rel_ab <- FALSE

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

# Pull min, mean, and max Akkermansia relative abundance for each subject.
akkermansia_range <- data.frame(value = c(), type = c(), subject = c())
for(s in 1:length(subjects)) {
  subject <- subjects[s]
  subj_counts <- counts[, which(metadata$ind_id == subject), drop = F]
  subj_props <- apply(subj_counts, 2, function(x) x/sum(x))
  akkermansia_range <- rbind(akkermansia_range,
                             data.frame(value = c(min(subj_props[akkermansia_idx,]),
                                                  mean(subj_props[akkermansia_idx,]),
                                                  max(subj_props[akkermansia_idx,])),
                                        type = c("min", "mean", "max"),
                                        subject = subject))
}

# Plot.
if(plot_subject_rel_ab) {
  akkermansia_range$idx <- rep(1:length(subjects), 3)
  p <- ggplot(akkermansia_range, aes(x = idx,
                                y = value,
                                color = type)) +
    geom_point() +
    geom_segment(aes(x = 1, xend = length(subjects),
                     y = 0.1/100, yend = 0.1/100), color = "black") +
    facet_grid(. ~ type) +
    xlab("subject index") +
    ylab("relative abundance (Akkermansia)")
  ggsave(file.path("output", "images", "Akkermansia_thresholds.png"),
         p,
         units = "in",
         dpi = 100,
         height = 5,
         width = 10)
}

if(filter_method == 1) {
  threshold <- 0.1/100
  keep_subjects <- as.numeric(akkermansia_range %>%
                              filter(type == "mean") %>%
                              filter(value > threshold) %>%
                              pull(subject))
  saveRDS(keep_subjects, file = file.path("data", "akkermansia_subjects.rds"))
} else if(filter_method == 2) {
  keep_subjects <- as.numeric(akkermansia_range %>%
                              filter(type == "mean") %>%
                              arrange(desc(value)) %>%
                              slice(1:10) %>%
                              pull(subject))
  saveRDS(keep_subjects, file = file.path("data", "akkermansia_subjects.rds"))
}  
