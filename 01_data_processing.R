# This script filters the data.

source("00_functions.R")

# ------------------------------------------------------------------------------
#   Global vars
# ------------------------------------------------------------------------------

plot_sample_series <- FALSE

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

# Pull per-group counts and metadata.
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
  return(list(cohort = cohort,
              counts = cohort_counts,
              subjects = ind_ids,
              visits = visits))
}

# ------------------------------------------------------------------------------
#   Parse data
# ------------------------------------------------------------------------------

data <- readRDS("data/phyloseq_r24_complete_16S_metadata_corrected.rds")
metadata <- sample_data(data)
counts <- otu_table(data)@.Data # 960 samples x ~151K taxa
tax <- tax_table(data)@.Data
rownames(tax) <- NULL # kill the huge sequence labels for now

# ------------------------------------------------------------------------------
#   Filter
# ------------------------------------------------------------------------------

# Filter to the cohort we care about (OB) using get_cohort_samples() function.

retain_samples <- metadata[metadata$group == "OB",]$sample_id
counts <- counts[rownames(counts) %in% retain_samples,]
metadata <- metadata[metadata$sample_id %in% retain_samples,]
colnames(counts) <- NULL

cat("Starting with",nrow(counts),"samples x",ncol(counts),"taxa\n")

# Remove all-zero taxa (almost 1/4).
present_taxa <- which(colSums(counts) > 0)
counts <- counts[,present_taxa]
tax <- tax[present_taxa,]

# There's a phenomenon here where some intermediate taxonomic levels are missing
# Replace these interstitial <NA>s with "(missing)" so we can differentiate them
# (otherwise the agglomeration I'm doing below truncates erroneously on these
# and causes problems).

# This will transform
#   Bacteria Firmicutes Clostridia Peptostreptococcales-Tissierellales <NA>
#     Finegoldia
# into
#   Bacteria Firmicutes Clostridia Peptostreptococcales-Tissierellales (Missing)
#     Finegoldia

for(i in 1:nrow(tax)) {
  idx <- which(is.na(tax[i,1:6]))
  if(length(idx) > 0) {
    first_na_idx <- min(idx)
    if(sum(is.na(tax[i,first_na_idx:6])) != length(first_na_idx:6)) {
      # Replace interstitial <NA>s with a string
      label_encountered <- FALSE
      for(j in 6:1) {
        if(!is.na(tax[i,j])) {
          label_encountered <- TRUE
        } else {
          if(label_encountered) {
            tax[i,j] <- "(Missing)"
          }
        }
      }
    }
  } # else fully defined
}

# Collapse to highest common taxonomic level.
agglomerated_counts <- NULL
for(tax_level in 1:6) {
  cat("Agglomerating at tax level:",colnames(tax)[tax_level],"\n")
  counts_subset <- as.data.frame(cbind(tax, t(counts)), stringsAsFactors = FALSE)
  counts_subset <- counts_subset[which(!is.na(counts_subset[,tax_level]) &
                                         is.na(counts_subset[,tax_level+1])),]
  copy_long <- pivot_longer(counts_subset,
                            !c("domain",
                               "phylum",
                               "class",
                               "order",
                               "family",
                               "genus"),
                            names_to = "sample",
                            values_to = "count")
  if(nrow(copy_long) > 0) {
    copy_long$count <- as.numeric(copy_long$count)
    result <- copy_long %>%
      group_by(domain, phylum, class, order, family, genus, sample) %>%
      summarize(total_counts = sum(count), .groups = 'drop')
    counts_subset <- as.data.frame(pivot_wider(result,
                                               names_from = "sample",
                                               values_from = "total_counts"))
    if(is.null(agglomerated_counts)) {
      agglomerated_counts <- counts_subset
    } else {
      agglomerated_counts <- rbind(agglomerated_counts, counts_subset)
    }
  }
}

# The pivots will have disordered the samples; get them back in the original
# (subject x visit) order.
orig_sample_order <- rownames(counts)
new_tax <- agglomerated_counts[,1:6]
new_counts <- agglomerated_counts[,7:ncol(agglomerated_counts)]
new_counts <- new_counts[,orig_sample_order]

counts <- new_counts
tax <- new_tax
# Note: `counts` and `tax` now have taxa as rows

# Remove taxa that aren't present (non-zero) in at least one of every subject's
# samples. We'll put this into an "other" category (the last row in the
# taxonomy).
retain_taxa <- c()
n_subjects <- length(unique(metadata$ind_id))
for(tax_idx in 1:nrow(counts)) {
  abundance <- unlist(counts[tax_idx,])
  names(abundance) <- metadata$ind_id
  d <- data.frame(count = unname(abundance), subject = metadata$ind_id)
  subject_no_observation <- d %>%
    group_by(subject) %>%
    summarize(total_observations = sum(count), .groups = 'drop') %>%
    filter(total_observations == 0) %>%
    nrow()
  if(subject_no_observation < n_subjects / 2) {
    retain_taxa <- c(retain_taxa, TRUE)
  } else {
    retain_taxa <- c(retain_taxa, FALSE)
  }
}

# Roll the rare-across-subjects taxa into the "other" category (last row).
new_counts <- counts[retain_taxa,]
new_counts <- rbind(new_counts, rowSums(counts[!retain_taxa,]))
new_tax <- tax[retain_taxa,]
new_tax <- rbind(new_tax, NA)

counts <- new_counts
tax <- new_tax

# Count features with max abundance < 1%.
props <- apply(counts, 2, function(x) x/sum(x))
max_relative_abundance <- apply(props, 1, max)
retain_taxa <- max_relative_abundance >= 0.005
cat(paste0("Retaining ", round((sum(retain_taxa) /
                                  length(max_relative_abundance))*100, 1),
           " percent of taxa (",sum(retain_taxa),")\n"))

# Mark "other" as being not retained, so that we can bundle it and the rare guys
# together.
retain_taxa[length(retain_taxa)] <- FALSE
new_counts <- rbind(counts[retain_taxa,], rowSums(counts[!retain_taxa,]))
retain_taxa[length(retain_taxa)] <- TRUE
new_tax <- tax[retain_taxa,]

counts <- new_counts
tax <- new_tax

# Clean up.
rownames(counts) <- NULL
rownames(tax) <- NULL

# Report totals.
retained_total <- sum(counts[1:(nrow(counts)-1),])
total <- sum(counts)
cat("Retained taxa account for",
    round((retained_total / total)*100, 1),
    "percent of total counts\n")

# Check that all taxa are unique at this point:
# shortnames <- apply(tax, 1, function(x) paste(x, collapse = ""))
# which(table(shortnames) > 1)

# Visualize these compositions quickly.

if(plot_sample_series) {
  sample_no <- table(metadata$ind_id)
  select_subj <- sample(names(sample_no[sample_no >= 4]), size = 1)
  counts_subset <- counts[,colnames(counts) %in%
                            metadata[metadata$ind_id == select_subj,]$sample_id]
  props <- collapse_below_minimum(counts_subset, tax, threshold = 0.01)
  palette <- generate_highcontrast_palette(nrow(props$data))
  # Strip down to taxa.
  tax_labels <- unname(sapply(1:(nrow(props$data)-1),
                              function(x) {
                                get_tax_label(x, props$tax)
                              }))
  rownames(props$data) <- c(tax_labels, "assorted low abundance")
  props$data <- rbind(sample_index = 1:ncol(props$data), props$data)
  plot_data <- pivot_longer(as.data.frame(t(props$data)),
                            !sample_index,
                            names_to = "taxon",
                            values_to = "relative_abundance")
  plot_data$taxon <- as.factor(plot_data$taxon)
  p <- ggplot(plot_data, aes(fill = taxon,
                             y = relative_abundance,
                             x = sample_index)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette)
  ggsave(file.path("output", "images", paste0("subject_series.png")),
         p,
         units = "in",
         dpi = 100,
         height = 6,
         width = 9)
}

# Save filtered data.
saveRDS(list(counts = counts, tax = tax, metadata = metadata),
        file = file.path("data", "processed_data.rds"))
