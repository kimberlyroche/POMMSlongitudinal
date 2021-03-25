# These are common/utility functions.

library(tidyverse)
library(phyloseq)
library(RColorBrewer)

dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)
dir.create(file.path("output", "fitted_models"), showWarnings = FALSE)

generate_highcontrast_palette <- function(S) {
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  sample(getPalette(S))
}

# This assumes a taxonomy as a taxa x hierarchy data.frame. Columns are
#   orders in the hierarchy (e.g. Phylum, ..., Genus).
# The last row of the taxonomy is is assumed to be an all-NA "other" category.
# This function pulls a short, readable label for this taxon based on the finest
#   level of identification.
get_tax_label <- function(taxon_idx, tax_map) {
  tax_pieces <- tax_map[taxon_idx,]
  resolved_levels <- which(!is.na(tax_pieces))
  if(length(resolved_levels) > 0) {
    level_idx <- max(resolved_levels)
    paste0(names(tax_pieces[level_idx]), " ", as.character(tax_pieces[[level_idx]]))
  }
  else {
    paste0("unknown taxon")
  }
}

# `data` is presumed to be a taxa x samples data.frame.
# This function collapses all taxa below some relative abundance (if parameter
#   `threshold` is < 1) or below some absolute abundance (if `threshold` > 1).
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
