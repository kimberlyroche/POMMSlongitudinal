# POMMSlongitudinal

[Explain what model does.]

## Files

This repo consist of analysis files numbered to indicate the order in which they should be run.

**00_functions.R** contains functions re-used through analyses.

**01_data_processing.R** contains code for parsing the count tables and taxonomy,
agglomerating sequence variants by shared taxonomic labels, and filtering out lowly abundant
taxa (into an 'Other' category).

**02_akkermansia_subjects.R** identifies individuals with high relative abundances of g. Akkermansia
and saves this list to a file.

**03_model_fitting.R** fits a `fido::pibble` model to the data, inferring correlations across
logratio bacterial abundances. The fitted model (including fits to permutations of the data
are saved to files.

**04_pull_hits.R** parses the fitted models to identify pairwise correlations between logratio
bacteria that are likely to be significant.

## Regarding permutation testing 

The model must be fit once to the data and then subsequently (many times) to
"permutations" of the data set (where taxa are shuffled across samples).

  In `03_model_fitting.R` set `permutation_test <- FALSE` and run the script
  once to generate the "canonical" output then set `permutation_test <- TRUE`
  and re-run the script as many times as you want samples from the null
  distribution to evaluate significance against. (These will be saved as
  specially named output files.)

(A) Baseline-only samples VS. (B) full series

  For (A), in `03_model_fitting.R` set
    baseline_only <- TRUE
    
  For (B), in `03_model_fitting.R` set
    baseline_only <- FALSE
    
(A) All subjects VS. (B) Akkermansia-rich subjects

  For (A), in `03_model_fitting.R` set
    filter_akkermansia_subjects <- FALSE
    
  For (B), in `03_model_fitting.R` set
    filter_akkermansia_subjects <- TRUE
    
(A) Akkermansia-only correlations VS. (B) all pairwise correlations between
                                          (CLR) taxa

  For (A), in `04_pull_hits.R` set
    filter_akkermansia_results <- TRUE
    
  For (B), in `04_pull_hits.R` set
    filter_akkermansia_results <- FALSE
