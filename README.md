# PDAC Prediction with Circulating miRNA

This repository contains the final project analysis pipeline for predicting pancreatic ductal adenocarcinoma (PDAC) from circulating miRNA measurements and demographics. The code is organized around one main modeling report and one downstream interpretation report.

## Main reports

### `PDACmiRNAAnalysis.rmd`

This is the primary modeling report. It:

- loads the GEO-derived normalized miRNA matrix and sample metadata
- applies the shared import-time `log2(x + 1)` transform
- creates the matched-set-aware training and holdout split
- runs blocked limma on the training set
- fits the six final model families
    - forward step logistic regression
    - lasso
    - elastic net
    - LDA
    - random forest
    - decision tree
- evaluates the final models on the holdout set
- writes the main summary tables and figures to `Tables/` and `Figures/`
- saves the downstream handoff object used by the second report

### `FinalAnalysis.rmd`

This is the downstream biology and annotation report. It does not rebuild the modeling pipeline. Instead, it reads the saved outputs from `PDACmiRNAAnalysis.rmd` and performs:

- downstream multiMiR queries for the ranked miRNA set
- validated target summaries
- disease/drug association summaries
- GO biological-process enrichment for validated target proteins
- merged summary tables that combine model-selection counts, limma statistics, and multiMiR results

## How the two reports interact

The handoff between the two reports is explicit.

`PDACmiRNAAnalysis.rmd` writes:

- `Data/Top-ranked consensus miRNA downstream inputs.rds`
    - a list containing:
        - `top_ranked_miRNA_details`
        - `merged_top_features`
- `Tables/Blocked limma filter results on full training set.csv`
- `Tables/Final Model Feature Counts.csv`

`FinalAnalysis.rmd` then reads those files and uses them to build the downstream query set. In practice:

1. The main report identifies the top-ranked non-tree consensus miRNA set.
2. It merges those with the top tree-model features into `merged_top_features`.
3. It saves both objects into `Data/Top-ranked consensus miRNA downstream inputs.rds`.
4. The downstream report loads that RDS, derives `multiMiR_query_id` values from the annotation field, and runs the multiMiR analysis on the merged query set.

`FinalAnalysis.rmd` also writes a reusable cache of the downstream query results to:

- `Data/Top-ranked consensus miRNA multiMiR results.rds`

That cache can be reused by setting `run_multimir_analysis <- FALSE` inside the downstream report.

## Recommended run order

Run the reports from the project root in this order:

1. `PDACmiRNAAnalysis.rmd`
2. `FinalAnalysis.rmd`

Example from an R session:

```r
rmarkdown::render("PDACmiRNAAnalysis.rmd")
rmarkdown::render("FinalAnalysis.rmd")
```

There is also a helper script for the main report:

```r
Rscript render_pdac_report.R
```

The main report is long-running on this dataset and can take about an hour to complete.

## Key directories

- `R/`: helper functions and model runners
- `Data/`: source data, cached limma results, and downstream RDS handoff files
- `Tables/`: exported CSV summaries
- `Figures/`: exported plots used by the rendered reports

## Main helper files

- `R/reportHelpers.R`: shared HTML-table and report-formatting helpers
- `R/limma_helpers.R`: shared blocked-limma ranking helpers
- `R/multi_mir_analysis.R`: reusable multiMiR query and summarization helpers
- `R/run_fold_*.R`: model-specific training helpers

## Package requirements

The main report depends on packages including:

- `readxl`
- `readr`
- `tidyr`
- `glmnet`
- `ggplot2`
- `pROC`
- `dplyr`
- `rpart`
- `randomForest`
- `limma`
- `mice`

The downstream report additionally relies on:

- `multiMiR`

Optional GO enrichment output in `FinalAnalysis.rmd` also uses:

- `clusterProfiler`
- `org.Hs.eg.db`

## Related report

`MultiMirAnalysis.rmd` is a related downstream report variant that uses the same `Data/Top-ranked consensus miRNA downstream inputs.rds` handoff pattern. The main end-to-end pair for this project is `PDACmiRNAAnalysis.rmd` followed by `FinalAnalysis.rmd`.

