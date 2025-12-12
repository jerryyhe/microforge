# microforge: Utilities for microbial genomics workflows

A compact R package with convenience utilities for qPCR standard curves, DESeq2 workflows, plate data tidying, and simple metagenomics processing. Designed to speed up common lab-data wrangling and small bioinformatics tasks.

## Highlights

- qPCR
  - Build standard curves and convert Cq values to copies: `calculate_std_curve()`, `calculate_copies()`.
- DESeq2 helpers
  - Load and format input, add labels, and run DESeq2 pipelines: `load_files()`, `files_to_array()`, `format_df()`, `add_labels()`, `run_deseq2()`.
- Plate utilities
  - Tidy 96- and 384-well plate exports: `tidy_plate96()`, `tidy_plate384()`.
- Metagenomics
  - Split taxonomy strings and compute relative abundance: `split_taxonomy()`, `calculate_relative_abundance()`.

See implementations in the `R/` directory (e.g., `R/qPCR_functions.R`, `R/deseq2_functions.R`, `R/plate_functions.R`, `R/metagenomics_functions.R`).

## Installation

From the repository root (recommended during development):

```r
# install devtools or remotes if needed
install.packages("devtools")
devtools::install_local(".")
# or load during development
devtools::load_all(".")
```

Or install directly from GitHub (replace `USERNAME`/`REPO` if published):

```r
remotes::install_github("USERNAME/microforge")
```

## Quick start

Load the package:

```r
library(microforge)
```

qPCR standard curve and convert Cq to copies:

```r
std <- calculate_std_curve(std_df)
samples_with_copies <- calculate_copies(sample_df, std)
```

Tidy plate data:

```r
p96 <- tidy_plate96(plate96_raw)
p384 <- tidy_plate384(plate384_raw)
```

Compute relative abundance (e.g., at family level):

```r
rel_ab <- calculate_relative_abundance(sqm_df, taxa_rank = "f")
```

Run DESeq2 pipeline:

```r
dds_res <- run_deseq2("dds.rds", "data/")
```

Note: refer to function help pages for argument details and examples (`?calculate_std_curve`, etc.).

## Testing & Vignettes

- Tests are in `tests/testthat/`; run them with:

```r
devtools::test()
```

- A vignette is available at `vignettes/intro.Rmd`. Build and view with:

```r
devtools::build_vignettes()
```

## Development & Contributing

- Follow existing style and add tests for new features (`testthat`).
- Update `NEWS.md` and the vignette when adding or changing behavior.
- Create small, focused PRs and include examples and unit tests.

## Documentation

Documentation is generated from roxygen2 comments and available in `man/`. Use `devtools::document()` to regenerate.

## License & Contact

- License: MIT (see `LICENSE`).
- Author: Jerry He — j.he@ubc.ca
- Package DESCRIPTION contains full metadata.

---

If you want, I can also add badges (CI, CRAN checks), a short usage vignette example, or more detailed function examples.