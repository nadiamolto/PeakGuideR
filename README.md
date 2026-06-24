<img src="man/figures/PeakGuideR_logo.png" alt="PeakGuideR logo showing stylized MSI peaks and spatial metabolomics annotation" align="right" width="280" style="margin-top: 25px;"/>
# PeakGuideR

<!-- badges: start -->
<!-- badges: end -->

**PeakGuideR** is an R package for evidence-based annotation support in spatial metabolomics peak data.

The package combines complementary evidence layers, including isotope morphology, carbon isotope-ratio support, elemental isotope-pattern support, adduct-family grouping, neutral-mass inference and database matching.

> **Important**  
> PeakGuideR does **not** provide definitive compound identification. Its outputs should be interpreted as putative annotation evidence that can guide downstream validation.



## Overview

PeakGuideR integrates several evidence layers for spatial metabolomics peak annotation:

| Evidence layer | Purpose |
|---|---|
| Isotope morphology | Detects isotope-related peaks using mass spacing and spatial similarity |
| Carbon isotope-ratio support | Evaluates whether M+1/M0 ratios are compatible with carbon isotope behaviour |
| Elemental isotope-pattern support | Evaluates elemental isotope-pattern evidence such as N, O, S, Cl and Br |
| Adduct-family grouping | Groups related peaks into putative adduct families |
| Neutral-mass inference | Infers consensus neutral masses from compatible adduct assignments |
| Database matching | Reports putative mass-matched candidate compounds |

---

## Installation

PeakGuideR can be installed from GitHub. Keep `build_vignettes = TRUE` to install the package with vignettes:

```r
remotes::install_github(
  "nadiamolto/PeakGuideR",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

The package vignettes can then be opened with:

```r
browseVignettes("PeakGuideR")
```

---

## Basic workflow

PeakGuideR accepts either a rMSI2 peak matrix object or a Cardinal MSI object.

### From a PeakGuideR peak matrix

```r
res <- run_peakguider_workflow(
  pkm = pkm,
  ion_mode = "pos",
  matrix = "HCCA"
)
```

### From a Cardinal MSI object

```r
res <- run_peakguider_workflow(
  pkm = MSImagingExperimentObject,
  ion_mode = "pos",
  matrix = "HCCA"
)
```

Internally, Cardinal objects are converted to peak-matrix format using `cardinal_to_peakmatrix()`.

The workflow returns a `peakguider_workflow` object:

```r
class(res)
names(res)
```

---

## Main outputs

| Output | Description |
|---|---|
| `res$feature_summary` | One row per m/z feature, summarizing isotope, EIPS and adduct-family evidence |
| `res$neutral_mass_candidates` | One row per inferred neutral mass and compound candidate |
| `res$relation_table` | Pairwise relationships between features |
| `res$adduct_families` | Inferred adduct families and neutral masses |
| `res$morph_results` | Isotope morphology results |
| `res$cir_results` | Carbon isotope-ratio validation results |
| `res$eips_results` | Elemental isotope-pattern support results |
| `res$adduct_edges` | Pairwise adduct-compatible relationships |

---

## Exploring the output

### Feature-level evidence

`feature_summary` summarizes evidence at the detected feature level.

```r
res$feature_summary |>
  dplyr::select(
    idx,
    mz,
    is_c13_m0,
    c13_m1_idx,
    c13_score,
    has_eips,
    eips_elements,
    has_adduct_family,
    adduct_roles,
    main_adduct_role,
    neutral_mass_consensus
  ) |>
  head()
```

Useful columns include:

| Column | Description |
|---|---|
| `idx` | Feature index |
| `mz` | Feature m/z value |
| `is_c13_m0` | Whether the feature is supported as a monoisotopic peak |
| `c13_m1_idx` | Feature index of the putative C13 M+1 partner |
| `c13_score` | Carbon isotope-ratio support score |
| `has_eips` | Whether elemental isotope-pattern support was detected |
| `eips_elements` | Supported isotope-pattern elements |
| `has_adduct_family` | Whether the feature belongs to an inferred adduct family |
| `adduct_roles` | Inferred adduct role or roles for the feature |
| `neutral_mass_consensus` | Inferred neutral mass associated with the main adduct family |

---

### Neutral-mass candidates

`neutral_mass_candidates` summarizes inferred neutral masses and database candidate matches.

```r
res$neutral_mass_candidates |>
  dplyr::select(
    neutral_mass_id,
    neutral_mass_consensus,
    inferred_adducts,
    candidate_source,
    candidate_db_id,
    candidate_name,
    candidate_formula,
    candidate_neutral_mass,
    candidate_ppm_error,
    has_standard_compound_match,
    has_standard_adduct_match
  ) |>
  dplyr::filter(!is.na(candidate_name))
```

Useful columns include:

| Column | Description |
|---|---|
| `neutral_mass_id` | Internal neutral-mass group identifier |
| `neutral_mass_consensus` | Inferred neutral mass |
| `inferred_adducts` | Adducts supporting the inferred neutral mass |
| `candidate_source` | Source database of the candidate match |
| `candidate_db_id` | Database identifier |
| `candidate_name` | Candidate compound name |
| `candidate_formula` | Candidate molecular formula |
| `candidate_neutral_mass` | Neutral monoisotopic mass from the database |
| `candidate_ppm_error` | Mass error between the inferred neutral mass and the candidate |
| `has_standard_compound_match` | Whether the candidate compound is present in the standard-adduct library |
| `has_standard_adduct_match` | Whether the inferred adduct is supported for that compound in the standard-adduct library |

---

### Summary of candidate matches

```r
res$neutral_mass_candidates |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_neutral_masses = dplyr::n_distinct(neutral_mass_id),
    n_rows_with_candidate = sum(!is.na(candidate_name)),
    n_neutral_masses_with_candidate = dplyr::n_distinct(
      neutral_mass_id[!is.na(candidate_name)]
    ),
    n_with_standard_compound_match = sum(has_standard_compound_match %in% TRUE),
    n_with_standard_adduct_match = sum(has_standard_adduct_match %in% TRUE)
  )
```

---

## Databases

PeakGuideR includes small example databases for testing and demonstration:

| Example database | Description |
|---|---|
| `compound_mass_database_example.rds` | Small compound mass database for testing and demonstration |
| `standards_adduct_library_example.rds` | Small standard-adduct library for testing and demonstration |

These example databases are automatically used when no external databases are provided.

For real analyses, users should download the full non-commercial annotation databases from Zenodo:

| Full database | Description |
|---|---|
| `compound_mass_database_noncommercial.rds` | Full non-commercial compound mass database |
| `standards_adduct_library_noncommercial.rds` | Full non-commercial standard-adduct library |

Zenodo record: `10.5281/zenodo.20705395`

These full databases are distributed separately because they include third-party compound metadata with non-commercial licensing restrictions.

---

### Loading full databases manually

After downloading the files from Zenodo, load them in R:

```r
compound_db <- readRDS("path/to/compound_mass_database_noncommercial.rds")
standards_db <- readRDS("path/to/standards_adduct_library_noncommercial.rds")
```

Then pass them to the workflow:

```r
res <- run_peakguider_workflow(
  pkm = pkm,
  ion_mode = "pos",
  matrix = "HCCA",
  compound_db = compound_db,
  standards_db = standards_db
)
```

Alternatively, if available in your installed version, use:

```r
dbs <- load_peakguider_databases(
  compound_db_path = "path/to/compound_mass_database_noncommercial.rds",
  standards_db_path = "path/to/standards_adduct_library_noncommercial.rds"
)

res <- run_peakguider_workflow(
  pkm = pkm,
  ion_mode = "pos",
  matrix = "HCCA",
  compound_db = dbs$compound_db,
  standards_db = dbs$standards_db
)
```

---

## Standard-adduct support

Standard-adduct support is currently matrix-specific. For `matrix = "HCCA"`, PeakGuideR can check whether a candidate compound and inferred adduct are present in the HCCA standard-adduct library.

In this context, HCCA refers to alpha-cyano-4-hydroxycinnamic acid. The standard-adduct library was generated using an HCCA-based solid ionic matrix combined with N,N-diethylaniline (DEA) deposited by low-temperature thermal evaporation.

> **Interpretation note**  
> Standard-adduct support should be interpreted as matrix- and preparation-specific evidence, not as universal adduct evidence for all MALDI matrices or deposition protocols.

---

## Important interpretation note

PeakGuideR outputs are evidence layers for annotation support.

A database match in `neutral_mass_candidates` means that the inferred neutral mass is compatible with a candidate compound within the selected ppm tolerance. It does not constitute definitive compound identification.

Experimental validation, MS/MS, standards or orthogonal evidence may be required for confident identification.

---

## Data licensing notice

PeakGuideR includes small example compound and standard-adduct databases for testing and demonstration. Full annotation databases are distributed separately for non-commercial research use.

These datasets include records derived from ChEBI, NORMAN and HMDB. ChEBI and NORMAN/SusDat records are distributed under CC BY 4.0 terms. HMDB-derived records are subject to CC BY-NC 4.0 non-commercial use restrictions according to the original HMDB licensing terms.

The databases are used only to retrieve putative mass-matched candidates and do not constitute compound identification. Users are responsible for ensuring that their use complies with the original providers' licenses.

---

## Citation

If you use PeakGuideR, please cite the package and the corresponding Zenodo database record when using the full annotation databases.

Zenodo database record: `https://doi.org/10.5281/zenodo.20705395`
