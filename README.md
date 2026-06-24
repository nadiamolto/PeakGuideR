# PeakGuideR

PeakGuideR is an R package for evidence-based annotation support in spatial 
metabolomics peak data.

The package combines complementary evidence layers, including isotope morphology,
carbon isotope-ratio support, elemental isotope-pattern support, adduct-family
grouping, neutral-mass inference and database matching.

PeakGuideR does **not** provide definitive compound identification. Its outputs
should be interpreted as putative annotation evidence that can guide downstream
validation.

## Installation

PeakGuideR can be installed from GitHub (keep "vignettes=TRUE" to install it 
with vignettes (recomended):

```r
remotes::install_github("nadiamolto/PeakGuideR", build_vignettes=TRUE, dependencies=TRUE)
```

The package vignettes can then opened with:

```r 
browseVignettes("PeakGuideR")
```

## Basic workflow

PeakGuideR accepts either a rMSI2 peak matrix object or a Cardinal MSI
object.

### From a PeakGuideR peak matrix

```r
res <- run_peakguider_workflow(
  pkm = pkm,
  ion_mode = "pos",
  matrix = "HCCA")
```

### From a Cardinal MSI object

```r
res <- run_peakguider_workflow(
  pkm = MSImagingExperimentObject,
  ion_mode = "pos",
  matrix = "HCCA")
```

Internally, Cardinal objects are converted to peak-matrix format
using `cardinal_to_peakmatrix()`.

The workflow returns a `peakguider_workflow` object:

```r
class(pg_res)
names(pg_res)
```

Main outputs include:

- `pg_res$feature_summary`: one row per m/z feature, summarizing isotope, EIPS and adduct-family evidence.
- `pg_res$neutral_mass_candidates`: one row per inferred neutral mass and compound candidate.
- `pg_res$relation_table`: pairwise relationships between features.
- `pg_res$adduct_families`: inferred adduct families and neutral masses.
- `pg_res$morph_results`, `pg_res$cir_results`, `pg_res$eips_results`, `pg_res$adduct_edges`: intermediate evidence tables regarding isotope and adduct module.

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

- `idx`: feature index.
- `mz`: feature m/z value.
- `is_c13_m0`: whether the feature is supported as a monoisotopic peak.
- `c13_m1_idx`: feature index of the putative C13 M+1 partner.
- `c13_score`: carbon isotope-ratio support score.
- `has_eips`: whether elemental isotope-pattern support was detected.
- `eips_elements`: supported isotope-pattern elements.
- `has_adduct_family`: whether the feature belongs to an inferred adduct family.
- `adduct_roles`: inferred adduct role or roles for the feature.
- `neutral_mass_consensus`: inferred neutral mass associated with the main adduct family.

### Neutral-mass candidates

`neutral_mass_candidates` summarizes inferred neutral masses and database
candidate matches.

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

- `neutral_mass_id`: internal neutral-mass group identifier.
- `neutral_mass_consensus`: inferred neutral mass.
- `inferred_adducts`: adducts supporting the inferred neutral mass.
- `candidate_source`: source database of the candidate match.
- `candidate_db_id`: database identifier.
- `candidate_name`: candidate compound name.
- `candidate_formula`: candidate molecular formula.
- `candidate_neutral_mass`: neutral monoisotopic mass from the database.
- `candidate_ppm_error`: mass error between the inferred neutral mass and the candidate.
- `has_standard_compound_match`: whether the candidate compound is present in the standard-adduct library.
- `has_standard_adduct_match`: whether the inferred adduct is supported for that compound in the standard-adduct library.

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

## Databases

PeakGuideR includes small example databases for testing and demonstration:

- `compound_mass_database_example.rds`
- `standards_adduct_library_example.rds`

These example databases are automatically used when no external databases are
provided.

For real analyses, users should download the full non-commercial annotation
databases from Zenodo:

- `compound_mass_database_noncommercial.rds`
- `standards_adduct_library_noncommercial.rds`

Zenodo record: `10.5281/zenodo.20705395`

These full databases are distributed separately because they include third-party
compound metadata with non-commercial licensing restrictions.

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
  eips_n_table = eips_n_table,
  eips_table = eips_table,
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
  eips_n_table = eips_n_table,
  eips_table = eips_table,
  compound_db = dbs$compound_db,
  standards_db = dbs$standards_db
)
```

### Standard-adduct support

Standard-adduct support is currently matrix-specific. For `matrix = "HCCA"`,
PeakGuideR can check whether a candidate compound and inferred adduct are present
in the HCCA standard-adduct library.

In this context, HCCA refers to alpha-cyano-4-hydroxycinnamic acid. The
standard-adduct library was generated using an HCCA-based solid ionic matrix
combined with N,N-diethylaniline (DEA) deposited by low-temperature thermal evaporation.

Standard-adduct support should therefore be interpreted as matrix- and
preparation-specific evidence, not as universal adduct evidence for all MALDI
matrices or deposition protocols.

## Important interpretation note

PeakGuideR outputs are evidence layers for annotation support.

A database match in `neutral_mass_candidates` means that the inferred neutral
mass is compatible with a candidate compound within the selected ppm tolerance.
It does not constitute definitive compound identification.

Experimental validation, MS/MS, standards or orthogonal evidence may be required
for confident identification.

## Data licensing notice

PeakGuideR includes small example compound and standard-adduct databases for
testing and demonstration. Full annotation databases are distributed separately
for non-commercial research use.

These datasets include records derived from ChEBI, NORMAN and HMDB.
ChEBI and NORMAN/SusDat records are distributed under CC BY 4.0 terms.
HMDB-derived records are subject to CC BY-NC 4.0 non-commercial use restrictions
according to the original HMDB licensing terms.

The databases are used only to retrieve putative mass-matched candidates and do
not constitute compound identification. Users are responsible for ensuring that
their use complies with the original providers' licenses.

## Citation

If you use PeakGuideR, please cite the package and the corresponding Zenodo
database record when using the full annotation databases.

Zenodo database record: `https://doi.org/10.5281/zenodo.20705395`
