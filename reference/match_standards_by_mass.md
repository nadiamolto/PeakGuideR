# Match inferred neutral masses against the standard-adduct library by mass

Retrieves standard-adduct library compounds within `ppm_tol` of each
inferred neutral mass, symmetrically to how `compound_db` is matched by
mass elsewhere in PeakGuideR.

Unlike the previous implementation, this function never joins on
identifiers hanging off a `compound_db` candidate row: `standards_db` is
treated as an independent evidence source that is retrieved purely by
neutral mass, ion mode and (optionally) matrix. Cross-source identity is
resolved afterwards by
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md).

Because `standards_db` typically has several adduct rows per compound
(all sharing the same `NEUTRAL_MONOISOTOPIC_MASS`), matches are
collapsed to one row per standard compound and per `neutral_mass_id`,
with `standard_db_adducts` listing every annotated adduct for that
compound (used downstream to compute `standard_adduct_recovery_score`).
When two or more standard-library compounds (isomers) fall within
tolerance of the same neutral mass, they are kept as separate rows
rather than merged.

## Usage

``` r
match_standards_by_mass(
  neutral_masses,
  standards_db,
  ion_mode = c("pos", "neg"),
  matrix = NULL,
  ppm_tol = 5
)
```

## Arguments

- neutral_masses:

  A data.frame with at least `neutral_mass_id` and
  `neutral_mass_consensus`.

- standards_db:

  Standard-adduct library data.frame. Must contain
  `NEUTRAL_MONOISOTOPIC_MASS` and `POLARITY`.

- ion_mode:

  `"pos"` or `"neg"`.

- matrix:

  Optional matrix name. When supplied and `standards_db` contains a
  `matrix` column, only rows whose `matrix` value contains `matrix`
  (case-insensitive) are kept. Use `NULL` to skip matrix filtering.

- ppm_tol:

  PPM tolerance for neutral-mass matching.

## Value

A data.frame with one row per `neutral_mass_id` and standard compound
candidate, with columns `neutral_mass_id`, `candidate_source`
(`"standards_db"`), `candidate_name`, `candidate_formula`,
`candidate_neutral_mass`, `candidate_ppm_error`, `candidate_db_id`,
`candidate_inchikey`, `candidate_smiles`, `candidate_hmdb`,
`candidate_chebi` and `standard_db_adducts` (semicolon-separated adducts
annotated for that compound in `standards_db`, after ion-mode/matrix
filtering).

## Examples

``` r
standards_db <- load_standards_adduct_library(quiet = TRUE)

neutral_masses <- data.frame(
  neutral_mass_id = 1L,
  neutral_mass_consensus = standards_db$NEUTRAL_MONOISOTOPIC_MASS[1]
)

match_standards_by_mass(
  neutral_masses,
  standards_db,
  ion_mode = "neg",
  matrix = "HCCA"
)
#> # A tibble: 1 × 12
#>   neutral_mass_id candidate_source candidate_name candidate_formula
#>             <int> <chr>            <chr>          <chr>            
#> 1               1 standards_db     CITRATE        C6H8O7           
#> # ℹ 8 more variables: candidate_neutral_mass <dbl>, candidate_ppm_error <dbl>,
#> #   candidate_db_id <chr>, candidate_inchikey <chr>, candidate_smiles <chr>,
#> #   candidate_hmdb <chr>, candidate_chebi <chr>, standard_db_adducts <chr>
```
