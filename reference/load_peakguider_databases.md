# Load external PeakGuideR databases

Loads full PeakGuideR annotation databases downloaded separately, for
example from Zenodo.

## Usage

``` r
load_peakguider_databases(compound_db_path, standards_db_path = NULL)
```

## Arguments

- compound_db_path:

  Path to `compound_mass_database_noncommercial.rds`.

- standards_db_path:

  Optional path to `standards_adduct_library_noncommercial.rds`.

## Value

A list with `compound_db` and `standards_db`.

## Examples

``` r
if (FALSE) { # \dontrun{
dbs <- load_peakguider_databases(
  compound_db_path = "path/to/compound_mass_database_noncommercial.rds",
  standards_db_path = "path/to/standards_adduct_library_noncommercial.rds"
)

data(example_pkm, package = "PeakGuideR")

res <- run_peakguider_workflow(
  pkm = example_pkm,
  ion_mode = "pos",
  matrix = "HCCA",
  compound_db = dbs$compound_db,
  standards_db = dbs$standards_db
)
} # }
```
