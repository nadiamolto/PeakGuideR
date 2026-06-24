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
