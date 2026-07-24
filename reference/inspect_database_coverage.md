# Inspect identifier coverage of a compound mass database

Reports, per `Source`, the percentage of rows with a non-empty
`StdInChIKey`, `StdInChI` and `SMILES`. This is useful before relying on
identifier-coverage assumptions (for example that a given `Source` never
has `StdInChIKey`), since coverage may differ between the bundled
example database and the full non-commercial database.

## Usage

``` r
inspect_database_coverage(compound_db)
```

## Arguments

- compound_db:

  A compound mass database data.frame.

## Value

A data.frame with one row per `Source`, columns `n` (row count) and
`pct_StdInChIKey`, `pct_StdInChI`, `pct_SMILES` (percentage of non-empty
values, rounded to 1 decimal).

## Examples

``` r
compound_db <- load_compound_mass_database(quiet = TRUE)
inspect_database_coverage(compound_db)
#> # A tibble: 3 × 5
#>   Source     n pct_StdInChIKey pct_StdInChI pct_SMILES
#>   <chr>  <int>           <dbl>        <dbl>      <dbl>
#> 1 CHEBI    403               0         89.1          0
#> 2 HMDB     247             100        100          100
#> 3 NORMAN   812             100        100          100
```
