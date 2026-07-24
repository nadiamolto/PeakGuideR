# Load the PeakGuideR compound mass database

Loads the example compound mass database included with PeakGuideR. The
full non-commercial database is distributed separately and can be
supplied manually to the workflow through the `compound_db` argument.

## Usage

``` r
load_compound_mass_database(quiet = FALSE)
```

## Arguments

- quiet:

  Logical. If `FALSE`, prints a message.

## Value

A data.frame.

## Examples

``` r
compound_db <- load_compound_mass_database(quiet = TRUE)
head(compound_db)
#> # A tibble: 6 × 11
#>    ...1 Source DB_ID Name             MolecularFormula MonoisotopicMass StdInChI
#>   <dbl> <chr>  <chr> <chr>            <chr>                       <dbl> <chr>   
#> 1    84 CHEBI  160   (1S,2S)-1-hydro… C6H8O7                       192. InChI=1…
#> 2   265 CHEBI  18350 1(2H)-Isoquinol… C9H7NO                       145. InChI=1…
#> 3   341 CHEBI  17472 1-Methyl-4-phen… C12H15NO                     189. InChI=1…
#> 4   362 CHEBI  16060 10-Deoxysarpagi… C19H22N2O                    294. InChI=1…
#> 5   486 CHEBI  915   2,4-Dihydroxyhe… C7H10O6                      190. InChI=1…
#> 6   571 CHEBI  17305 2-Dehydro-3-deo… C6H8O7                       192. InChI=1…
#> # ℹ 4 more variables: StdInChIKey <chr>, SMILES <chr>, Kegg <chr>, dup <lgl>
```
