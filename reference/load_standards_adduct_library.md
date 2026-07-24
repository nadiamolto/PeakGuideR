# Load the PeakGuideR standard adduct library

Loads the example standard adduct library included with PeakGuideR. The
full non-commercial standard adduct library is distributed separately
and can be supplied manually to the workflow through the `standards_db`
argument.

## Usage

``` r
load_standards_adduct_library(quiet = FALSE)
```

## Arguments

- quiet:

  Logical. If `FALSE`, prints a message.

## Value

A data.frame.

## Examples

``` r
standards_db <- load_standards_adduct_library(quiet = TRUE)
head(standards_db)
#> # A tibble: 6 × 22
#>   COMPOUND_ID COMPOUND Master_List_NAME POLARITY matrix          matrix_compound
#>   <chr>       <chr>    <chr>            <chr>    <chr>           <chr>          
#> 1 P1A5        CITRATE  CITRATE          neg      HCCA-DEA solid… alpha-cyano-4-…
#> 2 P1A5        CITRATE  CITRATE          neg      HCCA-DEA solid… alpha-cyano-4-…
#> 3 P1A5        CITRATE  CITRATE          neg      HCCA-DEA solid… alpha-cyano-4-…
#> 4 P1A5        CITRATE  CITRATE          neg      HCCA-DEA solid… alpha-cyano-4-…
#> 5 P1A5        CITRATE  CITRATE          neg      HCCA-DEA solid… alpha-cyano-4-…
#> 6 P1A5        CITRATE  CITRATE          pos      HCCA-DEA solid… alpha-cyano-4-…
#> # ℹ 16 more variables: matrix_additive <chr>, deposition_method <chr>,
#> #   adduct <chr>, mz_theo <dbl>, MOLECULAR_FORMULA <chr>,
#> #   NEUTRAL_MONOISOTOPIC_MASS <dbl>, SOLVENT <chr>, HMDB_clean <chr>,
#> #   CAS_ID <chr>, SMILES <chr>, KEGG_ID_CSID <chr>, METLIN_ID <chr>,
#> #   PC_CID <chr>, PC_SID <chr>, ChEBI <chr>, InCHIKey <chr>
```
