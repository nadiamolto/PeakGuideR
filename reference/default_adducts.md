# Default adduct definitions

Returns the default adduct definitions used by PeakGuideR for the
selected ion mode.

## Usage

``` r
default_adducts(ion_mode = c("pos", "neg"))
```

## Arguments

- ion_mode:

  Ion mode. Either `"pos"` or `"neg"`. If `NULL`, adducts for both ion
  modes are returned.

## Value

A data.frame containing adduct definitions, including adduct names, ion
mode and mass shifts relative to the neutral molecule.

## Examples

``` r
default_adducts("pos")
#>         name mode       mass sign
#> 1     [M+H]+  pos   1.007276    1
#> 2    [M+Na]+  pos  22.989218    1
#> 3   [M+NH4]+  pos  18.033823    1
#> 4     [M+K]+  pos  38.963158    1
#> 5 [M+H-H2O]+  pos -17.003289    1
#> 6 [M+2Na-H]+  pos  44.971160    1
#> 7  [M+2K-H]+  pos  76.919040    1
default_adducts("neg")
#>          name mode       mass sign
#> 8      [M-H]-  neg  -1.007276    1
#> 9     [M+Cl]-  neg  34.968853    1
#> 10 [M-H-H2O]-  neg -19.017841    1
default_adducts(NULL)
#>          name mode       mass sign
#> 1      [M+H]+  pos   1.007276    1
#> 2     [M+Na]+  pos  22.989218    1
#> 3    [M+NH4]+  pos  18.033823    1
#> 4      [M+K]+  pos  38.963158    1
#> 5  [M+H-H2O]+  pos -17.003289    1
#> 6  [M+2Na-H]+  pos  44.971160    1
#> 7   [M+2K-H]+  pos  76.919040    1
#> 8      [M-H]-  neg  -1.007276    1
#> 9     [M+Cl]-  neg  34.968853    1
#> 10 [M-H-H2O]-  neg -19.017841    1
```
