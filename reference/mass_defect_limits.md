# Reference limits of the mass defect for a neutral mass

Looks up, in the internal mass defect reference table, the limits of the
mass defect observed among HMDB and ChEBI compounds of the same integer
mass. The table has one row per integer mass; each neutral mass uses the
row of `floor(mass)`.

Masses outside the range of the table (55-1300 Da), integer masses
without reference compounds and missing masses return `NA` limits.

## Usage

``` r
mass_defect_limits(mass)
```

## Arguments

- mass:

  Numeric vector of neutral masses.

## Value

A data.frame with one row per element of `mass` and columns `lower` and
`upper`.
