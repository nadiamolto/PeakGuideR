# Mass defect of a neutral mass

Computes the mass defect of neutral masses as the decimal part of the
mass, `mass - floor(mass)`, a value in `[0, 1)`. The integer part,
`floor(mass)`, is the integer mass the neutral mass belongs to (e.g.
130.0 and 130.9 both belong to 130).

## Usage

``` r
mass_defect_value(mass)
```

## Arguments

- mass:

  Numeric vector of neutral masses.

## Value

Numeric vector of mass defects (same length as `mass`).
