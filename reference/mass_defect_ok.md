# Mass defect flag for neutral masses

Flags whether the mass defect of a neutral mass falls within the range
of mass defects observed among HMDB and ChEBI compounds of the same
integer mass. The flag is informative only: it does not change any
score.

The mass defect of each neutral mass is computed with
[`mass_defect_value()`](https://nadiamolto.github.io/PeakGuideR/reference/mass_defect_value.md)
and compared with the limits returned by
[`mass_defect_limits()`](https://nadiamolto.github.io/PeakGuideR/reference/mass_defect_limits.md).

## Usage

``` r
mass_defect_ok(neutral_mass)
```

## Arguments

- neutral_mass:

  Numeric vector of neutral masses.

## Value

A logical vector (same length as `neutral_mass`): `TRUE` when the mass
defect is within the limits, `FALSE` when it is outside, and `NA` when
the mass is missing, outside the range of the reference table (55-1300
Da) or belongs to an integer mass without reference compounds.
