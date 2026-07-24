# Standard-adduct recovery score for a single candidate row

Computes, for one standard-linked candidate compound, the fraction of
its `standards_db`-annotated adducts that were also detected
computationally for the same neutral mass: \$\$score = \|detected \cap
standard\| / \|standard\|\$\$

Returns `NA` (not `0`) when there is no standard-library compound within
tolerance for that neutral mass, since absence of an experimental
reference is not negative evidence.

When several standard-library compounds (isomers) fall within tolerance
of the same neutral mass, each is scored independently against the same
`detected_adducts`, which can help discriminate between isomers based on
which one's expected adducts were actually observed.

## Usage

``` r
compute_standard_adduct_recovery(standard_db_adducts, detected_adducts)
```

## Arguments

- standard_db_adducts:

  Character scalar: semicolon-separated adducts annotated in
  `standards_db` for this candidate compound, or `NA` if this candidate
  has no standards-library link.

- detected_adducts:

  Character vector of adducts detected computationally for the same
  `neutral_mass_id`.

## Value

A single numeric value in the range 0 to 1, or `NA_real_`.
