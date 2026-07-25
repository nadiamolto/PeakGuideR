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

`standards_db` can annotate adduct forms - most commonly dimers
(`[2M+H]+`, `[2M+Na]+`, ...) and zero-shift `[M]+` - that fall outside
the single-feature, monomer-vs-neutral-mass paradigm PeakGuideR's own
adduct detection uses
([`adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/adduct_candidates.md)/[`build_single_adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_single_adduct_candidates.md)
always solve `neutral_mass = mz - adduct_mass` for one feature; they
never pair two masses to reconstruct a dimer). Such forms can never
appear in `detected_adducts`, regardless of how complete the actual
recovery was, so counting them in the denominator would deflate the
score for reasons unrelated to detection quality. `testable_adducts`
restricts the comparison to adduct forms the detection pipeline could in
principle have produced; entries of `standard_db_adducts` outside that
set are dropped before scoring, not treated as missed.

## Usage

``` r
compute_standard_adduct_recovery(
  standard_db_adducts,
  detected_adducts,
  testable_adducts
)
```

## Arguments

- standard_db_adducts:

  Character scalar: semicolon-separated adducts annotated in
  `standards_db` for this candidate compound, or `NA` if this candidate
  has no standards-library link.

- detected_adducts:

  Character vector of adducts detected computationally for the same
  `neutral_mass_id`.

- testable_adducts:

  Character vector of adduct names the detection pipeline is
  structurally capable of producing for this `ion_mode` (typically
  `default_adducts(ion_mode)$name`). `standard_db_adducts` entries
  outside this set (e.g. dimers, `[M]+`) are excluded from both the
  numerator and the denominator before scoring. If the filtered set is
  empty, returns `NA_real_` (an impossible comparison, not zero
  recovery).

## Value

A single numeric value in the range 0 to 1, or `NA_real_`.
