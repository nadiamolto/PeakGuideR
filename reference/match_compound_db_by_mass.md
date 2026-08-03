# Match inferred neutral masses against a compound mass database

Retrieves `compound_db` candidates within `ppm_tol` of each inferred
neutral mass. `compound_db` is deduplicated by resolved compound
identity (see
[`resolve_compound_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_compound_identity.md)
and
[`deduplicate_compound_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/deduplicate_compound_identity.md))
before matching, so that compounds reported under more than one `Source`
(for example the same compound in both `HMDB` and `NORMAN`) are matched
only once per neutral mass.

## Usage

``` r
match_compound_db_by_mass(
  neutral_masses,
  compound_db,
  ppm_tol = 5,
  top_n = 10L
)
```

## Arguments

- neutral_masses:

  A data.frame with at least `neutral_mass_id` and
  `neutral_mass_consensus`.

- compound_db:

  Compound mass database. Must contain a numeric `MonoisotopicMass`
  column.

- ppm_tol:

  PPM tolerance for neutral-mass matching.

- top_n:

  Maximum number of candidates kept per `neutral_mass_id`, ranked by
  ascending ppm error. Use `NULL` to retain all candidates within
  `ppm_tol`. This cut is blind to `standards_db`: in a mass window
  crowded with more than `top_n` isomers/isobars, the one candidate that
  would have linked to a standards-library compound can be truncated
  away before
  [`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md)
  ever sees it.
  [`build_candidate_annotations()`](https://nadiamolto.github.io/PeakGuideR/reference/build_candidate_annotations.md)
  and
  [`build_single_adduct_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/build_single_adduct_candidates.md)
  compensate for this explicitly, via
  [`expand_compound_candidates_for_standards()`](https://nadiamolto.github.io/PeakGuideR/reference/expand_compound_candidates_for_standards.md),
  for the `neutral_mass_id` values where a standard was actually found -
  this function itself applies `top_n` unconditionally.

## Value

A data.frame with one row per `neutral_mass_id` and compound candidate:
`neutral_mass_id`, `candidate_source`, `candidate_db_id`,
`candidate_name`, `candidate_formula`, `candidate_neutral_mass`,
`candidate_inchi`, `candidate_inchikey`, `candidate_smiles`,
`candidate_kegg`, `candidate_mass`, `candidate_ppm_error`,
`compound_identity_id`.
