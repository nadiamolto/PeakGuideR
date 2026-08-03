# Resolve a per-row chemical identity key for isomer-ambiguity detection

Collapses candidate rows that refer to the same chemical identity before
`ambiguous_isomeric` is derived from their `priority_score`s, so that
two rows describing the *same* compound never look like two competing
candidates. This matters most for
`hypothesis_origin == "single_adduct_isolated"` rows: two adducts of one
compound that were each inferred independently (because they never
merged into a shared adduct family) otherwise produce two rows with an
identical `priority_score`, which a naive top-row-vs-second-row
comparison flags as ambiguous even though there is no second identity at
all.

Uses the same identifier hierarchy as
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md):
`InChIKey` -\> `InChI` -\> `SMILES` -\> database id (`broad_db_id` /
`standard_db_compound_id`) -\> name + formula, stopping at the first
non-missing identifier (falling back from the `broad_db_*` to the
`standard_db_*` columns at each level, since a row can carry either or
both depending on `source`). Rows with none of these populated (no
resolved database candidate at all) each receive their own unique key,
so they are never collapsed with one another.

## Usage

``` r
resolve_identity_key(combined)
```

## Arguments

- combined:

  A data.frame with `broad_db_inchikey`, `broad_db_inchi`,
  `broad_db_smiles`, `broad_db_id`, `broad_db_name`, `broad_db_formula`,
  `standard_db_inchikey`, `standard_db_smiles`,
  `standard_db_compound_id`, `standard_db_name`, `standard_db_formula`
  columns.

## Value

A character vector, the same length as `nrow(combined)`, with one
identity key per row.
