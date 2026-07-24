# Deduplicate a compound mass database by resolved compound identity

Collapses `compound_db` to one representative row per
`compound_identity_id` (see
[`resolve_compound_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_compound_identity.md)),
preferring the row with the most complete identifiers when several rows
share the same identity (for example the same compound reported under
`Source == "HMDB"` and `Source == "NORMAN"`).

## Usage

``` r
deduplicate_compound_identity(compound_db)
```

## Arguments

- compound_db:

  A compound mass database data.frame. If it does not already contain a
  `compound_identity_id` column,
  [`resolve_compound_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_compound_identity.md)
  is applied first.

## Value

A deduplicated data.frame with one row per `compound_identity_id`.
