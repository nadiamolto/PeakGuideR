# Resolve compound identity within a compound mass database

Assigns a deterministic `compound_identity_id` to each row of
`compound_db` by grouping rows that represent the same underlying
compound, so that downstream matching (for example against
`standards_db`) is never performed against duplicate rows for the same
compound.

Rows are grouped by exact `StdInChIKey` match first. Rows without a
usable `StdInChIKey` are then grouped by exact `StdInChI` match. Rows
lacking both identifiers each receive their own identity, since there is
no reliable identifier to group them on (molecular formula alone is not
treated as proof of identity; see
[`resolve_cross_source_identity()`](https://nadiamolto.github.io/PeakGuideR/reference/resolve_cross_source_identity.md)).

This does **not** use the existing `dup` column in `compound_db`, whose
derivation criteria are unknown and unverifiable; deduplication is
always recomputed deterministically from `StdInChIKey`/`StdInChI`.

## Usage

``` r
resolve_compound_identity(compound_db)
```

## Arguments

- compound_db:

  A compound mass database data.frame. If it contains `StdInChIKey`
  and/or `StdInChI` columns, these are used for grouping.

## Value

`compound_db` with an added integer column `compound_identity_id`.

## Examples

``` r
compound_db <- load_compound_mass_database(quiet = TRUE)
compound_db <- resolve_compound_identity(compound_db)
length(unique(compound_db$compound_identity_id))
#> [1] 1458
```
