# Resolve compound identity across a broad compound database and the standard-adduct library

Within each `neutral_mass_id`, links `compound_db`-derived candidates
(`compound_candidates`) with `standards_db`-derived candidates
(`standard_candidates`, typically from
[`match_standards_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_standards_by_mass.md))
using a hierarchy of identifiers, stopping at the first level that
produces a match: `InChIKey` -\> `InChI` -\> `SMILES` -\> shared
database ID (`HMDB` or `ChEBI`, only when the `compound_db` row's
`Source` matches that identifier type). Matching is 1:1: a given
candidate is linked to at most one counterpart, using the strongest
available identifier level.

A match on molecular formula alone (with or without an exact name match)
is **never** treated as identity. Such rows are kept as separate
`"broad_db_only"` / `"standards_only"` rows and cross-referenced through
`possible_link_formula_only` for manual review.

`identity_match_type` is bookkeeping only: it records which identifier
level produced a `"both"` fusion, but must not be used as a numeric
evidence score.

## Usage

``` r
resolve_cross_source_identity(compound_candidates, standard_candidates = NULL)
```

## Arguments

- compound_candidates:

  A data.frame of `compound_db`-derived candidates with at least
  `neutral_mass_id`, `candidate_source`, `candidate_db_id`,
  `candidate_name`, `candidate_formula`, `candidate_neutral_mass`,
  `candidate_ppm_error`, `candidate_inchi`, `candidate_inchikey`,
  `candidate_smiles`.

- standard_candidates:

  A data.frame of `standards_db`-derived candidates, typically the
  output of
  [`match_standards_by_mass()`](https://nadiamolto.github.io/PeakGuideR/reference/match_standards_by_mass.md),
  or `NULL`/empty if standard-adduct support is not being used.

## Value

A data.frame with one row per resolved candidate identity and columns:
`neutral_mass_id`, `source` (`"both"`, `"broad_db_only"` or
`"standards_only"`), `identity_match_type`,
`possible_link_formula_only`, `candidate_ppm_error`,
`candidate_neutral_mass`, `broad_db_name`, `broad_db_formula`,
`broad_db_inchikey`, `broad_db_inchi`, `broad_db_smiles`, `broad_db_id`,
`broad_db_source`, `standard_db_name`, `standard_db_formula`,
`standard_db_inchikey`, `standard_db_smiles`, `standard_db_compound_id`,
`standard_db_adducts`.

## Examples

``` r
# "Citric acid" and "Isocitric acid" share the formula C6H8O7, but only
# "Citric acid" shares a ChEBI id with the standard "CITRATE": formula
# alone must not fuse the two compounds.
compound_candidates <- data.frame(
  neutral_mass_id = c(1L, 1L),
  candidate_source = c("CHEBI", "CHEBI"),
  candidate_db_id = c("30769", "30887"),
  candidate_name = c("Citric acid", "Isocitric acid"),
  candidate_formula = c("C6H8O7", "C6H8O7"),
  candidate_neutral_mass = c(192.027, 192.027),
  candidate_ppm_error = c(0, 0),
  candidate_inchi = NA_character_,
  candidate_inchikey = NA_character_,
  candidate_smiles = NA_character_,
  stringsAsFactors = FALSE
)

standard_candidates <- data.frame(
  neutral_mass_id = 1L,
  candidate_source = "standards_db",
  candidate_name = "CITRATE",
  candidate_formula = "C6H8O7",
  candidate_neutral_mass = 192.027,
  candidate_ppm_error = 0,
  candidate_db_id = "P1A5",
  candidate_inchikey = "KRKNYBCHXYNGOX-UHFFFAOYSA-N",
  candidate_smiles = "OC(=O)CC(O)(CC(O)=O)C(O)=O",
  candidate_hmdb = "HMDB0000094",
  candidate_chebi = "30769",
  standard_db_adducts = "[M-H]-",
  stringsAsFactors = FALSE
)

resolve_cross_source_identity(compound_candidates, standard_candidates)
#>   neutral_mass_id        source identity_match_type possible_link_formula_only
#> 1               1          both           shared_id                       <NA>
#> 2               1 broad_db_only                <NA>                       <NA>
#>   candidate_ppm_error candidate_neutral_mass  broad_db_name broad_db_formula
#> 1                   0                192.027    Citric acid           C6H8O7
#> 2                   0                192.027 Isocitric acid           C6H8O7
#>   broad_db_inchikey broad_db_inchi broad_db_smiles broad_db_id broad_db_source
#> 1              <NA>           <NA>            <NA>       30769           CHEBI
#> 2              <NA>           <NA>            <NA>       30887           CHEBI
#>   standard_db_name standard_db_formula        standard_db_inchikey
#> 1          CITRATE              C6H8O7 KRKNYBCHXYNGOX-UHFFFAOYSA-N
#> 2             <NA>                <NA>                        <NA>
#>           standard_db_smiles standard_db_compound_id standard_db_adducts
#> 1 OC(=O)CC(O)(CC(O)=O)C(O)=O                    P1A5              [M-H]-
#> 2                       <NA>                    <NA>                <NA>
```
