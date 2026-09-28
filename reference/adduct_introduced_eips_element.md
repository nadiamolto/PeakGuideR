# EIPS-tracked element introduced by an adduct

Maps an adduct name (as in `default_adducts()$name`) to the elemental
isotope-pattern (EIPS) element it adds to the ion, if any. EIPS tracks
N, O, S, Cl and Br; only adducts that carry one of these into the ion
have an entry.

In the current
[`default_adducts()`](https://nadiamolto.github.io/PeakGuideR/reference/default_adducts.md)
table `[M+Cl]-` is the only such adduct (H, Na, K and NH4 add none of
the tracked elements), so a `35Cl`/`37Cl` isotope pattern on a feature
annotated as `[M+Cl]-` is explained by the adduct itself and is not
evidence that the neutral molecule contains chlorine. This is a fixed
internal lookup for the default adduct table, not a formula parser.

## Usage

``` r
adduct_introduced_eips_element(adduct)
```

## Arguments

- adduct:

  Character vector of adduct names.

## Value

Character vector the same length as `adduct`: the introduced element, or
`NA_character_` where the adduct introduces none.
