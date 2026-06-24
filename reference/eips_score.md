# Score elemental isotope-pattern candidates

Evaluates non-carbon isotope-pattern candidates by comparing the
observed isotope ratio against theoretical element-count ratios for N,
O, S, Cl and Br.

The function can optionally add formula support by comparing inferred
neutral masses with a database-derived EIPS lookup table.

Two optional reliability gates can be applied:

- For N/O/S, EIPS is only evaluated if the same `idx_M0` has a
  CIR-validated C13_M1 isotope.

- EIPS is not evaluated for pairs already detected as C13_M2 for an
  `idx_M0` with CIR-validated C13_M1.

For low-abundance N/O/S isotope channels, the function can require
independent C13 support for the same monoisotopic feature. This reduces
spurious elemental isotope assignments in noisy MSI peak tables.

In addition, inferred atom counts (`n_hat`) are bounded by the maximum
number of atoms observed in the database for that element. Candidates
exceeding that bound are discarded before final selection.

Expected isotope peak order for the elemental channel evaluated here. N
is evaluated at M+1; O, S, Cl and Br are evaluated at M+2.

EIPS does not perform the initial mass-difference search. It evaluates
elemental isotope candidates already detected by
[`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md).

## Usage

``` r
eips_score(
  result,
  pkm,
  eips_n_table = NULL,
  eips_table = NULL,
  ion_mode = c("pos", "neg"),
  adducts = NULL,
  min_score_final = 0.6,
  eips_rel_tol = 0.3,
  ratio_method = c("sum", "median", "mean"),
  min_quantile = 0.01,
  ppm_neutral = 5,
  n_window = 1L,
  top_n_hat = 5L,
  cir_df = NULL,
  morph_df = NULL,
  require_c13_for = c("N", "O", "S"),
  exclude_c13_m2 = TRUE,
  return_debug = FALSE
)
```

## Arguments

- result:

  Output of
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md)
  (must contain at least `idx_M0`, `idx_cand`, `score_final`, `element`,
  `iso_type`).

- pkm:

  Peak matrix list with `mass` and `intensity`.

- eips_n_table:

  Optional precomputed theoretical lookup table with columns `element`,
  `k`, `delta`, `n` and `R_theo`. If `NULL`, the internal PeakGuideR
  `eips_n_table` object is used.

- eips_table:

  Optional database-derived lookup table with columns `mz_mono`,
  `element` and `n_el`. If `NULL`, the internal PeakGuideR `eips_table`
  object is used.

- ion_mode:

  `"pos"` or `"neg"` (used only when formula support is enabled).

- adducts:

  Optional data.frame with adduct definitions. If `NULL`,
  [`default_adducts()`](https://nadiamolto.github.io/PeakGuideR/reference/default_adducts.md)
  is used. The table must contain `name`, `mode` and `mass`, where
  `mass` is the net m/z shift relative to the neutral mass.

- min_score_final:

  Minimum morphology score required to consider a pair.

- eips_rel_tol:

  Relative tolerance used as validity threshold and for ratio score
  scaling.

- ratio_method:

  `"sum"`, `"mean"`, or `"median"` aggregation of isotope ratios.
  `"sum"` computes sum(I_iso) / sum(I_M0) over selected pixels and is
  recommended for MSI isotope-ratio estimation.

- min_quantile:

  Quantile used for intensity masking. Default= 0.01.

- ppm_neutral:

  Ppm tolerance used when evaluating neutral mass consistency for
  elemental isotope-pattern support.

- n_window:

  Allowed deviation between inferred `n_hat` and formula `n_el` when
  formula support is evaluated.

- top_n_hat:

  Integer. Keep top-k best `n_hat` candidates by relative error.

- cir_df:

  Optional CIR results table. Must contain `idx_M0` and `is_valid_c13`.
  If supplied, enables the N/O/S gate.

- morph_df:

  Optional morphology table (typically the same object passed in
  `result`). Must contain `idx_M0`, `idx_cand`, `iso_type`. If supplied,
  enables exclusion of C13_M2 pairs.

- require_c13_for:

  Character vector of elements for which CIR-validated C13_M1 is
  required. Default: `c("N","O","S")`.

- exclude_c13_m2:

  Logical. If `TRUE`, excludes pairs already detected as C13_M2 for
  `idx_M0` values with CIR-valid C13_M1.

- return_debug:

  Logical. If `TRUE`, returns a list with final output and internal
  top-n candidates.

## Value

By default, a single `data.frame` / tibble with one row per evaluated
pair:

- `idx_M0`, `mz_M0`, `idx_cand`, `mz_iso`

- `iso_type`, `element`, `k`, `delta`

- `R_obs`, `n_hat`, `R_theo_hat`, `eips_rel_err`

- `is_valid_eips`, `score_eips`

- `has_formula_support`, `adduct_name`, `neutral_mass`, `mz_mono_db`,
  `n_el_db`, `neutral_err_ppm`

If `return_debug = TRUE`, returns a list with:

- `eips_validation`

- `top_n_candidates`
