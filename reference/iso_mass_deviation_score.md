# Mass deviation scoring for isotopic candidate pairs

Adds a continuous mass-accuracy score to isotopic pairs detected by
[`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md).
This does NOT search new candidates: it only scores the already-detected
pairs using their mass error (`mass_err_da` / `mass_err_ppm`).

This is useful when the morphology step uses a hard tolerance window
(pass/fail), but you still want a smooth penalty to rank candidates by
how close they are to the theoretical delta.

## Usage

``` r
iso_mass_deviation_score(
  pairs,
  tol_ppm = 5,
  tol_da = NULL,
  use = c("ppm", "da"),
  kernel = c("gaussian", "linear"),
  combine = TRUE,
  w_morph = 0.8
)
```

## Arguments

- pairs:

  Data.frame output of
  [`iso_morphology_candidates()`](https://nadiamolto.github.io/PeakGuideR/reference/iso_morphology_candidates.md).
  Must contain `score_final` and either `mass_err_ppm` or `mass_err_da`.

- tol_ppm:

  Numeric. PPM tolerance used as the reference window. Default `5`.
  (Typically the same `tol_ppm` used in `iso_morphology_candidates`.)

- tol_da:

  Numeric. Da tolerance used if scoring in Da (only used if `use="da"`).

- use:

  Character: `"ppm"` (default) or `"da"` to decide which error column to
  use.

- kernel:

  Character: `"gaussian"` (default) or `"linear"`.

- combine:

  Logical. If `TRUE`, compute `score_combined`. Default `TRUE`.

- w_morph:

  Numeric in range 0 to 1. Weight for morphology when combining. Default
  `0.8`.

## Value

The input `pairs` with additional columns:

- `mass_dev_err`: absolute mass error in chosen units (ppm or Da)

- `mass_dev_score`: continuous score in the range from 0 to 1.

- `score_combined`: optional combined score

## Details

The function expects `pairs` to contain `mass_err_ppm` (preferred) or
`mass_err_da`. The score is computed as either:

- **Gaussian (default)**: \$\$score = exp(-(err^2)/(2\*sigma^2))\$\$
  where `sigma = tol/2`.

- **Linear**: \$\$score = pmax(0, 1 - \|err\|/tol)\$\$

Optionally, a combined score can be produced: \$\$score\\combined =
w\\morph \* score\\final + (1-w\\morph) \* mass\\dev\\score\$\$

## Examples

``` r
if (FALSE) { # \dontrun{
pairs <- iso_morphology_candidates(pm, tol_ppm=5)
pairs2 <- iso_mass_deviation_score(pairs, tol_ppm=5, kernel="gaussian", w_morph=0.8)
head(pairs2)
} # }
```
