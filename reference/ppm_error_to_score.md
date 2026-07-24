# Convert a ppm mass error into a 0-1 plausibility score

Converts a candidate mass error (in ppm) into a `mass_error_score` in
the range 0 to 1, using the same Gaussian penalty as
[`gaussian_score()`](https://nadiamolto.github.io/PeakGuideR/reference/gaussian_score.md).

## Usage

``` r
ppm_error_to_score(candidate_ppm_error, ppm_tol)
```

## Arguments

- candidate_ppm_error:

  Numeric vector of non-negative ppm errors.

- ppm_tol:

  Numeric tolerance (ppm) controlling the width of the penalty.

## Value

Numeric vector in range 0 to 1.
