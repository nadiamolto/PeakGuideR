# Gaussian plausibility score from a non-negative error

This script is designed to make the scores compatible.

## Usage

``` r
gaussian_score(err, tol)
```

## Arguments

- err:

  Numeric vector of non-negative errors.

- tol:

  Numeric tolerance.controlling the width of the Gaussian penalty. With
  this formulation, the score is approximately 0.61 when `err = tol`.

## Value

Numeric vector in range 0 to 1.
