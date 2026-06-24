# CIR theoretical lookup table (100-1200 Da)

Precomputed M+1/M0 ratios for monoisotopic masses from reference
formulas.

## Usage

``` r
cir_table
```

## Format

### cir_table

A data frame with 1,101 rows and 2 columns:

- mz:

  monoisotopic m/z (100-1200 Da, 1 Da steps)

- R_theo:

  theoretical CIR ratio M/M+1

## Source

ChEBI/HMDB reference databases

## Details

Generated via GAM (k=20, gamma=1.2) fitted over carbon counts from ~184k
ChEBI/HMDB formulas using natural 13C abundance (p=0.0107).

Used by
[`cir_ratio()`](https://nadiamolto.github.io/PeakGuideR/reference/cir_ratio.md)
and
[`cir_score()`](https://nadiamolto.github.io/PeakGuideR/reference/cir_score.md)
for C13 isotope validation.
