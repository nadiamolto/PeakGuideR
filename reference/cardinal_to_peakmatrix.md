# Convert a Cardinal MSImagingExperiment object to a peak matrix

Converts a supported Cardinal MSI object into a PeakGuideR peak matrix
with an rMSIprocPeakMatrix-like structure.

## Usage

``` r
cardinal_to_peakmatrix(
  x,
  value = NULL,
  dataset_name = NULL,
  snr = NULL,
  area = NULL
)
```

## Arguments

- x:

  A Cardinal MSImagingExperiment object.

- value:

  Name of the imageData layer to extract. If NULL, "intensity" is used.

- dataset_name:

  Optional dataset name.

- snr:

  Optional SNR matrix.

- area:

  Optional area matrix.

## Value

An object with rMSIprocPeakMatrix-like structure. In addition to the
usual peak matrix fields, it contains `image_id`, a factor with one
entry per pixel giving the run/image the pixel comes from
(`Cardinal::run(x)`), or `NULL` if it could not be obtained. It is used
to plot each image separately (see
[`plot_ion_image()`](https://nadiamolto.github.io/PeakGuideR/reference/plot_ion_image.md))
and by `run_peakguider_workflow(multi_image = TRUE)`. The other fields,
including `numPixels` and `names`, are unchanged.

## Examples

``` r
if (FALSE) { # \dontrun{
# `cardinal_msi_data` is a Cardinal MSImagingExperiment object.
pkm <- cardinal_to_peakmatrix(cardinal_msi_data)

res <- run_peakguider_workflow(
  pkm = pkm,
  ion_mode = "pos",
  matrix = "HCCA"
)
} # }
```
