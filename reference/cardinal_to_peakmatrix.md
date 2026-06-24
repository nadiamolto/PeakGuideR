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

An object with rMSIprocPeakMatrix-like structure.
