# access_to_databases

## Compound databases and standard-adduct libraries

PeakGuideR can annotate inferred neutral masses using compound databases
and, optionally, standard-adduct libraries.

By default, the package includes small example databases so that the
workflow can be run without downloading any external files. These
example databases are intended for testing, examples and vignette
execution. They are not intended to be comprehensive reference
databases.

The example compound database is loaded automatically when
`compound_db = NULL`. It contains a small set of compounds with neutral
monoisotopic masses and basic metadata.

The example standard-adduct library is loaded automatically when
`standards_db = NULL` and `matrix = "HCCA"`. It contains a set of
matrix-specific standard-adduct entries used to add standard-support
evidence to candidate matches.

For larger analyses, users can either download the full PeakGuideR
reference databases or provide their own custom databases.

Full reference databases can be downloaded separately from Zenodo:

<https://doi.org/10.5281/zenodo.20705395>

These full databases are distributed separately from the R package
because they may contain metadata derived from third-party resources.
Users are responsible for complying with the licenses of the original
data sources.

After downloading the files, replace the paths below with the local
paths where the files were saved.

``` r

dbs <- load_peakguider_databases(
  compound_db_path = "~/Downloads/compound_mass_database.rds",
  standards_db_path = "~/Downloads/standards_adduct_library.rds"
)

compound_db <- dbs$compound_db
standards_db <- dbs$standards_db
```

The downloaded databases can then be passed to the workflow:

``` r

res <- run_peakguider_workflow(
  pkm = pkm,
  ion_mode = "pos",
  matrix = "HCCA",
  compound_db = compound_db,
  standards_db = standards_db
)
```
