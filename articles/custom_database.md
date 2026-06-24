# Using_custom_databases

## Using custom databases

PeakGuideR can be run with user-provided compound databases and
standard-adduct libraries. This allows users to analyse their data with
in-house databases, project-specific compound lists, or databases
generated from other sources.

### Custom compound database

A custom compound database can be supplied through the `compound_db`
argument. The table must be a `data.frame` and must contain, at minimum,
a numeric `MonoisotopicMass` column. This column is used to match
inferred neutral masses against candidate compounds.

Recommended columns are shown below.

| Column | Required | Description |
|----|---:|----|
| `MonoisotopicMass` | yes | Neutral monoisotopic mass used for ppm matching. |
| `Source` | no | Name of the database or source. |
| `DB_ID` | no | Compound identifier in the source database. |
| `Name` | no | Compound name. |
| `MolecularFormula` | no | Molecular formula. |
| `StdInChI` | no | Standard InChI. |
| `StdInChIKey` | no | Standard InChIKey. |
| `SMILES` | no | SMILES representation. |
| `Kegg` | no | KEGG identifier, if available. |

Column names are case-sensitive. For example, the mass column must be
named `MonoisotopicMass`.

``` r

my_compound_db <- data.frame(
  Source = c("custom", "custom"),
  DB_ID = c("C001", "C002"),
  Name = c("Citric acid", "L-carnitine"),
  MolecularFormula = c("C6H8O7", "C7H15NO3"),
  MonoisotopicMass = c(192.0270, 161.1052),
  stringsAsFactors = FALSE
)

res <- run_peakguider_workflow(
  pkm = pkm,
  ion_mode = "pos",
  matrix = "HCCA",
  compound_db = my_compound_db
)
```

### Custom standard-adduct library

A custom standard-adduct library can be supplied through the
`standards_db` argument. This table is optional. It is used only to add
standard-support evidence to candidate compound matches.

The table must be a `data.frame` and must contain, at minimum, the
columns `adduct` and `POLARITY`.

| Column | Required | Description |
|----|---:|----|
| `adduct` | yes | Adduct name, for example `[M+H]+`, `[M+Na]+` or `[M+K]+`. |
| `POLARITY` | yes | Ion mode associated with the adduct, usually `pos` or `neg`. |
| `COMPOUND_ID` | no | Internal or source-specific compound identifier. |
| `Master_List_NAME` | no | Standard compound name. |
| `HMDB_clean` | no | HMDB identifier, if available. |
| `ChEBI` | no | ChEBI identifier, if available. |
| `InCHIKey` | no | InChIKey used for compound matching. |
| `SMILES` | no | SMILES used for compound matching. |
| `MOLECULAR_FORMULA` | no | Molecular formula used for compound matching. |

``` r

my_standards_db <- data.frame(
  Master_List_NAME = c("Citric acid", "Citric acid"),
  MOLECULAR_FORMULA = c("C6H8O7", "C6H8O7"),
  adduct = c("[M+Na]+", "[M+K]+"),
  POLARITY = c("pos", "pos"),
  stringsAsFactors = FALSE
)

res <- run_peakguider_workflow(
  pkm = pkm,
  ion_mode = "pos",
  matrix = "HCCA",
  compound_db = my_compound_db,
  standards_db = my_standards_db
)
```

When using custom databases, users are responsible for ensuring that the
database content and identifiers are appropriate for their analysis and
comply with the licenses of the original data sources.
