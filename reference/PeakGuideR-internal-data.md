# Internal reference tables

Internal lookup tables used by the CIR and EIPS scoring functions and by
the mass defect flag.

## Details

- `cir_table`: carbon isotope-ratio reference table.

- `eips_table`: formula-level elemental isotope lookup table.

- `eips_n_table`: theoretical element-count isotope lookup table.

- `mass_defect_table`: limits of the mass defect (minimum and maximum)
  observed among HMDB and ChEBI compounds, one row per integer mass from
  55 to 1300 Da, used by the mass defect flag (built by
  `data-raw/mass_defect_table.R`).

These objects are internal implementation data and are not part of the
public PeakGuideR API.
