# Raw Influenza A wastewater data from Zurich

Wastewater measurements of Influenza A RNA from the WWTP Werdhölzli in
Zürich, Switzerland (September 2022–April 2023). Data are provided by
EAWAG.

## Usage

``` r
Influenza_A_Zurich
```

## Format

A list with the following elements:

- measurements:

  A `data.table` with columns `date`, `concentration` (gc/mL),
  `replicate_id`, `dilution`, `total_droplets`, and `positive_droplets`.

- flows:

  A `data.table` with columns `date` and `flow` (mL/day).

- units:

  A list with concentration and flow unit strings.
