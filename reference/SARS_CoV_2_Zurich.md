# Raw SARS-CoV-2 wastewater data from Zurich

Wastewater measurements of SARS-CoV-2 RNA and case counts from the WWTP
Werdhölzli in Zürich, Switzerland (January–April 2022). Data are
provided by EAWAG under the CC BY 4.0 license.

## Usage

``` r
SARS_CoV_2_Zurich
```

## Format

A list with the following elements:

- measurements:

  A `data.table` with columns `date` and `concentration` (gc/mL).

- flows:

  A `data.table` with columns `date` and `flow` (mL/day).

- cases:

  A `data.table` with columns `date` and `cases` (cases/day).

- units:

  A list with concentration and flow unit strings.

## Source

<https://sensors-eawag.ch/sars/>
