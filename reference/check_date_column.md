# Check and convert date column to Date type

Internal helper function to validate that a date column is of type Date
and convert it if necessary. Handles POSIXct/POSIXlt conversion and
parsing of character/factor columns.

## Usage

``` r
check_date_column(data, date_col_name)
```

## Arguments

- data:

  A data.table containing the date column.

- date_col_name:

  The name of the date column (for error messages).

## Value

The data.table with the date column validated and converted to Date.
