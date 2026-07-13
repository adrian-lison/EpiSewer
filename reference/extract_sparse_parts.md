# Extract the sparse representation of a dgRMatrix

This function extracts the non-zero values and their corresponding row
and column indices from a dgRMatrix, which is a sparse matrix format in
R. It returns a list containing the non-zero values, their column
indices (1-based), and the row pointers (1-based).

## Usage

``` r
extract_sparse_parts(A)
```

## Arguments

- A:

  A dgRMatrix object representing a sparse matrix.

## Value

A list with three components:

- `w`: A numeric vector of non-zero values from the matrix.

- `v`: An integer vector of column indices corresponding to the non-zero
  values (1-based).

- `u`: An integer vector of row pointers indicating the start of each
  row in the non-zero values vector (1-based).
