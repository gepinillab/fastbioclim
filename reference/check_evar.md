# Standardize Environmental Variable Input Format

Converts data frames, matrices, or vectors into a consistent matrix
format with preserved column/row names for downstream processing.

## Usage

``` r
check_evar(evar)
```

## Arguments

- evar:

  Input environmental variable data. Can be: - data.frame: Converted to
  matrix - vector: Transposed to row matrix - matrix: Preserved with
  names checked

## Value

Matrix with column names preserved. If input was a vector, returns 1-row
matrix with vector names as column names.
