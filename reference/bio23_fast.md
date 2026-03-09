# bio23_fast: Radiation Seasonality (CV)

Calculates the Coefficient of Variation (CV) of solar radiation. The
formula used is: \`(StandardDeviation / (Mean + 1)) \* 100\`. (The "+1"
is added to the mean to avoid division by zero).

## Usage

``` r
bio23_fast(srad, n_units, cell)
```

## Arguments

- srad:

  A numeric \*\*matrix\*\* of solar radiation values. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units (e.g., 12 months).

- n_units:

  A single \*\*integer\*\* representing the number of temporal units
  (e.g., 12).

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`srad\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio23" (Solar Radiation Seasonality)
and "cell".
