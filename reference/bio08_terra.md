# bio08_terra: Mean Temperature of Wettest Period

Calculates mean temperature of the period with the highest precipitation
sum.

## Usage

``` r
bio08_terra(tmp, wettest_period)
```

## Arguments

- tmp:

  spatRaster of temperature period sums.

- wettest_period:

  spatRaster indicating the index (1-based) of the wettest period for
  each cell.

## Value

spatRaster with "bio08".
