# bio15_terra: Precipitation Seasonality (CV)

Calculates coefficient of variation in precipitation across units.

## Usage

``` r
bio15_terra(prcp)
```

## Arguments

- prcp:

  Matrix containing precipitation values for each unit.

## Value

spatRaster with "bio15".

## Note

The "1 +" is to avoid strange CVs for areas where mean rainfaill is \<
1)
