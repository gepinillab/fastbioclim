# `fastbioclim`: Scalable Derivation of Climate Variables

<!-- badges: start -->
[![R-CMD-check](https://github.com/gepinillab/fastbioclim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gepinillab/fastbioclim/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`fastbioclim` is an R package for efficiently deriving standard bioclimatic and custom summary variables from large-scale climate raster data. It is designed to overcome the memory limitations of traditional approaches by intelligently switching between processing backends.

## Overview

Working with large climate datasets often presents a major challenge: the data is too large to fit into memory. `fastbioclim` solves this problem by providing a powerful and unified interface with a dual-backend architecture:

1.  **In-Memory (`"terra"`)**: For datasets that fit comfortably in RAM, `fastbioclim` uses the highly optimized `terra` package for maximum speed.
2.  **In-Disk (`"tiled"`)**: For massive datasets that exceed available RAM, it automatically switches to a memory-safe tiled workflow. This backend uses `exactextractr` and `Rfast` to process the data chunk by chunk, ensuring that even continent-scale analyses can run on a standard computer.

The core of the package is a set of "smart" wrapper functions—`derive_bioclim()` and `derive_statistics()`—that automatically select the best backend, providing a seamless experience for the user.

## Key Features

*   **Smart `auto` Method**: Automatically detects data size and chooses the optimal processing workflow (in-memory vs. tiled).
*   **Scalability**: Robustly handles raster datasets of any size, from small study areas to global extents.
*   **Comprehensive Bioclimatic Variables**: Calculates the full set of 35 standard bioclimatic variables, including those based on temperature, precipitation, solar radiation, and moisture.
*   **Custom Statistics**: Provides a flexible `derive_statistics()` function to compute custom summaries (mean, max, min, standard deviation, etc.) for any climate variable.
*   **Flexible Time Periods**: Easily define custom periods (e.g., weeks, semesters) instead of being limited to standard quarters.
*   **Advanced Control**: Supports static index rasters for advanced time-series analysis and defining specific seasonal periods.

## Installation

You can install the development version of `fastbioclim` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("gepinillab/fastbioclim")
```

## Usage

The package provides two primary, user-facing functions:

*   **`derive_bioclim()`**: For calculating the standard set of 19-35 bioclimatic variables.
*   **`derive_statistics()`**: For calculating custom summary statistics on any variable.

### Example: A Quick Demonstration

This example demonstrates the core functionality using simple, self-contained mock data.

```r
library(fastbioclim)
library(terra)

# 1. Create mock monthly climate data (12 layers) for a small area
set.seed(123) # for reproducibility
r <- rast(nrows = 10, ncols = 10, nlyr = 12, 
          xmin = -79, xmax = -78, ymin = 13, ymax = 14, 
          crs = "EPSG:4326")

tmin_rast <- setValues(r, rnorm(ncell(r) * nlyr(r), mean = 15, sd = 2))
tmax_rast <- tmin_rast + rnorm(ncell(r) * nlyr(r), mean = 10, sd = 1)
prcp_rast <- setValues(r, rpois(ncell(r) * nlyr(r), lambda = 150))

# 2. Derive standard bioclimatic variables
# The function will automatically use the fast "terra" method for this small dataset
output_dir_bioclim <- file.path(tempdir(), "bioclim_output")

bioclim_vars <- derive_bioclim(
  bios = 1:19,
  tmin = tmin_rast,
  tmax = tmax_rast,
  prcp = prcp_rast,
  output_dir = output_dir_bioclim,
  overwrite = TRUE
)

print(bioclim_vars)
plot(bioclim_vars[[c("bio01", "bio12")]])
```
<img src="man/figures/README-bioclim_plot.png" width="80%" />

```r
# 3. Derive custom summary statistics for a different variable (e.g., wind speed)
wind_rast <- setValues(r, rnorm(ncell(r) * nlyr(r), mean = 5, sd = 2))
output_dir_custom <- file.path(tempdir(), "custom_output")

custom_stats <- derive_statistics(
  variable = wind_rast,
  stats = c("mean", "max", "stdev"),
  prefix_variable = "wind",
  output_dir = output_dir_custom,
  overwrite = TRUE
)

print(custom_stats)
plot(custom_stats)
```
<img src="man/figures/README-custom_plot.png" width="80%" />

### Handling Large-Scale Data (The Tiled Workflow)

The real power of `fastbioclim` shines with large datasets. The `method = "auto"` setting in `derive_bioclim()` and `derive_statistics()` handles this automatically.

When the wrapper function detects that the input rasters are too large to fit in memory, it seamlessly switches to the tiled workflow.

**Important Requirement:** For the tiled workflow to function, your input `SpatRaster` objects must be pointing to files on disk, not held entirely in memory.

```r
# Conceptual example for large, file-based rasters
tmin_on_disk <- rast("path/to/your/large_tmin.tif")
tmax_on_disk <- rast("path/to/your/large_tmax.tif")
prcp_on_disk <- rast("path/to/your/large_prcp.tif")

# The call is identical. `derive_bioclim` will detect the large file size
# and automatically use the memory-safe tiled method.
large_scale_vars <- derive_bioclim(
  bios = 1:19,
  tmin = tmin_on_disk,
  tmax = tmax_on_disk,
  prcp = prcp_on_disk,
  output_dir = "./final_bioclim_rasters"
)
```

## Under Active Development

This R package is currently under active development. While it is functional, it may contain bugs or undergo changes to the API.

Contributions, bug reports, and feature requests are highly encouraged. Please [open an issue](https://github.com/gepinillab/fastbioclim/issues) on our GitHub repository to provide feedback.
