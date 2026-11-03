# Check for issue #28: an out-of-range static period index returns the
# min_idx/max_idx columns as data in bios_fast.
#
# Run with: Rscript dev/issue-checks/issue-28.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #28: Out-of-range static period index returns min_idx/max_idx values as data in bios_fast\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(28)
dir_in <- tempfile("in_")
dir.create(dir_in)
template <- rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3, crs = "EPSG:4326")
write_units <- function(vals, prefix) {
  f <- file.path(dir_in, sprintf("%s_%02d.tif", prefix, 1:12))
  for (i in 1:12) {
    r <- rast(template)
    values(r) <- vals[, i]
    writeRaster(r, f[i], datatype = "FLT8S", overwrite = TRUE)
  }
  rast(f)
}
tmin_v <- matrix(rnorm(ncell(template) * 12, 10, 5), ncell(template))
tmin <- write_units(tmin_v, "tmin")
tmax <- write_units(tmin_v + 10, "tmax")
prcp <- write_units(matrix(runif(ncell(template) * 12, 0, 100), ncell(template)), "prcp")

# Static warmest period 11 everywhere. With circular = FALSE and 12 units there are
# only 10 periods, so 11 is out of range.
wp <- rast(template)
values(wp) <- 11
wp_file <- file.path(dir_in, "warmest_period.tif")
writeRaster(wp, wp_file, datatype = "INT2S", overwrite = TRUE)
wp <- rast(wp_file)

run <- function(method) {
  res <- suppressWarnings(suppressMessages(
    derive_bioclim(bios = c(10, 18), tmin = tmin, tmax = tmax, prcp = prcp, circular = FALSE,
                   warmest_period = wp, method = method, output_dir = tempfile(),
                   verbose = FALSE)
  ))
  suppressMessages(sf::sf_use_s2(old_s2))
  res
}
terra_v <- values(run("terra"))
tiled_v <- values(run("tiled"))

cat("EXPECTED: bio10 and bio18 NA in every cell on both backends\n")
cat("GOT:      cell 1 terra bio10 =", terra_v[1, "bio10"], "bio18 =", terra_v[1, "bio18"],
    "| tiled bio10 =", tiled_v[1, "bio10"], "bio18 =", tiled_v[1, "bio18"], "\n")
ok <- all(is.na(terra_v)) && all(is.na(tiled_v))
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
