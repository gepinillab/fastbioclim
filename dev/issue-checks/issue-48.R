# Check for issue #48: bioclim_terra() forwards every static index in `...` to
# each bio function and stops with "unused argument".
#
# Run with: Rscript dev/issue-checks/issue-48.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #48: bioclim_terra() passes every static index through ... and crashes with 'unused argument'\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(48)
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
write_index <- function(value, name) {
  r <- rast(template)
  values(r) <- value
  f <- file.path(dir_in, paste0(name, ".tif"))
  writeRaster(r, f, datatype = "INT2S", overwrite = TRUE)
  rast(f)
}
tmin_v <- matrix(rnorm(ncell(template) * 12, 10, 5), ncell(template))
tmax_v <- tmin_v + 10
tmin <- write_units(tmin_v, "tmin")
tmax <- write_units(tmax_v, "tmax")
# Constant indices: warmest unit 3, coldest unit 7, warmest period 2 (units 2:4).
warmest_unit <- write_index(3, "warmest_unit")
coldest_unit <- write_index(7, "coldest_unit")
warmest_period <- write_index(2, "warmest_period")

run <- function(method) {
  tryCatch({
    res <- suppressWarnings(suppressMessages(
      derive_bioclim(bios = c(5, 6, 10), tmin = tmin, tmax = tmax,
                     warmest_unit = warmest_unit, coldest_unit = coldest_unit,
                     warmest_period = warmest_period, method = method,
                     output_dir = tempfile(), verbose = FALSE)
    ))
    values(res)
  }, error = function(e) conditionMessage(e), finally = suppressMessages(sf::sf_use_s2(old_s2)))
}
terra_v <- run("terra")
tiled_v <- run("tiled")

expected <- cbind(bio05 = tmax_v[, 3], bio06 = tmin_v[, 7],
                  bio10 = rowMeans((tmin_v[, 2:4] + tmax_v[, 2:4]) / 2))
matches <- function(v) is.matrix(v) && isTRUE(all.equal(unname(v[, c("bio05", "bio06", "bio10")]),
                                                         unname(expected), tolerance = 1e-5))
describe <- function(v) if (is.matrix(v)) paste("values, cell 1 =", paste(signif(v[1, ], 6), collapse = ", ")) else paste("error:", v)

cat("EXPECTED: bio05, bio06 and bio10 from the static indices on both backends; cell 1 =",
    paste(signif(expected[1, ], 6), collapse = ", "), "\n")
cat("GOT:      terra", describe(terra_v), "\n")
cat("          tiled", describe(tiled_v), "\n")
cat("RESULT:", if (matches(terra_v) && matches(tiled_v)) "PASS" else "FAIL", "\n")
