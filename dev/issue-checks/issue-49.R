# Check for issue #49: derive_statistics(method = "terra") fails with any static
# index because the geometry check compares layer counts.
#
# Run with: Rscript dev/issue-checks/issue-49.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #49: derive_statistics(method = \"terra\") fails with any static index: compareGeom number of layers does not match\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(49)
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
prcp_v <- matrix(runif(ncell(template) * 12, 0, 100), ncell(template))
tmin <- write_units(tmin_v, "tmin")
prcp <- write_units(prcp_v, "prcp")
max_unit <- write_index(4, "max_unit")
max_period <- write_index(2, "max_period")
max_interactive <- write_index(5, "max_interactive")

run <- function(method) {
  tryCatch({
    res <- suppressWarnings(suppressMessages(
      derive_statistics(tmin, stats = c("max", "max_period"), inter_variable = prcp,
                        inter_stats = "max_inter", max_unit = max_unit,
                        max_period = max_period, max_interactive = max_interactive,
                        method = method, output_prefix = "tmin", output_dir = tempfile(),
                        verbose = FALSE)
    ))
    values(res)
  }, error = function(e) conditionMessage(e), finally = suppressMessages(sf::sf_use_s2(old_s2)))
}
terra_v <- run("terra")
tiled_v <- run("tiled")

# Selected unit, and the mean over the 3-unit circular period starting at the index.
# Layers are compared by name because the backends order them differently.
expected <- cbind(tmin_max = tmin_v[, 4],
                  tmin_max_period = rowMeans(tmin_v[, 2:4]),
                  tmin_inter_high = rowMeans(tmin_v[, 5:7]))
matches <- function(v) is.matrix(v) && setequal(colnames(v), colnames(expected)) &&
  isTRUE(all.equal(unname(v[, colnames(expected)]), unname(expected), tolerance = 1e-5))
describe <- function(v) if (is.matrix(v)) paste("values, cell 1 =", paste(colnames(v), signif(v[1, ], 6), sep = "=", collapse = ", ")) else paste("error:", v)

cat("EXPECTED: values from the static indices on both backends; cell 1 =",
    paste(colnames(expected), signif(expected[1, ], 6), sep = "=", collapse = ", "), "\n")
cat("GOT:      terra", describe(terra_v), "\n")
cat("          tiled", describe(tiled_v), "\n")
cat("RESULT:", if (matches(terra_v) && matches(tiled_v)) "PASS" else "FAIL", "\n")
