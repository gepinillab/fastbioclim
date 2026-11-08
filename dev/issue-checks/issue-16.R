# Check for issue #16: a static period index cannot replace the variable that
# would otherwise select the period, on either backend.
#
# Run with: Rscript dev/issue-checks/issue-16.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #16: Static period indices can't replace the data they stand in for, in either backend\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(16)
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
prcp_v <- matrix(runif(ncell(template) * 12, 0, 100), ncell(template))
tmin <- write_units(tmin_v, "tmin")
tmax <- write_units(tmax_v, "tmax")
prcp <- write_units(prcp_v, "prcp")
# Constant period 2 = units 2:4 (period_length 3).
period <- write_index(2, "period")

# Case A: bio08 (mean temperature of the wettest period) with wettest_period and no prcp.
# Case B: bio18 (precipitation of the warmest period) with warmest_period and no tmin/tmax.
run <- function(method, case) {
  tryCatch({
    res <- suppressWarnings(suppressMessages(
      if (case == "A") {
        derive_bioclim(bios = 8, tmin = tmin, tmax = tmax, wettest_period = period,
                       method = method, output_dir = tempfile(), verbose = FALSE)
      } else {
        derive_bioclim(bios = 18, prcp = prcp, warmest_period = period,
                       method = method, output_dir = tempfile(), verbose = FALSE)
      }
    ))
    values(res)[, 1]
  }, error = function(e) conditionMessage(e), finally = suppressMessages(sf::sf_use_s2(old_s2)))
}
expected <- list(A = rowMeans((tmin_v[, 2:4] + tmax_v[, 2:4]) / 2),
                 B = rowSums(prcp_v[, 2:4]))
describe <- function(v) if (is.numeric(v)) paste("cell 1 =", signif(v[1], 6)) else paste("error:", v)
ok <- TRUE
for (case in c("A", "B")) {
  for (method in c("terra", "tiled")) {
    got <- run(method, case)
    pass <- is.numeric(got) && isTRUE(all.equal(got, expected[[case]], tolerance = 1e-5))
    ok <- ok && pass
    cat(sprintf("case %s, %-5s: EXPECTED cell 1 = %s | GOT %s\n", case, method,
                signif(expected[[case]][1], 6), describe(got)))
  }
}
cat("EXPECTED: bio08 and bio18 computed from the static period index on both backends, without the variable it replaces\n")
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
