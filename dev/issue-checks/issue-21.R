# Check for issue #21: NA handling differs between terra and tiled backends.
#
# Expected policy (from the issue): an output is NA when any month of the data it
# is computed from is NA, and both backends give the same result.
#
# Run with: Rscript dev/issue-checks/issue-21.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #21: NA handling differs between terra and tiled backends\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(1)
dir_in <- tempfile("in_")
dir.create(dir_in)
template <- rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:4326")
n <- ncell(template)
season <- sin(2 * pi * (0:11) / 12)
tmin <- t(sapply(seq_len(n), function(i) 5 + 8 * season + rnorm(12)))
tmax <- tmin + 10
prcp <- t(sapply(seq_len(n), function(i) 80 + 60 * cos(2 * pi * (0:11) / 12 + i) + runif(12, 0, 20)))
prcp[1, 12] <- NA   # cell 1: December precipitation missing
tmin[2, 1] <- NA    # cell 2: January tmin missing (first layer the tiled backend reads)

write_units <- function(vals, prefix) {
  f <- file.path(dir_in, sprintf("%s_%02d.tif", prefix, 1:12))
  for (i in 1:12) {
    r <- rast(template)
    values(r) <- vals[, i]
    writeRaster(r, f[i], datatype = "FLT8S", overwrite = TRUE)
  }
  rast(f)
}
tmin_r <- write_units(tmin, "tmin")
tmax_r <- write_units(tmax, "tmax")
prcp_r <- write_units(prcp, "prcp")

bios <- c(1, 5, 12, 16, 18)
run <- function(method) {
  res <- tryCatch(
    suppressWarnings(suppressMessages(
      derive_bioclim(bios = bios, tmin = tmin_r, tmax = tmax_r, prcp = prcp_r,
                     method = method, tile_degrees = 2, output_dir = tempfile(),
                     verbose = FALSE)
    )),
    error = function(e) e
  )
  suppressMessages(sf::sf_use_s2(old_s2))
  res
}
terra_res <- run("terra")
tiled_res <- run("tiled")

# Which outputs should be NA in cells 1 and 2:
# cell 1 (prcp NA): bio12, bio16 and bio18 use prcp; bio01 and bio05 do not.
# cell 2 (tmin NA): bio01 and bio18 use tavg = (tmin + tmax) / 2; the others do not.
expected_na <- rbind(
  cell1 = c(bio01 = FALSE, bio05 = FALSE, bio12 = TRUE, bio16 = TRUE, bio18 = TRUE),
  cell2 = c(bio01 = TRUE, bio05 = FALSE, bio12 = FALSE, bio16 = FALSE, bio18 = TRUE)
)
fmt <- function(m) paste(apply(m, 1, function(r) paste(names(r)[r], collapse = "+")),
                         collapse = " ; ")

if (inherits(terra_res, "error") || inherits(tiled_res, "error")) {
  cat("EXPECTED: NA cells", fmt(expected_na), "on both backends\n")
  cat("GOT:      error:", conditionMessage(if (inherits(tiled_res, "error")) tiled_res else terra_res), "\n")
  cat("RESULT: FAIL\n")
} else {
  tv <- values(terra_res)[1:2, , drop = FALSE]
  dv <- values(tiled_res)[1:2, , drop = FALSE]
  rownames(tv) <- rownames(dv) <- c("cell1", "cell2")
  terra_na <- is.na(tv)
  tiled_na <- is.na(dv)
  all_t <- values(terra_res)
  all_d <- values(tiled_res)
  agree <- identical(is.na(all_t), is.na(all_d)) &&
    isTRUE(all.equal(all_t[!is.na(all_t)], all_d[!is.na(all_d)], tolerance = 1e-5))
  cat("EXPECTED: NA in", fmt(expected_na), "| backends agree everywhere\n")
  cat("GOT:      terra NA in", fmt(terra_na), "| tiled NA in", fmt(tiled_na),
      "| backends agree:", agree, "\n")
  cat("          cell 1 terra:", paste(colnames(tv), signif(tv[1, ], 5), sep = "=", collapse = " "), "\n")
  cat("          cell 1 tiled:", paste(colnames(dv), signif(dv[1, ], 5), sep = "=", collapse = " "), "\n")
  ok <- agree && identical(unname(terra_na), unname(expected_na)) &&
    identical(unname(tiled_na), unname(expected_na))
  cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
}
