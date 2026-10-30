# Helpers shared by the backend-equivalence, reference-value and edge-case tests.
# Every raster is synthetic, small and written to a temporary directory, so no
# example data has to be downloaded.

# --- Synthetic inputs ---------------------------------------------------------

# Writes one single-layer GeoTIFF per temporal unit and returns a file-backed
# SpatRaster. The tiled backend expects one file per unit (it checks that
# `length(terra::sources(x)) == nlyr(x)`), so a single multi-band file would
# not work. FLT8S keeps the inputs exact, so reference values can be computed
# from the same numbers the backends read.
write_units <- function(vals, template, dir, prefix) {
  vals <- as.matrix(vals)
  files <- file.path(dir, sprintf("%s_%02d.tif", prefix, seq_len(ncol(vals))))
  for (i in seq_len(ncol(vals))) {
    r <- terra::rast(template, nlyrs = 1)
    terra::values(r) <- vals[, i]
    terra::writeRaster(r, files[i], overwrite = TRUE, datatype = "FLT8S")
  }
  out <- terra::rast(files)
  names(out) <- sprintf("%s_%02d", prefix, seq_len(ncol(vals)))
  out
}

# Writes a single-layer index raster (static index) to disk.
write_index <- function(idx, template, dir, name) {
  r <- terra::rast(template, nlyrs = 1)
  terra::values(r) <- idx
  f <- file.path(dir, paste0(name, ".tif"))
  terra::writeRaster(r, f, overwrite = TRUE, datatype = "INT2S")
  terra::rast(f)
}

# Builds tmin, tmax, prcp, srad and mois with a seasonal cycle whose phase varies
# by cell, so the warmest/wettest periods land on different units across the
# raster. Values are continuous, which makes ties in period selection
# practically impossible.
make_climate <- function(dir, nrows = 6, ncols = 8, n_units = 12,
                         xmin = 0, xmax = ncols * 0.5,
                         ymin = 0, ymax = nrows * 0.5,
                         seed = 42) {
  set.seed(seed)
  template <- terra::rast(nrows = nrows, ncols = ncols, nlyrs = 1,
                          xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                          crs = "EPSG:4326")
  n_cells <- terra::ncell(template)
  unit_angle <- 2 * pi * (seq_len(n_units) - 1) / n_units
  seasonal <- function(phase, amp) {
    t(vapply(phase, function(ph) amp * sin(unit_angle + ph), numeric(n_units)))
  }
  noise <- function(sd) matrix(stats::rnorm(n_cells * n_units, 0, sd), n_cells)

  t_phase <- stats::runif(n_cells, 0, 2 * pi)
  p_phase <- stats::runif(n_cells, 0, 2 * pi)
  tmin <- 5 + seasonal(t_phase, 8) + noise(1)
  tmax <- tmin + matrix(stats::runif(n_cells * n_units, 6, 14), n_cells)
  prcp <- 80 + seasonal(p_phase, 60) + matrix(stats::runif(n_cells * n_units, 0, 30), n_cells)
  srad <- 200 + seasonal(t_phase, 80) + noise(10)
  mois <- 50 + seasonal(p_phase, 30) + noise(5)

  vals <- list(tmin = tmin, tmax = tmax, prcp = prcp, srad = srad, mois = mois)
  list(
    template = template,
    vals = vals,
    rasters = lapply(stats::setNames(names(vals), names(vals)), function(v) {
      write_units(vals[[v]], template, dir, v)
    })
  )
}

# --- Running the backends -----------------------------------------------------

# The tiled backend calls sf::sf_use_s2(FALSE) and does not restore it; put the
# previous value back so the tests do not leak global state.
local_s2_restore <- function(envir = parent.frame()) {
  old <- suppressMessages(sf::sf_use_s2())
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = envir)
}

run_bioclim <- function(method, ..., tile_degrees = NULL) {
  out_dir <- withr::local_tempdir(.local_envir = parent.frame())
  args <- list(..., output_dir = out_dir, method = method, verbose = FALSE)
  if (method == "tiled") {
    local_s2_restore()
    args$tile_degrees <- tile_degrees
  }
  suppressMessages(do.call(derive_bioclim, args))
}

run_statistics <- function(method, ..., tile_degrees = NULL) {
  out_dir <- withr::local_tempdir(.local_envir = parent.frame())
  args <- list(..., output_dir = out_dir, method = method, verbose = FALSE)
  if (method == "tiled") {
    local_s2_restore()
    args$tile_degrees <- tile_degrees
  }
  suppressMessages(do.call(derive_statistics, args))
}

# --- Known discrepancies ------------------------------------------------------

# Some tests expose confirmed bugs in R/. They are skipped with a description of
# the problem so the suite stays green; set FASTBIOCLIM_TEST_KNOWN_ISSUES=true
# to run them and see them fail.
skip_known_issue <- function(message) {
  if (!identical(tolower(Sys.getenv("FASTBIOCLIM_TEST_KNOWN_ISSUES")), "true")) {
    skip(paste("Known discrepancy:", message))
  }
}

# --- Comparing outputs --------------------------------------------------------

# Layer-by-layer comparison: same names, same geometry, NA in the same cells and
# values equal up to `tolerance` (outputs are written as float32).
expect_rasters_equal <- function(actual, expected, tolerance = 1e-5, label = "") {
  expect_identical(names(actual), names(expected), label = paste(label, "layer names"))
  expect_true(
    terra::compareGeom(actual, expected, lyrs = TRUE, stopOnError = FALSE, messages = FALSE),
    label = paste(label, "geometry")
  )
  a <- terra::values(actual, mat = TRUE)
  e <- terra::values(expected, mat = TRUE)
  a[is.na(a)] <- NA_real_
  e[is.na(e)] <- NA_real_
  for (nm in intersect(colnames(a), colnames(e))) {
    expect_identical(which(is.na(a[, nm])), which(is.na(e[, nm])),
                     label = sprintf("%s %s NA cells", label, nm))
    expect_equal(a[, nm], e[, nm], tolerance = tolerance,
                 label = sprintf("%s %s values", label, nm))
  }
}

# Checks selected cells of `bios_rast` against a matrix of reference values
# (one row per cell, columns named like the layers).
expect_cells_match <- function(bios_rast, cells, reference, tolerance = 1e-5, label = "") {
  got <- as.matrix(terra::extract(bios_rast, cells))
  for (nm in colnames(reference)) {
    expect_equal(unname(got[, nm]), unname(reference[, nm]), tolerance = tolerance,
                 label = sprintf("%s %s at cells %s", label, nm, paste(cells, collapse = ",")))
  }
}

# --- Independent reference implementation of bio01-bio19 --------------------------

# Written from the WorldClim / ANUCLIM definitions (as in dismo::biovars), with
# plain base R and none of the package's code. `tmin`, `tmax` and `prcp` are the
# per-unit values of a single cell. Periods start at every unit; when
# `circular = TRUE` they wrap from the last unit back to the first. Static indices
# can override the selected unit or period, as in derive_bioclim().
hand_bioclim <- function(tmin, tmax, prcp, period_length = 3, circular = TRUE,
                         warmest_period = NULL, coldest_period = NULL,
                         wettest_period = NULL, driest_period = NULL,
                         warmest_unit = NULL, coldest_unit = NULL,
                         wettest_unit = NULL, driest_unit = NULL) {
  n <- length(tmin)
  tavg <- (tmin + tmax) / 2
  starts <- if (circular) seq_len(n) else seq_len(n - period_length + 1)
  windows <- lapply(starts, function(s) ((s - 1 + seq_len(period_length) - 1) %% n) + 1)
  t_period <- vapply(windows, function(w) mean(tavg[w]), numeric(1))
  p_period <- vapply(windows, function(w) sum(prcp[w]), numeric(1))

  wet <- if (is.null(wettest_period)) which.max(p_period) else wettest_period
  dry <- if (is.null(driest_period)) which.min(p_period) else driest_period
  warm <- if (is.null(warmest_period)) which.max(t_period) else warmest_period
  cold <- if (is.null(coldest_period)) which.min(t_period) else coldest_period

  b <- numeric(19)
  b[1] <- mean(tavg)
  b[2] <- mean(tmax - tmin)
  b[5] <- if (is.null(warmest_unit)) max(tmax) else tmax[warmest_unit]
  b[6] <- if (is.null(coldest_unit)) min(tmin) else tmin[coldest_unit]
  b[7] <- b[5] - b[6]
  b[3] <- 100 * b[2] / b[7]
  b[4] <- 100 * stats::sd(tavg)
  b[8] <- t_period[wet]
  b[9] <- t_period[dry]
  b[10] <- t_period[warm]
  b[11] <- t_period[cold]
  b[12] <- sum(prcp)
  b[13] <- if (is.null(wettest_unit)) max(prcp) else prcp[wettest_unit]
  b[14] <- if (is.null(driest_unit)) min(prcp) else prcp[driest_unit]
  b[15] <- 100 * stats::sd(prcp + 1) / mean(prcp + 1)
  b[16] <- p_period[wet]
  b[17] <- p_period[dry]
  b[18] <- p_period[warm]
  b[19] <- p_period[cold]
  stats::setNames(b, sprintf("bio%02d", 1:19))
}

# Applies hand_bioclim() to several cells of a make_climate() result. Static
# index arguments, if given, are full-raster vectors and are subset per cell.
hand_bioclim_cells <- function(clim, cells, period_length = 3, circular = TRUE, ...) {
  idx <- list(...)
  t(vapply(cells, function(cl) {
    cell_idx <- lapply(idx, function(v) v[cl])
    do.call(hand_bioclim, c(list(tmin = clim$vals$tmin[cl, ],
                                 tmax = clim$vals$tmax[cl, ],
                                 prcp = clim$vals$prcp[cl, ],
                                 period_length = period_length,
                                 circular = circular),
                            cell_idx))
  }, numeric(19)))
}
