# Edge cases: NA cells, a single-cell raster, extents that do not divide evenly
# into tiles, custom periods and static indices. The terra checks are fast and
# run everywhere; anything that runs the tiled backend is skipped on CRAN.

all_inputs <- function(clim) {
  r <- clim$rasters
  list(tmin = r$tmin, tmax = r$tmax, prcp = r$prcp, srad = r$srad, mois = r$mois)
}

# --- NA cells -----------------------------------------------------------------

na_cells <- c(1L, 5L, 20L, 48L)  # a corner, a tile edge, an interior cell, the last cell

make_na_climate <- function(dir) {
  clim <- make_climate(withr::local_tempdir(.local_envir = parent.frame()), nrows = 6, ncols = 8)
  for (v in names(clim$vals)) {
    clim$vals[[v]][na_cells, ] <- NA
    clim$rasters[[v]] <- write_units(clim$vals[[v]], clim$template, dir, v)
  }
  clim
}

test_that("terra backend: cells that are NA in every layer give NA in every bio", {
  clim <- make_na_climate(withr::local_tempdir())
  res <- do.call(run_bioclim, c(list("terra", bios = 1:35), all_inputs(clim)))

  vals <- terra::values(res, mat = TRUE)
  for (nm in colnames(vals)) {
    expect_identical(which(is.na(vals[, nm])), na_cells, label = paste(nm, "NA cells"))
  }
  ok_cells <- setdiff(seq_len(48), na_cells)[1:3]
  expect_cells_match(res, ok_cells, hand_bioclim_cells(clim, ok_cells), label = "terra, NA raster")
})

test_that("tiled backend: cells that are NA in every layer match the terra backend", {
  skip_on_cran()
  clim <- make_na_climate(withr::local_tempdir())
  terra_res <- do.call(run_bioclim, c(list("terra", bios = 1:35), all_inputs(clim)))
  tiled_res <- do.call(run_bioclim, c(list("tiled", bios = 1:35, tile_degrees = 1), all_inputs(clim)))

  expect_rasters_equal(tiled_res, terra_res, label = "all-NA cells")
})

test_that("backends agree when only some layers of a cell are NA", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  dir <- withr::local_tempdir()
  # Cell 10: one month of precipitation missing.
  clim$vals$prcp[10, 5] <- NA
  # Cell 30: the first month of tmin missing. tmin_01 is the first layer the tiled
  # backend reads, so this is also the layer it uses to decide which cells to keep.
  clim$vals$tmin[30, 1] <- NA
  clim$rasters$prcp <- write_units(clim$vals$prcp, clim$template, dir, "prcp")
  clim$rasters$tmin <- write_units(clim$vals$tmin, clim$template, dir, "tmin")

  terra_res <- do.call(run_bioclim, c(list("terra", bios = 1:35), all_inputs(clim)))
  tiled_res <- do.call(run_bioclim, c(list("tiled", bios = 1:35, tile_degrees = 1), all_inputs(clim)))

  expect_rasters_equal(tiled_res, terra_res, label = "partially-NA cells")
})

# --- Single-cell raster -----------------------------------------------------------

test_that("terra backend handles a single-cell raster", {
  clim <- make_climate(withr::local_tempdir(), nrows = 1, ncols = 1)
  res <- do.call(run_bioclim, c(list("terra", bios = 1:35), all_inputs(clim)))

  expect_equal(dim(res), c(1, 1, 35))
  expect_cells_match(res, 1, hand_bioclim_cells(clim, 1), label = "terra, single cell")
})

test_that("tiled backend handles a single-cell raster", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 1, ncols = 1)
  terra_res <- do.call(run_bioclim, c(list("terra", bios = 1:35), all_inputs(clim)))
  tiled_res <- do.call(run_bioclim, c(list("tiled", bios = 1:35, tile_degrees = 5), all_inputs(clim)))

  expect_rasters_equal(tiled_res, terra_res, label = "single cell bioclim")
  expect_cells_match(tiled_res, 1, hand_bioclim_cells(clim, 1), label = "tiled, single cell")

  stats <- c("mean", "max", "min", "sum", "cv_cli", "max_period", "min_period")
  terra_stats <- run_statistics("terra", variable = clim$rasters$srad, stats = stats)
  tiled_stats <- run_statistics("tiled", variable = clim$rasters$srad, stats = stats,
                                tile_degrees = 5)
  expect_rasters_equal(tiled_stats[[sort(names(tiled_stats))]],
                       terra_stats[[sort(names(terra_stats))]], label = "single cell statistics")
})

# --- Extents that do not divide evenly into tiles ---------------------------------

test_that("backends agree when the extent does not divide evenly into tiles", {
  skip_on_cran()
  # 9 x 15 cells at 0.5 degrees = 7.5 x 4.5 degrees, with an origin off the tile grid.
  clim <- make_climate(withr::local_tempdir(), nrows = 9, ncols = 15,
                       xmin = -3.25, xmax = 4.25, ymin = 10.5, ymax = 15, seed = 3)
  terra_res <- do.call(run_bioclim, c(list("terra", bios = 1:35), all_inputs(clim)))

  # 2-degree tiles leave partial tiles on the right and top edges; 1.3-degree
  # tiles also cut through cells.
  for (td in c(2, 1.3)) {
    tiled_res <- do.call(run_bioclim, c(list("tiled", bios = 1:35, tile_degrees = td),
                                        all_inputs(clim)))
    expect_rasters_equal(tiled_res, terra_res, label = sprintf("bioclim, tile_degrees = %s", td))
  }

  stats <- c("mean", "max", "min", "sum", "cv_cli", "max_period", "min_period")
  terra_stats <- run_statistics("terra", variable = clim$rasters$mois,
                                inter_variable = clim$rasters$tmax, stats = stats,
                                inter_stats = c("max_inter", "min_inter"))
  tiled_stats <- run_statistics("tiled", variable = clim$rasters$mois,
                                inter_variable = clim$rasters$tmax, stats = stats,
                                inter_stats = c("max_inter", "min_inter"), tile_degrees = 1.3)
  expect_rasters_equal(tiled_stats[[sort(names(tiled_stats))]],
                       terra_stats[[sort(names(terra_stats))]],
                       label = "statistics, tile_degrees = 1.3")
})

# --- Custom periods -----------------------------------------------------------------

custom_period_cases <- list(
  list(n_units = 12, period_length = 4, circular = FALSE),
  list(n_units = 12, period_length = 1, circular = TRUE),
  list(n_units = 8, period_length = 2, circular = TRUE),
  list(n_units = 8, period_length = 3, circular = FALSE)
)

test_that("terra backend matches the hand reference for custom periods", {
  for (case in custom_period_cases) {
    clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8, n_units = case$n_units)
    res <- run_bioclim("terra", bios = 1:19, tmin = clim$rasters$tmin, tmax = clim$rasters$tmax,
                       prcp = clim$rasters$prcp, period_length = case$period_length,
                       circular = case$circular)
    cells <- c(1, 17, 48)
    expected <- hand_bioclim_cells(clim, cells, period_length = case$period_length,
                                   circular = case$circular)
    expect_cells_match(res, cells, expected,
                       label = sprintf("terra, n_units = %d, period_length = %d, circular = %s",
                                       case$n_units, case$period_length, case$circular))
  }
})

test_that("backends agree for custom periods", {
  skip_on_cran()
  for (case in custom_period_cases) {
    clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8, n_units = case$n_units)
    label <- sprintf("n_units = %d, period_length = %d, circular = %s",
                     case$n_units, case$period_length, case$circular)
    common <- c(list(bios = 1:35, period_length = case$period_length, circular = case$circular),
                all_inputs(clim))
    terra_res <- do.call(run_bioclim, c(list("terra"), common))
    tiled_res <- do.call(run_bioclim, c(list("tiled", tile_degrees = 1), common))
    expect_rasters_equal(tiled_res, terra_res, label = paste("bioclim,", label))

    stat_args <- list(variable = clim$rasters$srad, inter_variable = clim$rasters$prcp,
                      stats = c("max_period", "min_period"),
                      inter_stats = c("max_inter", "min_inter"),
                      period_length = case$period_length, circular = case$circular)
    terra_stats <- do.call(run_statistics, c(list("terra"), stat_args))
    tiled_stats <- do.call(run_statistics, c(list("tiled", tile_degrees = 1), stat_args))
    expect_rasters_equal(tiled_stats[[sort(names(tiled_stats))]],
                         terra_stats[[sort(names(terra_stats))]],
                         label = paste("statistics,", label))
  }
})

# --- Static indices -----------------------------------------------------------------

# Random static indices: unit indices in 1..12 and period indices in 1..12
# (12 circular quarters).
make_static_indices <- function(clim, dir, names, n_max = 12, seed = 99) {
  set.seed(seed)
  n_cells <- terra::ncell(clim$template)
  vals <- lapply(stats::setNames(names, names), function(nm) sample.int(n_max, n_cells, replace = TRUE))
  list(
    vals = vals,
    rasters = lapply(stats::setNames(names, names), function(nm) {
      write_index(vals[[nm]], clim$template, dir, nm)
    })
  )
}

period_index_names <- c("warmest_period", "coldest_period", "wettest_period", "driest_period")
unit_index_names <- c("warmest_unit", "coldest_unit", "wettest_unit", "driest_unit")
static_cells <- c(1, 12, 29, 48)

test_that("terra backend applies static period indices", {
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  idx <- make_static_indices(clim, withr::local_tempdir(), period_index_names)
  bios <- c(8:11, 16:19)

  res <- do.call(run_bioclim, c(list("terra", bios = bios, tmin = clim$rasters$tmin,
                                     tmax = clim$rasters$tmax, prcp = clim$rasters$prcp),
                                idx$rasters))
  expected <- do.call(hand_bioclim_cells, c(list(clim, static_cells), idx$vals))
  expect_cells_match(res, static_cells, expected[, sprintf("bio%02d", bios)],
                     label = "terra, static period indices")
})

test_that("tiled backend applies static period indices", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  idx <- make_static_indices(clim, withr::local_tempdir(), period_index_names)
  bios <- c(8:11, 16:19)
  common <- c(list(bios = bios, tmin = clim$rasters$tmin, tmax = clim$rasters$tmax,
                   prcp = clim$rasters$prcp), idx$rasters)

  terra_res <- do.call(run_bioclim, c(list("terra"), common))
  tiled_res <- do.call(run_bioclim, c(list("tiled", tile_degrees = 1), common))
  expected <- do.call(hand_bioclim_cells, c(list(clim, static_cells), idx$vals))

  expect_cells_match(tiled_res, static_cells, expected[, sprintf("bio%02d", bios)],
                     label = "tiled, static period indices")
  expect_rasters_equal(tiled_res, terra_res, label = "static period indices")
})

test_that("terra backend applies static unit and period indices together", {
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  idx <- make_static_indices(clim, withr::local_tempdir(), c(unit_index_names, period_index_names))

  skip_known_issue(paste(
    "bioclim_terra() forwards every static index in `...` to bio05_terra(), bio06_terra(),",
    "bio13_terra(), etc., which accept only their own index, so the call fails with",
    "'unused argument'. This happens whenever any static index is supplied and bios",
    "include 3, 5, 6, 7, 13, 14, 21, 22, 29 or 30 (e.g. bios = 1:19 with warmest_period)."
  ))
  res <- do.call(run_bioclim, c(list("terra", bios = 1:19, tmin = clim$rasters$tmin,
                                     tmax = clim$rasters$tmax, prcp = clim$rasters$prcp),
                                idx$rasters))
  expected <- do.call(hand_bioclim_cells, c(list(clim, static_cells), idx$vals))
  expect_cells_match(res, static_cells, expected, label = "terra, unit + period indices")
})

test_that("tiled backend applies static unit and period indices together", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  idx <- make_static_indices(clim, withr::local_tempdir(), c(unit_index_names, period_index_names))

  res <- do.call(run_bioclim, c(list("tiled", bios = 1:19, tile_degrees = 1,
                                     tmin = clim$rasters$tmin, tmax = clim$rasters$tmax,
                                     prcp = clim$rasters$prcp),
                                idx$rasters))
  expected <- do.call(hand_bioclim_cells, c(list(clim, static_cells), idx$vals))
  expect_cells_match(res, static_cells, expected, label = "tiled, unit + period indices")
})

test_that("terra backend ignores objects named like static indices in the global environment", {
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  # A user who built a static index earlier may still have it in the workspace.
  # Only indices passed to derive_bioclim() should be used.
  stray <- write_index(rep(1, terra::ncell(clim$template)), clim$template,
                       withr::local_tempdir(), "warmest_period")
  assign("warmest_period", stray, envir = globalenv())
  withr::defer(rm("warmest_period", envir = globalenv()))

  skip_known_issue(paste(
    "bioclim_terra() detects static indices with exists(\"warmest_period\") (and similar),",
    "which also searches enclosing environments, so an object with that name in the",
    "global environment silently replaces the computed warmest period."
  ))
  res <- run_bioclim("terra", bios = c(10, 18), tmin = clim$rasters$tmin,
                     tmax = clim$rasters$tmax, prcp = clim$rasters$prcp)
  expected <- hand_bioclim_cells(clim, static_cells)[, c("bio10", "bio18")]
  expect_cells_match(res, static_cells, expected, label = "terra, stray global warmest_period")
})

# Hand reference for derive_statistics() with static indices: values of the
# selected unit, or the mean over the selected 3-unit circular period.
statistics_static_setup <- function(dir_data, dir_idx) {
  clim <- make_climate(dir_data, nrows = 6, ncols = 8)
  idx <- make_static_indices(clim, dir_idx,
                             c("max_unit", "min_unit", "max_period", "min_period",
                               "max_interactive", "min_interactive"))
  args <- c(list(variable = clim$rasters$srad, inter_variable = clim$rasters$prcp,
                 stats = c("max", "min", "max_period", "min_period"),
                 inter_stats = c("max_inter", "min_inter"), output_prefix = "srad"),
            idx$rasters)
  srad <- clim$vals$srad
  period_mean <- function(cl, p) mean(srad[cl, ((p - 1 + 0:2) %% 12) + 1])
  expected <- t(vapply(static_cells, function(cl) c(
    srad_max = srad[cl, idx$vals$max_unit[cl]],
    srad_min = srad[cl, idx$vals$min_unit[cl]],
    srad_max_period = period_mean(cl, idx$vals$max_period[cl]),
    srad_min_period = period_mean(cl, idx$vals$min_period[cl]),
    srad_inter_high = period_mean(cl, idx$vals$max_interactive[cl]),
    srad_inter_low = period_mean(cl, idx$vals$min_interactive[cl])
  ), numeric(6)))
  list(args = args, expected = expected)
}

test_that("derive_statistics(): terra backend applies static indices", {
  setup <- statistics_static_setup(withr::local_tempdir(), withr::local_tempdir())

  skip_known_issue(paste(
    "stats_terra() checks geometry with compareGeom(lyrs = TRUE) across the variable",
    "and the single-layer static index rasters, so any static index makes it stop with",
    "'[compareGeom] number of layers does not match'."
  ))
  terra_res <- do.call(run_statistics, c(list("terra"), setup$args))
  expect_cells_match(terra_res, static_cells, setup$expected,
                     label = "terra statistics, static indices")
})

test_that("derive_statistics(): tiled backend applies static indices", {
  skip_on_cran()
  setup <- statistics_static_setup(withr::local_tempdir(), withr::local_tempdir())

  tiled_res <- do.call(run_statistics, c(list("tiled", tile_degrees = 1), setup$args))
  expect_cells_match(tiled_res, static_cells, setup$expected,
                     label = "tiled statistics, static indices")
})
