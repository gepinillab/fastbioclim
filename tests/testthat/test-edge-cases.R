# Edge cases: NA cells, a single-cell raster, extents that do not divide evenly
# into tiles, custom periods and static indices. The terra checks are fast and
# run everywhere; anything that runs the tiled backend is skipped on CRAN.
# The rasters here are 4 x 3 degrees, so `tile_degrees = 2` gives 2 x 2 = 4 tiles.
# Tests about tiling geometry use `geometry_bios`, one bio from each code path,
# because the tiled backend's run time grows with the number of output layers.

# Temperature, combined temperature/precipitation, precipitation CV,
# precipitation periods, radiation periods and moisture periods.
geometry_bios <- c(1, 3, 8, 15, 18, 24, 27, 33, 35)

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
  tiled_res <- do.call(run_bioclim, c(list("tiled", bios = 1:35, tile_degrees = 2), all_inputs(clim)))

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
  tiled_res <- do.call(run_bioclim, c(list("tiled", bios = 1:35, tile_degrees = 2), all_inputs(clim)))

  expect_rasters_equal(tiled_res, terra_res, label = "partially-NA cells")
})

# NA policy (#21): a bio is NA in a cell when any unit of any variable it is
# computed from is NA there; bios that do not use that variable are unaffected.
test_that("a bio is NA where any unit of a variable it uses is NA (#21)", {
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  dir <- withr::local_tempdir()
  clim$vals$prcp[10, 12] <- NA  # cell 10: December precipitation missing
  clim$vals$tmin[30, 1] <- NA   # cell 30: January tmin missing
  clim$vals$srad[40, 6] <- NA   # cell 40: June radiation missing
  clim$rasters$prcp <- write_units(clim$vals$prcp, clim$template, dir, "prcp")
  clim$rasters$tmin <- write_units(clim$vals$tmin, clim$template, dir, "tmin")
  clim$rasters$srad <- write_units(clim$vals$srad, clim$template, dir, "srad")

  res <- do.call(run_bioclim, c(list("terra", bios = 1:35), all_inputs(clim)))
  na_bios <- function(cell) {
    v <- terra::values(res)[cell, ]
    sort(names(v)[is.na(v)])
  }
  bios_using <- function(var) {
    sort(sprintf("bio%02d", Filter(function(b) var %in% bio_inputs(b), 1:35)))
  }
  # tavg = (tmin + tmax) / 2, so a missing tmin also makes tavg missing.
  expect_identical(na_bios(10), bios_using("prcp"))
  expect_identical(na_bios(30), sort(union(bios_using("tmin"), bios_using("tavg"))))
  expect_identical(na_bios(40), bios_using("srad"))
  expect_true("bio12" %in% na_bios(10))
  expect_false("bio01" %in% na_bios(10))

  # A static warmest period replaces tavg as the selector of bio18, so a missing
  # tmin no longer makes bio18 NA.
  wp <- write_index(rep(5L, terra::ncell(clim$template)), clim$template, dir, "warmest_period")
  res_wp <- run_bioclim("terra", bios = c(10, 18), tmin = clim$rasters$tmin,
                        tmax = clim$rasters$tmax, prcp = clim$rasters$prcp, warmest_period = wp)
  v <- terra::values(res_wp)[30, ]
  expect_true(is.na(v[["bio10"]]))
  expect_false(is.na(v[["bio18"]]))
})

test_that("derive_statistics(): backends agree when only some layers of a cell are NA (#21)", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  dir <- withr::local_tempdir()
  clim$vals$srad[c(10, 30), c(5, 1)] <- NA  # variable: one unit missing in cells 10 and 30
  clim$vals$prcp[20, 7] <- NA               # interactive variable: one unit missing in cell 20
  srad <- write_units(clim$vals$srad, clim$template, dir, "srad")
  prcp <- write_units(clim$vals$prcp, clim$template, dir, "prcp")
  stats <- c("mean", "max", "min", "sum", "stdev", "cv_cli", "max_period", "min_period")
  args <- list(variable = srad, inter_variable = prcp, stats = stats,
               inter_stats = c("max_inter", "min_inter"), output_prefix = "srad")

  terra_res <- do.call(run_statistics, c(list("terra"), args))
  tiled_res <- do.call(run_statistics, c(list("tiled", tile_degrees = 2), args))
  expect_rasters_equal(tiled_res[[sort(names(tiled_res))]], terra_res[[sort(names(terra_res))]],
                       label = "statistics, partially-NA cells")

  v <- terra::values(terra_res)
  expect_true(all(is.na(v[c(10, 30), ])))
  # Cell 20: only the interactive statistics use prcp.
  expect_identical(sort(colnames(v)[is.na(v[20, ])]), c("srad_inter_high", "srad_inter_low"))
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
  terra_res <- do.call(run_bioclim, c(list("terra", bios = geometry_bios), all_inputs(clim)))

  # 2.6-degree tiles give a 3 x 2 grid with partial tiles on the right and top
  # edges, and their edges fall inside cells (cell edges are at -3.25 + 0.5k).
  tiled_res <- do.call(run_bioclim, c(list("tiled", bios = geometry_bios, tile_degrees = 2.6),
                                      all_inputs(clim)))
  expect_rasters_equal(tiled_res, terra_res, label = "bioclim, tile_degrees = 2.6")

  stats <- c("mean", "max", "min", "sum", "cv_cli", "max_period", "min_period")
  terra_stats <- run_statistics("terra", variable = clim$rasters$mois,
                                inter_variable = clim$rasters$tmax, stats = stats,
                                inter_stats = c("max_inter", "min_inter"))
  tiled_stats <- run_statistics("tiled", variable = clim$rasters$mois,
                                inter_variable = clim$rasters$tmax, stats = stats,
                                inter_stats = c("max_inter", "min_inter"), tile_degrees = 2.6)
  expect_rasters_equal(tiled_stats[[sort(names(tiled_stats))]],
                       terra_stats[[sort(names(terra_stats))]],
                       label = "statistics, tile_degrees = 2.6")
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
  # Two of the cases above, between them covering a non-12 unit count and a
  # non-circular window. Only period-based bios depend on these settings; these
  # cover every period selector (wettest, driest, warmest, coldest, most and
  # least moist) and every period aggregate (temperature, radiation and moisture
  # means, precipitation sums).
  period_bios <- c(8, 9, 16, 17, 26, 27, 32, 33)
  for (case in custom_period_cases[c(1, 3)]) {
    clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8, n_units = case$n_units)
    label <- sprintf("n_units = %d, period_length = %d, circular = %s",
                     case$n_units, case$period_length, case$circular)
    common <- c(list(bios = period_bios, period_length = case$period_length, circular = case$circular),
                all_inputs(clim))
    terra_res <- do.call(run_bioclim, c(list("terra"), common))
    tiled_res <- do.call(run_bioclim, c(list("tiled", tile_degrees = 2), common))
    expect_rasters_equal(tiled_res, terra_res, label = paste("bioclim,", label))

    stat_args <- list(variable = clim$rasters$srad, inter_variable = clim$rasters$prcp,
                      stats = c("max_period", "min_period"),
                      inter_stats = c("max_inter", "min_inter"),
                      period_length = case$period_length, circular = case$circular)
    terra_stats <- do.call(run_statistics, c(list("terra"), stat_args))
    tiled_stats <- do.call(run_statistics, c(list("tiled", tile_degrees = 2), stat_args))
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
  tiled_res <- do.call(run_bioclim, c(list("tiled", tile_degrees = 2), common))
  expected <- do.call(hand_bioclim_cells, c(list(clim, static_cells), idx$vals))

  expect_cells_match(tiled_res, static_cells, expected[, sprintf("bio%02d", bios)],
                     label = "tiled, static period indices")
  expect_rasters_equal(tiled_res, terra_res, label = "static period indices")
})

test_that("terra backend applies static unit and period indices together", {
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  idx <- make_static_indices(clim, withr::local_tempdir(), c(unit_index_names, period_index_names))

  res <- do.call(run_bioclim, c(list("terra", bios = 1:19, tmin = clim$rasters$tmin,
                                     tmax = clim$rasters$tmax, prcp = clim$rasters$prcp),
                                idx$rasters))
  expected <- do.call(hand_bioclim_cells, c(list(clim, static_cells), idx$vals))
  expect_cells_match(res, static_cells, expected, label = "terra, unit + period indices")
})

test_that("terra backend passes each static index only to the bio that uses it (#48)", {
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  hand_names <- c(unit_index_names, "warmest_period", "wettest_period")
  other_names <- c("high_rad_unit", "low_rad_unit", "high_mois_unit", "low_mois_unit")
  idx <- make_static_indices(clim, withr::local_tempdir(), c(hand_names, other_names))
  # Every call supplies indices that belong to other bios than the ones requested.
  bios <- c(5, 6, 10, 13, 14)
  res <- do.call(run_bioclim, c(list("terra", bios = bios, tmin = clim$rasters$tmin,
                                     tmax = clim$rasters$tmax, prcp = clim$rasters$prcp),
                                idx$rasters))
  expected <- do.call(hand_bioclim_cells, c(list(clim, static_cells), idx$vals[hand_names]))
  expect_cells_match(res, static_cells, expected[, sprintf("bio%02d", bios)],
                     label = "terra, indices of other bios supplied")

  res <- do.call(run_bioclim, c(list("terra", bios = c(21, 22, 29, 30), srad = clim$rasters$srad,
                                     mois = clim$rasters$mois), idx$rasters))
  pick <- function(v, index) vapply(static_cells, function(cl) clim$vals[[v]][cl, idx$vals[[index]][cl]], numeric(1))
  expected <- cbind(bio21 = pick("srad", "high_rad_unit"), bio22 = pick("srad", "low_rad_unit"),
                    bio29 = pick("mois", "high_mois_unit"), bio30 = pick("mois", "low_mois_unit"))
  expect_cells_match(res, static_cells, expected, label = "terra, unit indices of other bios supplied")
})

test_that("tiled backend applies static unit and period indices together", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  idx <- make_static_indices(clim, withr::local_tempdir(), c(unit_index_names, period_index_names))

  res <- do.call(run_bioclim, c(list("tiled", bios = 1:19, tile_degrees = 2,
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

  terra_res <- do.call(run_statistics, c(list("terra"), setup$args))
  expect_cells_match(terra_res, static_cells, setup$expected,
                     label = "terra statistics, static indices")
})

test_that("derive_statistics(): terra backend checks the geometry and layers of static indices (#49)", {
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  one_layer <- write_index(rep(2, terra::ncell(clim$template)), clim$template,
                           withr::local_tempdir(), "max_unit")
  # An index with more than one layer is rejected.
  expect_error(
    run_statistics("terra", variable = clim$rasters$srad, stats = "max",
                   max_unit = clim$rasters$prcp),
    "exactly one layer"
  )
  # An index on a different grid is rejected.
  other_grid <- write_index(rep(2, 12), terra::rast(nrows = 3, ncols = 4, xmin = 0, xmax = 4,
                                                    ymin = 0, ymax = 3, crs = "EPSG:4326"),
                            withr::local_tempdir(), "max_unit")
  expect_error(run_statistics("terra", variable = clim$rasters$srad, stats = "max",
                              max_unit = other_grid))
  # A valid single-layer index selects that unit.
  res <- run_statistics("terra", variable = clim$rasters$srad, stats = "max",
                        max_unit = one_layer, output_prefix = "srad")
  expect_cells_match(res, static_cells, cbind(srad_max = clim$vals$srad[static_cells, 2]),
                     label = "terra statistics, constant max_unit")
})

test_that("derive_statistics(): tiled backend applies static indices", {
  skip_on_cran()
  setup <- statistics_static_setup(withr::local_tempdir(), withr::local_tempdir())

  tiled_res <- do.call(run_statistics, c(list("tiled", tile_degrees = 2), setup$args))
  expect_cells_match(tiled_res, static_cells, setup$expected,
                     label = "tiled statistics, static indices")
})

# --- Non-rectangular user_region (#11) ----------------------------------------------

test_that("define_translate() drops source cells outside the target grid (#11)", {
  # 10-column source grid; the 4 x 4 target starts at source row 4, column 3.
  translate <- define_translate(ncol_src = 10L, ncol_tgt = 4L, row_offset = 3L,
                                col_offset = 2L, nrow_tgt = 4L)
  # Source cells at (row, col): (4, 3), (4, 6), (4, 7), (7, 3), (7, 6), (8, 3).
  expect_identical(translate(c(33, 36, 37, 63, 66, 73)), c(1L, 4L, NA, 13L, 16L, NA))
})

test_that("tiled backends match terra with a non-rectangular user_region (#11)", {
  skip_on_cran()
  # Every layer holds the source cell id, so a value written to the wrong cell
  # shows up as a different number. The L-shaped region cuts through cells on
  # all sides, so exact_extract() returns cells that crop(mask = TRUE) drops.
  template <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                          crs = "EPSG:4326")
  id <- write_units(matrix(seq_len(100), 100, 12), template, withr::local_tempdir(), "id")
  region <- sf::st_as_sf(sf::st_union(
    sf::st_as_sfc(sf::st_bbox(c(xmin = 2.3, ymin = 4, xmax = 6.3, ymax = 6.3), crs = sf::st_crs(4326))),
    sf::st_as_sfc(sf::st_bbox(c(xmin = 4, ymin = 2, xmax = 6.3, ymax = 4), crs = sf::st_crs(4326)))
  ))

  terra_stats <- run_statistics("terra", variable = id, stats = "mean", user_region = region)
  tiled_stats <- run_statistics("tiled", variable = id, stats = "mean", user_region = region,
                                tile_degrees = 2)
  expect_rasters_equal(tiled_stats, terra_stats, label = "statistics, L-shaped region")

  terra_bio <- run_bioclim("terra", bios = 12, prcp = id, user_region = region)
  tiled_bio <- run_bioclim("tiled", bios = 12, prcp = id, user_region = region, tile_degrees = 2)
  expect_rasters_equal(tiled_bio, terra_bio, label = "bioclim, L-shaped region")

  avg_args <- list(index = rep(1:6, 2), output_names = "id", user_region = region)
  terra_avg <- do.call(run_calculate, c(list(calculate_average, "terra", id), avg_args))
  tiled_avg <- do.call(run_calculate, c(list(calculate_average, "tiled", id, tile_degrees = 2),
                                        avg_args))
  expect_rasters_equal(tiled_avg, terra_avg, label = "calculate_average, L-shaped region")

  roll_args <- list(window_size = 2, freq = 6, output_prefix = "id", user_region = region)
  terra_roll <- do.call(run_calculate, c(list(calculate_roll, "terra", id), roll_args))
  tiled_roll <- do.call(run_calculate, c(list(calculate_roll, "tiled", id, tile_degrees = 2),
                                         roll_args))
  expect_rasters_equal(tiled_roll, terra_roll, label = "calculate_roll, L-shaped region")
})

# --- Layer subsets of multi-band files (#27) -------------------------------------------

# Writes `vals` (one column per unit) into bands 13-24 of a single 24-band file,
# with unrelated decoy values in bands 1-12, and returns the 12-layer subset.
write_band_subset <- function(vals, template, dir, prefix) {
  decoy <- matrix(-999, nrow(vals), 12)
  r <- terra::rast(template, nlyrs = 24)
  terra::values(r) <- cbind(decoy, vals)
  f <- file.path(dir, paste0(prefix, "_multiband.tif"))
  terra::writeRaster(r, f, overwrite = TRUE, datatype = "FLT8S")
  terra::rast(f)[[13:24]]
}

test_that("layer_sources() gives one source per layer (#27)", {
  clim <- make_climate(withr::local_tempdir(), nrows = 2, ncols = 3)
  sub <- write_band_subset(clim$vals$prcp, clim$template, withr::local_tempdir(), "prcp")
  src <- layer_sources(sub)
  expect_length(src, 12)
  expect_true(all(startsWith(src, "vrt://")))
  expect_equal(unname(terra::values(terra::rast(src))), unname(terra::values(sub)))
  # Single-band files keep their plain paths.
  expect_identical(layer_sources(clim$rasters$prcp), terra::sources(clim$rasters$prcp))
})

test_that("tiled backends read the selected layers of a multi-band file (#27)", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  dir <- withr::local_tempdir()
  sub <- Map(write_band_subset, clim$vals, list(clim$template), list(dir), names(clim$vals))

  bios <- c(1, 5, 8, 12, 16, 20, 28)
  terra_bio <- run_bioclim("terra", bios = bios, tmin = sub$tmin, tmax = sub$tmax,
                           prcp = sub$prcp, srad = sub$srad, mois = sub$mois)
  tiled_bio <- run_bioclim("tiled", bios = bios, tmin = sub$tmin, tmax = sub$tmax,
                           prcp = sub$prcp, srad = sub$srad, mois = sub$mois, tile_degrees = 2)
  expect_rasters_equal(tiled_bio, terra_bio, label = "bioclim, band subsets")
  hand_bios <- sprintf("bio%02d", c(1, 5, 8, 12, 16))
  expect_cells_match(tiled_bio, c(1, 48), hand_bioclim_cells(clim, c(1, 48))[, hand_bios],
                     label = "tiled bioclim, band subsets")

  # A static index stored as one band of a multi-band file.
  idx <- terra::rast(clim$template, nlyrs = 2)
  terra::values(idx) <- cbind(1, rep(c(3, 7), length.out = terra::ncell(idx)))
  idx_file <- file.path(dir, "indices.tif")
  terra::writeRaster(idx, idx_file, datatype = "INT2S")
  max_unit <- terra::rast(idx_file)[[2]]
  stat_args <- list(variable = sub$srad, stats = c("mean", "max"), max_unit = max_unit,
                    output_prefix = "srad")
  tiled_stats <- do.call(run_statistics, c(list("tiled", tile_degrees = 2), stat_args))
  expect_equal(unname(terra::values(tiled_stats[["srad_mean"]])[, 1]), rowMeans(clim$vals$srad),
               tolerance = 1e-5)
  expected_max <- clim$vals$srad[cbind(seq_len(48), rep(c(3, 7), length.out = 48))]
  expect_equal(unname(terra::values(tiled_stats[["srad_max"]])[, 1]), expected_max, tolerance = 1e-5)

  # A full-extent user_region, because the tiled calculate_*() backends need one (#12).
  region <- terra::vect(terra::ext(sub$prcp), crs = terra::crs(sub$prcp))
  terra_avg <- run_calculate(calculate_average, "terra", sub$prcp, index = rep(1:4, 3),
                             output_names = "p")
  tiled_avg <- run_calculate(calculate_average, "tiled", sub$prcp, index = rep(1:4, 3),
                             output_names = "p", user_region = region, tile_degrees = 2)
  expect_rasters_equal(tiled_avg, terra_avg, label = "calculate_average, band subset")
  terra_roll <- run_calculate(calculate_roll, "terra", sub$prcp, window_size = 2, freq = 4,
                              output_prefix = "p")
  tiled_roll <- run_calculate(calculate_roll, "tiled", sub$prcp, window_size = 2, freq = 4,
                              output_prefix = "p", user_region = region, tile_degrees = 2)
  expect_rasters_equal(tiled_roll, terra_roll, label = "calculate_roll, band subset")
})

# --- Out-of-range static period indices (#28) ------------------------------------------

test_that("period extractors return NA for out-of-range period indices (#28)", {
  # Two cells, three periods, then the min_idx and max_idx columns.
  periods <- cbind(period_1 = c(10, 20), period_2 = c(11, 21), period_3 = c(12, 22),
                   min_idx = c(1, 1), max_idx = c(3, 3))
  idx <- c(4, 2)  # 4 would point at min_idx
  expect_warning(res <- bio10_fast(periods, idx, period_length = 3, cell = 1:2), "out of bounds")
  expect_identical(unname(res[, "bio10"]), c(NA_real_, 21))
  expect_warning(res <- bio18_fast(periods, c(5, 0), cell = 1:2), "out of bounds")
  expect_identical(unname(res[, "bio18"]), c(NA_real_, NA_real_))
  expect_identical(unname(bio16_fast(periods, c(NA, 1), cell = 1:2)[, "bio16"]), c(NA_real_, 20))
})

test_that("backends give NA for a static period index beyond the number of periods (#28)", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  # circular = FALSE with 12 units and period_length 3 gives 10 periods; 11 is out of range.
  wp <- write_index(rep(11L, terra::ncell(clim$template)), clim$template,
                    withr::local_tempdir(), "warmest_period")
  common <- list(bios = c(10, 18, 26, 34), tmin = clim$rasters$tmin, tmax = clim$rasters$tmax,
                 prcp = clim$rasters$prcp, srad = clim$rasters$srad, mois = clim$rasters$mois,
                 circular = FALSE, warmest_period = wp)
  terra_res <- do.call(run_bioclim, c(list("terra"), common))
  tiled_res <- suppressWarnings(do.call(run_bioclim, c(list("tiled", tile_degrees = 2), common)))
  expect_true(all(is.na(terra::values(terra_res))))
  expect_true(all(is.na(terra::values(tiled_res))))
})

# --- CRS of tiled outputs (#31) ---------------------------------------------------------

test_that("tiled outputs keep the input CRS, including its EPSG code (#31)", {
  skip_on_cran()
  set.seed(31)
  # 4 x 5 cells of 10 km in EPSG:3035 (ETRS89 / LAEA Europe), whose PROJ.4 string
  # does not round-trip to the EPSG code.
  template <- terra::rast(nrows = 4, ncols = 5, xmin = 4000000, xmax = 4050000,
                          ymin = 3000000, ymax = 3040000, crs = "EPSG:3035")
  vals <- matrix(stats::runif(20 * 12, 0, 100), 20)
  x <- write_units(vals, template, withr::local_tempdir(), "prcp")
  code <- function(r) terra::crs(r, describe = TRUE)$code
  region <- terra::vect(terra::ext(x), crs = terra::crs(x))

  outputs <- list(
    bioclim = run_bioclim("tiled", bios = 12, prcp = x, tile_degrees = 1e6),
    statistics = run_statistics("tiled", variable = x, stats = "mean", tile_degrees = 1e6),
    average = run_calculate(calculate_average, "tiled", x, index = rep(1:4, 3),
                            output_names = "p", user_region = region, tile_degrees = 1e6),
    roll = run_calculate(calculate_roll, "tiled", x, window_size = 2, freq = 4,
                         output_prefix = "p", user_region = region, tile_degrees = 1e6)
  )
  for (nm in names(outputs)) {
    expect_identical(code(outputs[[nm]]), "3035", label = paste(nm, "EPSG code"))
    expect_true(terra::same.crs(outputs[[nm]], x), label = paste(nm, "same.crs"))
  }
})
