# Checks that the in-memory ("terra") and on-disk ("tiled") backends give the
# same results. The rasters are 12 x 16 cells at 0.5 degrees (8 x 6 degrees), so
# `tile_degrees = 2` splits them into 4 x 3 = 12 tiles and `tile_degrees = 4`
# into 2 x 2 = 4 tiles. The tiled backend calls gc() once per tile and once per
# output layer, which dominates its run time on rasters this small, so only the
# all-35-bios test uses 12 tiles.

test_that("derive_bioclim(): terra and tiled backends agree for all 35 bios", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 12, ncols = 16)
  inputs <- clim$rasters

  terra_res <- run_bioclim("terra", bios = 1:35, tmin = inputs$tmin, tmax = inputs$tmax,
                           prcp = inputs$prcp, srad = inputs$srad, mois = inputs$mois)
  tiled_res <- run_bioclim("tiled", bios = 1:35, tmin = inputs$tmin, tmax = inputs$tmax,
                           prcp = inputs$prcp, srad = inputs$srad, mois = inputs$mois,
                           tile_degrees = 2)

  expect_equal(terra::nlyr(terra_res), 35)
  expect_identical(names(terra_res), sprintf("bio%02d", 1:35))
  expect_rasters_equal(tiled_res, terra_res, label = "bioclim tiled vs terra")
})

test_that("derive_bioclim(): backends agree when tavg is supplied directly", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 12, ncols = 16, seed = 7)
  inputs <- clim$rasters
  tavg <- write_units((clim$vals$tmin + clim$vals$tmax) / 2 + 0.3,
                      clim$template, withr::local_tempdir(), "tavg")
  bios <- c(1, 2, 4, 8:11, 18, 19)

  terra_res <- run_bioclim("terra", bios = bios, tmin = inputs$tmin, tmax = inputs$tmax,
                           tavg = tavg, prcp = inputs$prcp)
  tiled_res <- run_bioclim("tiled", bios = bios, tmin = inputs$tmin, tmax = inputs$tmax,
                           tavg = tavg, prcp = inputs$prcp, tile_degrees = 4)

  expect_rasters_equal(tiled_res, terra_res, label = "bioclim (tavg supplied)")
})

test_that("derive_statistics(): terra and tiled backends agree", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 12, ncols = 16)
  inputs <- clim$rasters
  stats <- c("mean", "max", "min", "sum", "cv_cli", "max_period", "min_period")
  inter_stats <- c("max_inter", "min_inter")

  for (period_stats in c("mean", "sum")) {
    terra_res <- run_statistics("terra", variable = inputs$srad, inter_variable = inputs$prcp,
                                stats = stats, inter_stats = inter_stats,
                                period_stats = period_stats, output_prefix = "srad")
    tiled_res <- run_statistics("tiled", variable = inputs$srad, inter_variable = inputs$prcp,
                                stats = stats, inter_stats = inter_stats,
                                period_stats = period_stats, output_prefix = "srad",
                                tile_degrees = 4)

    expect_setequal(names(terra_res),
                    c("srad_mean", "srad_max", "srad_min", "srad_sum", "srad_cv",
                      "srad_max_period", "srad_min_period", "srad_inter_high", "srad_inter_low"))
    # The two backends return the layers in different orders, so align by name.
    expect_rasters_equal(tiled_res[[sort(names(tiled_res))]], terra_res[[sort(names(terra_res))]],
                         label = sprintf("statistics (period_stats = '%s')", period_stats))
  }
})

# 'stdev' is checked against an independent reference on each backend: the
# sample standard deviation, the convention derive_bioclim() uses for bio04 and
# bio31 on both backends. If both backends match it, they match each other.
stdev_reference <- function(clim, cells) {
  matrix(apply(clim$vals$srad[cells, ], 1, stats::sd), ncol = 1,
         dimnames = list(NULL, "srad_stdev"))
}

test_that("derive_statistics(): 'stdev' on the terra backend is the sample standard deviation", {
  clim <- make_climate(withr::local_tempdir(), nrows = 12, ncols = 16)
  terra_res <- run_statistics("terra", variable = clim$rasters$srad, stats = "stdev",
                              output_prefix = "srad")
  cells <- c(1, 50, 192)

  expect_cells_match(terra_res, cells, stdev_reference(clim, cells),
                     label = "terra stdev vs sample sd")
})

test_that("derive_statistics(): 'stdev' on the tiled backend is the sample standard deviation", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 12, ncols = 16)
  tiled_res <- run_statistics("tiled", variable = clim$rasters$srad, stats = "stdev",
                              output_prefix = "srad", tile_degrees = 4)
  cells <- c(1, 50, 192)

  expect_cells_match(tiled_res, cells, stdev_reference(clim, cells),
                     label = "tiled stdev vs sample sd")
})

test_that("derive_statistics(): 'stdev' is the same on both backends (#22)", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 12, ncols = 16)
  terra_res <- run_statistics("terra", variable = clim$rasters$srad, stats = "stdev",
                              output_prefix = "srad")
  tiled_res <- run_statistics("tiled", variable = clim$rasters$srad, stats = "stdev",
                              output_prefix = "srad", tile_degrees = 4)
  expect_rasters_equal(tiled_res, terra_res, label = "stdev tiled vs terra")
  expect_equal(unname(terra::values(terra_res)[, 1]), apply(clim$vals$srad, 1, stats::sd),
               tolerance = 1e-5)
})
