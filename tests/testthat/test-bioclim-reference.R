# Checks bio01-bio19 against values computed by hand (hand_bioclim() in
# helper-backends.R, base R only) for a few cells: the corners, an interior cell
# and cells on tile edges.

reference_cells <- c(1, 8, 21, 28, 41, 48)

test_that("terra backend matches hand-computed bio01-bio19", {
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  res <- run_bioclim("terra", bios = 1:19, tmin = clim$rasters$tmin,
                     tmax = clim$rasters$tmax, prcp = clim$rasters$prcp)

  expected <- hand_bioclim_cells(clim, reference_cells)
  expect_cells_match(res, reference_cells, expected, label = "terra")
})

test_that("tiled backend matches hand-computed bio01-bio19", {
  skip_on_cran()
  clim <- make_climate(withr::local_tempdir(), nrows = 6, ncols = 8)
  # 4 x 3 degree raster with 1-degree tiles: 12 tiles, and every reference cell
  # touches a tile edge.
  res <- run_bioclim("tiled", bios = 1:19, tmin = clim$rasters$tmin,
                     tmax = clim$rasters$tmax, prcp = clim$rasters$prcp,
                     tile_degrees = 1)

  expected <- hand_bioclim_cells(clim, reference_cells)
  expect_cells_match(res, reference_cells, expected, label = "tiled")
})

test_that("hand reference agrees with a worked example", {
  # Sanity check of the reference itself on a cell small enough to verify on paper.
  tmin <- c(1, 2, 4, 7, 10, 13, 15, 14, 11, 7, 4, 2)
  tmax <- tmin + c(8, 8, 9, 10, 11, 12, 12, 12, 11, 10, 9, 8)
  prcp <- c(90, 80, 60, 40, 20, 10, 5, 10, 30, 60, 70, 100)
  b <- hand_bioclim(tmin, tmax, prcp)

  expect_equal(b[["bio05"]], 27)                     # July tmax: 15 + 12
  expect_equal(b[["bio06"]], 1)                      # January tmin
  expect_equal(b[["bio07"]], 26)
  expect_equal(b[["bio12"]], 575)
  expect_equal(b[["bio13"]], 100)
  expect_equal(b[["bio14"]], 5)
  expect_equal(b[["bio16"]], 100 + 90 + 80)          # Dec-Jan-Feb (wraps around)
  expect_equal(b[["bio17"]], 10 + 5 + 10)            # Jun-Jul-Aug
  expect_equal(b[["bio19"]], 100 + 90 + 80)          # coldest quarter is Dec-Jan-Feb
  expect_equal(b[["bio18"]], 10 + 5 + 10)            # warmest quarter is Jun-Jul-Aug
  expect_equal(b[["bio10"]], mean(((tmin + tmax) / 2)[6:8]))
})
