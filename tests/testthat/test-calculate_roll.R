# The base code for this test were created with Gemini Pro 2.5 (Jan 2025) using 
# as an input the function and asking to generated test. It was later modified and checked manually by
# Gonzalo E. Pinilla-Buitrago.
# tests/testthat/test-calculate_roll.R

# --- Test Setup ---

# Create a temporary directory for all test outputs
test_output_dir <- file.path(tempdir(), "test_calc_roll_mockery")
dir.create(test_output_dir, showWarnings = FALSE, recursive = TRUE)

# Create a consistent test raster: 3 years of monthly data (36 layers)
set.seed(456)
in_mem_rast <- rast(nrows = 10, ncols = 10, nlyr = 36, vals = 1:(100*36))

# Create a file-backed version for tiled workflow tests
file_backed_path <- file.path(test_output_dir, "source_raster_roll.tif")
writeRaster(in_mem_rast, file_backed_path, overwrite = TRUE)
file_backed_rast <- rast(file_backed_path)

# Create a clipping region
poly_wkt <- "POLYGON((-50 -50, 50 -50, 50 50, -50 50, -50 -50))"
clipping_region <- vect(poly_wkt, crs = crs(in_mem_rast))

# --- Tests Start Here ---

test_that("Input validation works correctly (no mocks needed)", {
  # Test for missing window_size
  expect_error(
    calculate_roll(x = in_mem_rast, freq = 12),
    "'window_size' is required."
  )
  
  # Test for non-divisible layers
  bad_rast <- in_mem_rast[[1:35]] # 35 layers is not divisible by 12
  expect_error(
    calculate_roll(x = bad_rast, window_size = 2, freq = 12),
    "Total layers .* is not divisible by frequency"
  )
  
  # Test for window_size too large
  expect_error(
    calculate_roll(x = in_mem_rast, window_size = 4, freq = 12), # 4-year window on 3 years of data
    "window_size .* cannot be larger than total cycles"
  )
  
  # Test for tiled workflow with in-memory raster
  expect_error(
    calculate_roll(x = in_mem_rast, window_size = 2, freq = 12, method = "tiled"),
    "The 'tiled' workflow requires the input SpatRaster to point to a file on disk."
  )
})


test_that("Terra (in-memory) workflow is correctly dispatched", {
  m_roll_terra <- mock("path/to/terra_result.tif")
  m_rast <- mock(rast())
  
  stub(calculate_roll, "roll_terra", m_roll_terra)
  stub(calculate_roll, "terra::rast", m_rast)

  calculate_roll(
    x = in_mem_rast,
    window_size = 2,
    freq = 12,
    method = "terra",
    verbose = FALSE
  )
  
  expect_called(m_roll_terra, 1)
  
  # Check that key arguments are passed correctly to the helper
  args <- mock_args(m_roll_terra)[[1]]
  expect_equal(args$window_size, 2)
  expect_equal(args$freq, 12)
  expect_s4_class(args$x, "SpatRaster")
})


test_that("`auto` method chooses in-memory workflow when data fits", {
  # Mocks are self-contained for this specific test case
  m_roll_terra <- mock("path_terra.tif")
  m_rast <- mock(rast())
  m_mem_info_fits <- mock(c(fits_mem = 1))
  
  stub(calculate_roll, "roll_terra", m_roll_terra)
  stub(calculate_roll, "terra::rast", m_rast)
  stub(calculate_roll, "terra::mem_info", m_mem_info_fits)
  # We don't even need to mock roll_fast here, as it shouldn't be called

  calculate_roll(
    x = file_backed_rast, 
    window_size = 2, 
    freq = 12, 
    method = "auto", 
    verbose = FALSE
  )
  
  expect_called(m_mem_info_fits, 1)
  expect_called(m_roll_terra, 1)
})

test_that("`auto` method chooses tiled workflow when data does not fit", {
  # A completely new set of mocks for this test case
  m_roll_fast <- mock("path_temp")
  m_write_layers <- mock("path_tiled.tif")
  m_rast <- mock(rast())
  m_mem_info_too_big <- mock(c(fits_mem = 0))

  stub(calculate_roll, "roll_fast", m_roll_fast)
  stub(calculate_roll, "write_layers", m_write_layers)
  stub(calculate_roll, "terra::rast", m_rast)
  stub(calculate_roll, "terra::mem_info", m_mem_info_too_big)
  # We don't need to mock roll_terra here

  calculate_roll(
    x = file_backed_rast, 
    window_size = 2, 
    freq = 12, 
    method = "auto", 
    verbose = FALSE
  )
  
  expect_called(m_mem_info_too_big, 1)
  expect_called(m_roll_fast, 1)
  expect_called(m_write_layers, 1) # Don't forget to check this!
})

test_that("`overwrite = FALSE` stops execution if files exist", {
  out_dir <- file.path(test_output_dir, "overwrite_test_roll")
  dir.create(out_dir, showWarnings = FALSE)
  
  # Manually create one of the files the function is expected to generate
  # For 3 years of data (cycles 1,2,3) with window_size=2, the first window is 1-2.
  # The first unit is 1. Freq is 12.
  # The first filename will be output_w1-2_u01.tif (default template)
  fake_output_file <- file.path(out_dir, "output_w1-2_u01.tif")
  file.create(fake_output_file)

  expect_error(
    calculate_roll(
      x = in_mem_rast,
      window_size = 2,
      freq = 12,
      output_dir = out_dir,
      overwrite = FALSE,
      verbose = FALSE
    ),
    "Output files already exist"
  )
})


test_that("`name_template` correctly generates filenames", {
  m_roll_terra <- mock("path.tif")
  m_rast <- mock(rast())
  stub(calculate_roll, "roll_terra", m_roll_terra)
  stub(calculate_roll, "terra::rast", m_rast)
  
  calculate_roll(
    x = in_mem_rast,
    window_size = 2,
    freq = 12,
    output_prefix = "test_run",
    name_template = "run_{prefix}_window_{start_window}_month_{idx_unit}",
    method = "terra",
    verbose = FALSE
  )
  
  args <- mock_args(m_roll_terra)[[1]]
  # The `output_names_list` is a list of character vectors
  expect_type(args$output_names_list, "list")
  
  # Check the first name of the first window
  first_window_names <- args$output_names_list[[1]]
  # Total cycles = 3. `start_units` is `c(1,2)`. `start_window` for first element is "1". `end_window` is "2".
  # `idx_unit` for first element is "01".
  expected_first_name <- "run_test_run_window_1_month_01"
  expect_equal(first_window_names[1], expected_first_name)
})


test_that("Clipping with `user_region` is passed to helpers", {
  m_rast <- mock(rast(), cycle = TRUE) # cycle=TRUE for two calls
  stub(calculate_roll, "terra::rast", m_rast)
  
  # --- Terra method ---
  m_roll_terra <- mock("path.tif")
  stub(calculate_roll, "roll_terra", m_roll_terra)
  
  calculate_roll(
    x = in_mem_rast, window_size = 2, freq = 12, method = "terra",
    user_region = clipping_region, verbose = FALSE
  )
  
  args_terra <- mock_args(m_roll_terra)[[1]]
  # The raster passed to the helper should be the cropped version
  expect_lt(prod(dim(args_terra$x)[1:2]), prod(dim(in_mem_rast)[1:2]))

  # --- Tiled method ---
  m_roll_fast <- mock("temp.qs")
  m_write_layers <- mock("path.tif")
  stub(calculate_roll, "roll_fast", m_roll_fast)
  stub(calculate_roll, "write_layers", m_write_layers)

  calculate_roll(
    x = file_backed_rast, window_size = 2, freq = 12, method = "tiled",
    user_region = clipping_region, verbose = FALSE
  )
  
  # The user_region object should be passed directly to roll_fast
  args_fast <- mock_args(m_roll_fast)[[1]]
  expect_false(is.null(args_fast$user_region))
  expect_s4_class(args_fast$user_region, "SpatVector")
  
  # Check that our final rast mock was called twice
  expect_called(m_rast, 2)
})

# --- Cleanup ---
unlink(test_output_dir, recursive = TRUE)
