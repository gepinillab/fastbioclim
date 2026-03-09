# The base code for this test were created with Gemini Pro 2.5 (Jan 2025) using 
# as an input the function and asking to generated test. It was later modified and checked manually by
# Gonzalo E. Pinilla-Buitrago.
# # tests/testthat/test-calculate_average.R

# # tests/testthat/test-calculate_average.R

# Create a temporary directory for all test outputs
test_output_dir <- file.path(tempdir(), "test_calc_avg_mockery")
dir.create(test_output_dir, showWarnings = FALSE, recursive = TRUE)

# Create a consistent test raster
set.seed(123)
in_mem_rast <- rast(nrows = 10, ncols = 10, nlyr = 12, vals = 1:1200)

# Create a file-backed version for tiled workflow tests
file_backed_path <- file.path(test_output_dir, "source_raster.tif")
writeRaster(in_mem_rast, file_backed_path, overwrite = TRUE)
file_backed_rast <- rast(file_backed_path)

# Create a sample index and clipping region
season_index <- rep(1:4, each = 3)
poly_wkt <- "POLYGON((-50 -50, 50 -50, 50 50, -50 50, -50 -50))"
clipping_region <- vect(poly_wkt, crs = crs(in_mem_rast))


# --- Tests Start Here ---

test_that("Input validation works correctly (no mocks needed)", {
  expect_error(
    calculate_average(x = 1:10, index = season_index),
    "Input 'x' must be a SpatRaster."
  )
  expect_error(
    calculate_average(x = in_mem_rast, index = 1:5),
    regexp = "Length of 'index'"
  )
  # This error is checked *before* helpers are called, so it's a good unit test.
  expect_error(
    calculate_average(x = in_mem_rast, index = season_index, method = "tiled"),
    "The 'tiled' workflow requires the input SpatRaster to point to a file on disk."
  )
})


test_that("Terra (in-memory) workflow is correctly dispatched", {
  # Mock 1: The internal helper function.
  # It will return a fake file path, as before.
  m_avg_terra <- mock("path/to/terra_result.tif")
  
  # Mock 2: The final terra::rast() call that causes the error.
  # We make it return a dummy, empty SpatRaster object.
  m_rast <- mock(rast()) 
  
  # We use stub() to replace BOTH functions for the duration of this test.
  stub(calculate_average, "aggregate_terra", m_avg_terra)
  stub(calculate_average, "terra::rast", m_rast) # <-- THE FIX

  # Now, we run the function. It will complete without error.
  calculate_average(
    x = in_mem_rast,
    index = season_index,
    method = "terra",
    verbose = FALSE,
    overwrite = FALSE # Explicitly test the default
  )
  
  # Assert that our MOCK for the helper was called exactly once.
  expect_called(m_avg_terra, 1)
  
  # Assert that our MOCK for terra::rast was also called once with the fake path.
  expect_called(m_rast, 1)
  expect_equal(mock_args(m_rast)[[1]][[1]], "path/to/terra_result.tif")
  
  # We can still inspect the arguments passed to the helper mock.
  # This confirms the wrapper is passing data correctly to the helper.
  args <- mock_args(m_avg_terra)[[1]]
  expect_s4_class(args$x, "SpatRaster")
  expect_equal(args$index, season_index)
  # The default for overwrite is FALSE, so let's test that.
  expect_false(args$overwrite) 
})


test_that("Tiled (out-of-core) workflow is correctly dispatched", {
  # Mock 1 & 2: The internal helpers for the tiled workflow
  m_avg_fast <- mock("path/to/temp_qs")
  m_write_layers <- mock("path/to/final_result.tif")
  
  # Mock 3: The final terra::rast() call that causes the error
  m_rast <- mock(rast()) # Return a dummy SpatRaster
  
  # Stub all three functions
  stub(calculate_average, "aggregate_fast", m_avg_fast)
  stub(calculate_average, "write_layers", m_write_layers)
  stub(calculate_average, "terra::rast", m_rast)

  calculate_average(
    x = file_backed_rast,
    index = season_index,
    method = "tiled",
    verbose = FALSE,
    overwrite = TRUE
  )
  
  # Assert that the helpers were called in order
  expect_called(m_avg_fast, 1)
  expect_called(m_write_layers, 1)
  
  # Assert that the final rast() call was made
  expect_called(m_rast, 1)

  # Check that the output of the first helper was passed as input to the second
  args_write <- mock_args(m_write_layers)[[1]]
  expect_equal(args_write$input_dir, "path/to/temp_qs")

  # Check that the output of the second helper was passed to the final rast() call
  args_rast <- mock_args(m_rast)[[1]]
  expect_equal(args_rast[[1]], "path/to/final_result.tif")
})


test_that("`auto` method correctly chooses workflow based on mocked mem_info", {
  # Mock the helpers for both workflows
  m_avg_terra <- mock("terra_path.tif")
  m_avg_fast <- mock("temp_path")
  m_write_layers <- mock("tiled_path.tif")
  
  # Mock the final failing call
  m_rast <- mock(rast())
  stub(calculate_average, "terra::rast", m_rast)

  # --- Case 1: Mock mem_info to say data FITS in memory ---
  m_mem_info_fits <- mock(c(fits_mem = 1))
  
  stub(calculate_average, "terra::mem_info", m_mem_info_fits)
  stub(calculate_average, "aggregate_terra", m_avg_terra)
  stub(calculate_average, "aggregate_fast", m_avg_fast) 

  calculate_average(
    x = file_backed_rast,
    index = season_index,
    method = "auto",
    verbose = FALSE,
    overwrite = TRUE
  )
  
  expect_called(m_mem_info_fits, 1)
  expect_called(m_avg_terra, 1)
  expect_called(m_avg_fast, 0)
  expect_called(m_rast, 1) # Check it was called
  expect_equal(mock_args(m_rast)[[1]][[1]], "terra_path.tif") # Check with correct path


  # --- Case 2: Mock mem_info to say data DOES NOT FIT in memory ---
  # Reset call counts on mocks by re-creating them for clarity
  m_avg_terra <- mock("terra_path.tif")
  m_avg_fast <- mock("temp_path")
  m_write_layers <- mock("tiled_path.tif")
  m_rast <- mock(rast())
  
  m_mem_info_too_big <- mock(c(fits_mem = 0))
  
  stub(calculate_average, "terra::mem_info", m_mem_info_too_big)
  stub(calculate_average, "aggregate_terra", m_avg_terra)
  stub(calculate_average, "aggregate_fast", m_avg_fast)
  stub(calculate_average, "write_layers", m_write_layers)
  stub(calculate_average, "terra::rast", m_rast)

  calculate_average(
    x = file_backed_rast,
    index = season_index,
    method = "auto",
    verbose = FALSE,
    overwrite = TRUE
  )
  
  expect_called(m_mem_info_too_big, 1)
  expect_called(m_avg_terra, 0)
  expect_called(m_avg_fast, 1)
  expect_called(m_write_layers, 1)
  expect_called(m_rast, 1) # Check it was called
  expect_equal(mock_args(m_rast)[[1]][[1]], "tiled_path.tif") # Check with correct path
})


test_that("`overwrite = FALSE` stops execution if files exist (no mocks needed)", {
  out_dir <- file.path(test_output_dir, "overwrite_test")
  dir.create(out_dir, showWarnings = FALSE)
  
  # Create a fake output file
  fake_output_file <- file.path(out_dir, "avg_unit_01.tif")
  file.create(fake_output_file)

  expect_error(
    calculate_average(
      x = in_mem_rast,
      index = season_index,
      output_dir = out_dir,
      overwrite = FALSE,
      verbose = FALSE
    ),
    "Output files already exist and `overwrite` is FALSE."
  )
})

test_that("Clipping with `user_region` is passed to helpers", {
  # Mock the final failing call. 
  # Set cycle = TRUE so it can be called more than once.
  m_rast <- mock(rast(), cycle = TRUE) # <-- THE FIX
  stub(calculate_average, "terra::rast", m_rast)

  # --- Terra method ---
  # Re-stub the terra helper for this specific call
  m_avg_terra <- mock("path.tif")
  stub(calculate_average, "aggregate_terra", m_avg_terra)
  
  calculate_average(
    x = in_mem_rast,
    index = season_index,
    method = "terra",
    user_region = clipping_region,
    verbose = FALSE,
    overwrite = TRUE
  )
  
  args <- mock_args(m_avg_terra)[[1]]
  expect_lt(prod(dim(args$x)[1:2]), prod(dim(in_mem_rast)[1:2]))

  # --- Tiled method ---
  # Re-stub the tiled helpers for this specific call
  m_avg_fast <- mock("temp.qs")
  m_write_layers <- mock("path.tif")
  stub(calculate_average, "aggregate_fast", m_avg_fast)
  stub(calculate_average, "write_layers", m_write_layers)

  calculate_average(
    x = file_backed_rast,
    index = season_index,
    method = "tiled",
    user_region = clipping_region,
    verbose = FALSE
  )
  
  args_fast <- mock_args(m_avg_fast)[[1]]
  expect_false(is.null(args_fast$user_region))
  expect_s4_class(args_fast$user_region, "SpatVector")
  
  # Now we can correctly assert that our rast mock was called twice
  expect_called(m_rast, 2)
})

# --- Cleanup ---
unlink(test_output_dir, recursive = TRUE)