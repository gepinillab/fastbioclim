#' Export Bioclimatic Variables to GeoTIFF
#'
#' Assembles intermediate `.qs` files into a full vector in memory, performs
#' cell ID mapping if necessary, and writes final GeoTIFF rasters.
#'
#' @param input_dir Character string. Path to directory with intermediate `.qs` files and `template_info.qs`
#' @param output_dir Character string. Path for the final GeoTIFF files.
#' @param file_pattern Character string. A prefix or base pattern to identify the
#'        groups of `.qs` files. For example, if files are "bio01_1.qs", "bio01_2.qs",
#'        "bio02_1.qs", etc., `file_pattern` would be "bio". If files are
#'        "var_mean_1.qs", "var_sum_1.qs", etc., `file_pattern` could be "var".
#'        The function will attempt to extract the full variable name (e.g., "bio01", "var_mean")
#'        from the filenames. Default is "bio".
#' @param gdal_opt Character vector. GDAL creation options for the output GeoTIFF files.
#'        These options control compression, threading, and other advanced features.
#'        See the GDAL documentation for a full list of options for the GeoTIFF driver.
#'        The default is `c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS")`,
#'        which provides good, lossless compression. To disable compression, use `NULL`.
#' @param overwrite Logical. If `TRUE`, any existing GeoTIFF files in the `output_dir`
#'        with the same name will be overwritten. If `FALSE` (the default), the function
#'        will throw an error if it attempts to write to a file that already exists,
#'        which is a safety measure to prevent accidental data loss.
#' @return None. Writes GeoTIFF files to `output_dir`.
#' @author Luis Osorio-Olvera, Gonzalo E. Pinilla-Buitrago
#' @keywords internal
#' @import terra purrr stringr qs rio

write_layers <- function(input_dir, 
                         output_dir,
                         file_pattern = "bio",
                         gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
                         overwrite = FALSE) {

  # --- 1. Input Validation and Load Template Info ---
  if (!dir.exists(input_dir)) stop("Input directory not found: ", input_dir)
  template_info_file <- file.path(input_dir, "template_info.qs")
  if (!file.exists(template_info_file)) stop("template_info.qs not found in ", input_dir)
  
  template_info <- tryCatch({ 
    qs::qread(template_info_file) 
  }, error = function(e) { 
    stop("Failed read template: ",e$message) 
  })
  
  if (!is.list(template_info) || !all(c("original_geom", "target_geom") %in% names(template_info))) 
    stop("template_info.qs error.")
    
  validate_geom <- function(geom_list, name) { 
    if (!is.list(geom_list) || !all(c("extent","dimensions","crs") %in% names(geom_list))) 
      stop(name," fmt err.")
    if (length(geom_list$extent) != 4 || length(geom_list$dimensions) != 2 || !is.character(geom_list$crs)) 
      stop(name," fmt err.") 
  }
  validate_geom(template_info$original_geom, "original_geom")
  validate_geom(template_info$target_geom, "target_geom")
  
  # --- 2. Find and Organize Intermediate Files ---
  qs_file_regex <- paste0("^", file_pattern, "([^[:digit:]_]|_)+[^_]+_\\d+\\.qs$")
  all_qs_paths <- list.files(input_dir, pattern = "\\.qs$", full.names = TRUE)
  all_qs_paths <- all_qs_paths[!basename(all_qs_paths) %in% c("template_info.qs")]
  if (length(all_qs_paths) == 0) {
    stop("No intermediate '.qs' files found in ", input_dir, " (excluding template_info.qs).")
  }
  
  var_base_names <- stringr::str_extract(basename(all_qs_paths), "^.*(?=_\\d+\\.qs$)")
  if (!is.null(file_pattern) && nchar(file_pattern) > 0) {
      # Ensure var_base_names start with the file_pattern
      valid_indices <- startsWith(var_base_names, file_pattern)
      all_qs_paths <- all_qs_paths[valid_indices]
      var_base_names <- var_base_names[valid_indices]
  }
  
  if (length(all_qs_paths) == 0) {
    stop(paste0("No intermediate '.qs' files found matching the pattern starting with '", file_pattern, "'."))
  }
  qs_paths_df <- data.frame(paths = all_qs_paths, names = var_base_names, stringsAsFactors = FALSE)
  qs_paths_df <- qs_paths_df[!is.na(qs_paths_df$names), ] # Remove if extraction failed

  if (nrow(qs_paths_df) == 0) {
    stop(paste0("No valid intermediate files found to process for pattern '", file_pattern, "'."))
  }

  qs_paths_dfL <- split(qs_paths_df, qs_paths_df$names)
  if (file_pattern == "bio") {
    writing_order <- as.numeric(stringr::str_extract(pattern = "\\d{1,2}", names(qs_paths_dfL))) |> order()
  } else {
    writing_order <- order(names(qs_paths_dfL))
  }
  qs_paths_dfL <- qs_paths_dfL[writing_order]
  gc()
  
  # --- 3. Reconstruct Template Rasters ---
  create_template <- function(g) {
    terra::rast(xmin = g$extent[1], xmax = g$extent[2], 
                ymin = g$extent[3], ymax = g$extent[4], 
                nrows = g$dimensions[1], ncols = g$dimensions[2], 
                crs = g$crs)
  }
  
  original_template <- create_template(template_info$original_geom)
  target_template <- create_template(template_info$target_geom)
  
  if (terra::ncell(target_template) == 0) stop("Target template empty.")
  if (terra::ncell(original_template) == 0) stop("Original template empty.")
  
  crop_occurred <- !all(terra::ext(original_template) == terra::ext(target_template)) || 
                   !all(dim(original_template)[1:2] == dim(target_template)[1:2])
  
  n_target_cells <- terra::ncell(target_template)
  
  # --- 4. Process Each Bioclim Variable ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  seq_along(qs_paths_dfL) |> purrr::walk(function(x) {
    qs_paths_df_base <- qs_paths_dfL[[x]]
    var_paths <- qs_paths_df_base$paths
    current_qs_name <- qs_paths_df_base$names[1]
    
    rvals <- rep(NA_real_, n_target_cells)
    
    for (i in seq_along(var_paths)) {
      bioval <- tryCatch(rio::import(var_paths[i]), error = function(e) NULL)
      if (is.null(bioval) || !("cell" %in% colnames(bioval)) || nrow(bioval) == 0) {
        next
      }
      cellID <- bioval$cell
      rvals[cellID] <- bioval[[1]]
    }
    
    # Write to GeoTIFF
    outRast <- terra::rast(target_template, names = current_qs_name)
    output_file <- file.path(output_dir, paste0(current_qs_name, ".tif"))
    
    message("Writing GeoTIFF: ", output_file)
    write_success <- FALSE
    
    tryCatch({
      
      n_target_cols <- terra::ncol(outRast)
      n_target_rows <- terra::nrow(outRast)
      step_size_write <- min(ifelse(nchar(n_target_rows) == 1, 1, 16), n_target_rows)
      comienzo <- seq(1, n_target_rows, by = step_size_write)
      
      terra::writeStart(outRast, filename = output_file, overwrite = TRUE,
                        gdal = gdal_opt)
      
      pb_write <- utils::txtProgressBar(min = 0, max = length(comienzo), style = 3, width = 50)
      
      for (i in seq_along(comienzo)) {
        start_row_block <- comienzo[i]
        num_rows_in_block <- min(step_size_write, n_target_rows - start_row_block + 1)
        if (num_rows_in_block <= 0) break
        
        start_cell_block <- terra::cellFromRowCol(outRast, start_row_block, 1)
        end_cell_block <- terra::cellFromRowCol(outRast, start_row_block + num_rows_in_block - 1, n_target_cols)
        values_for_block <- rvals[start_cell_block:end_cell_block]
        
        terra::writeValues(x = outRast, v = values_for_block,
                           start = start_row_block, 
                           nrows = num_rows_in_block)
        
        utils::setTxtProgressBar(pb_write, i)
      }
      close(pb_write)
      write_success <- TRUE
    }, error = function(e) {
      warning("Write error ", output_file, ": ", e$message, call. = FALSE)
    }, finally = {
      try(terra::writeStop(outRast), silent = TRUE)
    })
    
    if (!write_success) {
      if (file.exists(output_file)) file.remove(output_file)
      warning(current_qs_name, " failed write.")
    }
    gc()
  })
  
  # --- 5. Final Cleanup ---
  # Enable the debug mode: Sys.setenv(BIOCLIM_DEBUG_RAW_VARS = "TRUE")
  clean_temporary_files <- identical(toupper(Sys.getenv("BIOCLIM_DEBUG_KEEP_TEMP_FILES")), "TRUE")
  
  if (clean_temporary_files) {
    message("DEBUG MODE: Writing raw variable tiles because BIOCLIM_DEBUG_KEEP_TEMP_FILES is set to TRUE.")
    message("Intermediate files kept: ", input_dir)
  } else {
    unlink(input_dir, recursive = TRUE, force = TRUE)
  }
  
  created_tif <- file.path(normalizePath(output_dir), paste0(names(qs_paths_dfL), ".tif"))
  return(created_tif)
}