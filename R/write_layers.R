#' Export Bioclimatic Variables to GeoTIFF
#'
#' Assembles intermediate `.qs` files into a full vector in memory, performs
#' cell ID mapping if necessary, and writes final GeoTIFF rasters.
#'
#' @param biovardir Character string. Path to directory with intermediate `.qs` files and `template_info.qs`
#' @param save_dir Character string. Path for the final GeoTIFF files. Default: "bioclimatic".
#' @param clean_temporary_files Logical. Delete `biovardir` after writing? Default: FALSE.
#'
#' @return None. Writes GeoTIFF files to `save_dir`.
#' @author Luis Osorio-Olvera, Gonzalo E. Pinilla-Buitrago
#'
#' @export
#' @import terra purrr stringr qs rio

write_layers <- function(biovardir, save_dir = "bioclimatic",
                         clean_temporary_files = FALSE){

  # --- 1. Input Validation and Load Template Info ---
  if (!dir.exists(biovardir)) stop("Dir not found: ", biovardir)
  template_info_file <- file.path(biovardir, "template_info.qs")
  if (!file.exists(template_info_file)) stop("template_info.qs not found.")
  
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
  
  # --- 2. Find and Organize Intermediate BIO Files ---
  bio_paths <- list.files(biovardir, pattern = "bio\\d{1,2}_.*\\.qs$", full.names = TRUE)
  if (length(bio_paths) == 0) stop("No intermediate bioclim '.qs' files found.")
  
  bios_names <- stringr::str_extract(pattern = "bio\\d{1,2}", basename(bio_paths))
  bio_paths_df <- data.frame(bio_paths, bios_names)
  bio_paths_df <- bio_paths_df[!is.na(bio_paths_df$bios_names), ]
  if (nrow(bio_paths_df) == 0) stop("No valid intermediate files found.")
  
  bio_paths_dfL <- split(bio_paths_df, bio_paths_df$bios_names)
  writing_order <- as.numeric(stringr::str_extract(pattern = "\\d{1,2}", names(bio_paths_dfL))) |> order()
  bio_paths_dfL <- bio_paths_dfL[writing_order]
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
                   !all(terra::dim(original_template)[1:2] == terra::dim(target_template)[1:2])
  
  n_target_cells <- terra::ncell(target_template)
  
  # --- 4. Process Each Bioclim Variable ---
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  seq_along(bio_paths_dfL) |> purrr::walk(function(x) {
    bio_paths_df_base <- bio_paths_dfL[[x]]
    biovar_paths <- bio_paths_df_base$bio_paths
    current_bio_name <- bio_paths_df_base$bios_names[1]
    
    message(paste("\nAssembling", current_bio_name, "in memory..."))
    pb <- utils::txtProgressBar(min = 0, max = length(biovar_paths), style = 3, width = 50)
    rvals <- rep(NA_real_, n_target_cells)
    
    for (i in seq_along(biovar_paths)) {
      bioval <- tryCatch(rio::import(biovar_paths[i]), error = function(e) NULL)
      if (is.null(bioval) || !("cell" %in% colnames(bioval)) || nrow(bioval) == 0) {
        utils::setTxtProgressBar(pb, i)
        next
      }
      cellID <- bioval$cell
      rvals[cellID] <- bioval[[1]]
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Write to GeoTIFF
    outRast <- terra::rast(target_template, names = current_bio_name)
    output_file <- file.path(save_dir, paste0(current_bio_name, ".tif"))
    
    message("Writing GeoTIFF: ", output_file)
    write_success <- FALSE
    
    tryCatch({
      
      n_target_cols <- terra::ncol(outRast)
      n_target_rows <- terra::nrow(outRast)
      step_size_write <- min(ifelse(nchar(n_target_rows) == 1, 1, 16), n_target_rows)
      comienzo <- seq(1, n_target_rows, by = step_size_write)
      
      terra::writeStart(outRast, filename = output_file, overwrite = TRUE,
                        gdal = c("COMPRESS=DEFLATE", "NBITS=16", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"))
      
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
    
    if (write_success) {
      message("\n", current_bio_name, " written successfully.")
    } else {
      if (file.exists(output_file)) file.remove(output_file)
      warning(current_bio_name, " failed write.")
    }
    
    gc()
  })
  
  # --- 5. Final Cleanup ---
  if (clean_temporary_files) {
    unlink(biovardir, recursive = TRUE, force = TRUE)
    message("Temp cleaned: ", biovardir)
  } else {
    message("Intermediate files kept: ", biovardir)
  }
  
  message("\nAll processed GeoTIFFs written to: ", normalizePath(save_dir))
  return(invisible(NULL))
}