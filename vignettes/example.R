tmin_path <- list.files(
  "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmin",
  #"~/Downloads/tmin/",
  pattern = "\\.tif$",
  full.names = TRUE
)

tmax_path <- list.files(
  "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmax",
  #"~/Downloads/tmax",
  pattern = "\\.tif$",
  full.names = TRUE
)

prec_path <- list.files(
  "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  #"~/Downloads/prcp/",
  pattern = "\\.tif$",
  full.names = TRUE
)
future::plan("multisession", workers = 4) 
# [v10] 497 + 710 = 1207
# [v11] 453 + 716 = 1169
# [v12] 426 + 678 = 1104
# [fastbioclim] XXX + 622 = XXX
# future::plan("sequential") 
col <- AOI::aoi_get(country = "Colombia")
world <- geodata::world(resolution = 1, level = 0, path = tempdir(), 
  version = "latest")
tictoc::tic("Global calculation")
bioclim_directory_path <- fastbioclim::bioclim_vars(bios = 1:19,
                                                    n_units = 12,
                                                    tmin_path = tmin_path,
                                                    tmax_path = tmax_path,
                                                    prec_path = prec_path,
                                                    user_region = world, 
                                                    write_raw_vars = FALSE,
                                                    temp_dir = "/Users/Gonzalo/bioclim_qs")
tictoc::toc()
gc()
future::plan("sequential") 
tictoc::tic("Write raster")
write_layers(biovardir = bioclim_directory_path,
  save_dir = "/Users/Gonzalo/bioclim_r",
  clean_temporary_files = FALSE)
tictoc::toc()
r <- rast("/Users/Gonzalo/bioclim_r/bio1.tif")
plot(r)
