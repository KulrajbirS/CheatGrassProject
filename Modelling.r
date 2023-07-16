# Required libraries
library(lidR)
library(future)
library(tidyverse)
library(sf)
library(terra)
library(ClimateNAr) 
# library(rgdal) # old package
# library(raster) # Old package, use terra instead
# library(ClimateNA) # Not sure where this package is from...
# library(sp) # Old package, use sf instead

# Load the LiDAR data
### Note: I moved the LAS file to a folder called "01_raw"
# lidar_data <- readLAS("LiDAR Points.las")

raw_dir <- file.path("./01_raw")
raw_las <- list.files(raw_dir, pattern = ".las$", full.names = TRUE)

### I moved this from below, it can stay up at the top as a placeholder to
### get the information on the resolution of the NDVI

# Load the NDVI
### Use the "terra" package instead (i.e.: the "rast()" function)
# ndvi <- raster("Multispectral_ndvi.tif")

ndvi <- rast(list.files(raw_dir, pattern = "ndvi.tif$", full.names = TRUE))

# Now read in the las file
lidar_data <- readLAS(raw_las)

# Warning messages:
# 1: Invalid data: 64813732 points with a return number equal to 0 found. 
# 2: Invalid data: 64813732 points with a number of returns equal to 0 found.

# These can be overcome by changing all of the column values for "ReturnNumber"
# and "NumberOfReturns" to 1L (the "L" at the end represents an integer value)
lidar_data$ReturnNumber <- 1L
lidar_data$NumberOfReturns <- 1L

# Save the LAS file, then reload it as a LAS Catalog for tiling
writeLAS(lidar_data, file.path(raw_dir, "lidar_data_fixed.las"))
lidar_data <- readLAScatalog(file.path(raw_dir, "lidar_data_fixed.las"), chunk_size = 105, chunk_buffer = 0)

# Tile the raw file (makes it easier on any given PC). First, create the
# directory where the tiles will be stored
tile_dir <- file.path("./02_tile")
dir.create(tile_dir, showWarnings = FALSE)

# Adjust where the tiles will be created (i.e.: the alignment) and define the 
# output file names using the coordinates of the bottom left corner of each tile
opt_chunk_alignment(lidar_data) <- c(
  floor(mean(c(lidar_data$Min.X, lidar_data$Max.X))), 
  floor(mean(c(lidar_data$Min.Y, lidar_data$Max.Y))))
opt_output_files(lidar_data) <- file.path(tile_dir, "{XLEFT}_{YBOTTOM}")

# Set up parallel environment (the plan() function is from the future package,
# whereas the set_lidr_threads() function is from the lidR package. Both have
# their use cases, but when in doubt future parallelization will almost always
# allow for faster computations, though you will be throttled by the amount of 
# RAM you have on your PC, so adjust these at your discretion)
plan(multisession, workers = availableCores())
set_lidr_threads(1)
lidar_tiles <- catalog_retile(lidar_data)

# Transform points to WGS84

# I was unable to get this function to work within a reasonable amount of time,
# so I believe that I will need to just reprocess the imagery in Pix4D and have
# a separate point cloud from that, or more simply produce a DEM using the 
# original BC Albers LiDAR and then transform that to 4326. This will be 
# slightly less accurate, but for demonstration purposes this will suffice.

# Create a DEM in BC Albers
lidar_tiles <- readLAScatalog(tile_dir, chunk_buffer = 10.5)

# Need to classify which points are noise, remove the noise, and classify ground
# points as well. Can accomplish this in a function to be called in catalog_map():

ctg_clean <- function(las) {
  las <- classify_noise(las, ivf(res = 4, n = 15))
  las <- filter_poi(las, Classification != LASNOISE)
  las <- classify_ground(las, csf())
  return(las)
}

# Create output directory where cleaned tiles will be saved to
clean_dir <- file.path("./03_cleaned")
dir.create(clean_dir, showWarnings = FALSE)
opt_output_files(lidar_tiles) <- file.path(clean_dir, "{*}")

# Set up parallel environment and run function
plan(multisession, workers = availableCores() / 2)
set_lidr_threads(2)
lidar_clean <- catalog_map(lidar_tiles, ctg_clean)

# Create DEM using lidR package functions
lidar_clean <- readLAScatalog(clean_dir, chunk_buffer = 10.5)
dem <- rasterize_terrain(lidar_clean, res = res(ndvi)[1], algorithm = tin())

# Reproject DEM to WGS84 projection and save as a .asc file
dem_wgs84 <- project(dem, "EPSG:4326", threads = TRUE)
writeRaster(dem_wgs84, "dem_wgs84.asc", overwrite = TRUE, NAflag = -9999)
plan(sequential)

# Fix line ending type (terra package saves .asc files with Unix line endings,
# but ClimateBC/NA requires Windows line endings). This isn't documented anywhere,
# it took me a long time to learn these tricks!
txt <- readLines("dem_wgs84.asc")
txt <- gsub("-9999.0", "-9999", txt)
f <- file("dem_wgs84.asc", open = "wb")
cat(txt, file = f, sep = "\r\n")
close(f)

# Call ClimateBC from command line, need to normalize file paths in order to 
# work properly
climate_dir <- file.path("./04_climate_wgs84")
dir.create(climate_dir, showWarnings = FALSE)

ClimateNA_cmdLine(
  exe = "ClimateBC_v7.40.exe",
  wkDir = normalizePath("C:/Users/matth/Downloads/Climatebc_v740/"),
  period = "Normal_1991_2020.nrm", MSY = "Y",
  inputFile = normalizePath("C:/Users/matth/Documents/GitHub/CheatGrassProject/dem_wgs84.asc"),
  outputFile = normalizePath("C:/Users/matth/Documents/GitHub/CheatGrassProject/04_climate_wgs84/")
)

# Use climateNAr function "rasterStack" to access rasters and convert to 
# necessary numerical values
layer_ids <- list.files(file.path(climate_dir, "Normal_1991_2020Y"))
climate_layers_wgs84 <- rast(ClimateNAr::rasterStack(
  file.path(climate_dir, "Normal_1991_2020Y"), varList = layer_ids, 
  rType = "grid", vConvert = TRUE))

# Project back to BC Albers
climate_out <- file.path("./05_climate_tif")
dir.create(climate_out, showWarnings = FALSE)
climate_layers <- project(climate_layers_wgs84, "EPSG:3005", 
                          filename = paste0(climate_out, "/Normal_1991_2020Y_", names(climate_layers), ".tif"))

# # Write these as .tif files
# writeRaster(climate_layers, paste0(climate_out, "/Normal_1991_2020Y_", names(climate_layers), ".tif"),
#             names = paste0(climate_out, "/Normal_1991_2020Y_", names(climate_layers)))

# Create a DEM for the WGS84 projection and write it to a .asc file
### This is not how to create a DEM, see above.

# dem <- raster::raster(lidar_data_wgs84)
# writeRaster(dem, filename="dem.asc", format="ascii", NAflag=-9999)

# Use ClimateBC to generate climate variables
# clim_variables <- ClimateNA::climateNA("dem.asc")
# clim_variables <- as.data.frame(lapply(clim_variables, function(x) as.numeric(as.character(x))))

# Transform the rasters to BC Albers projection
# clim_variables_albers <- projectRaster(clim_variables, crs = "+init=epsg:3005")


# Create a DEM from the BC Albers LiDAR, using the same cell size and extents as the NDVI and MS data
### Already created above
# dem_albers <- raster::raster(lidar_data)
# res(dem_albers) <- res(ndvi)
# extent(dem_albers) <- extent(ndvi)

# Resample all layers to fit within the bounds of the NDVI
dem_albers <- resample(dem, ndvi, method = "bilinear", threads = TRUE)
climate_layers <- resample(climate_layers, ndvi, method = "bilinear", threads = TRUE)

# Load the field point data
### I see here that you are using the sp package. Please use the sf package
### instead as the sp package will be deprecated.
field_points <- read.csv("Survey.csv", stringsAsFactors = FALSE) |> 
  janitor::clean_names() |> 
  rename(cover = x_cover) |> 
  mutate(cover = gsub("<|>", "", cover)) |> 
  separate_wider_delim(cover, "-", names_sep = "-", too_few = "align_end") |> 
  rename(cover = last_col()) |> 
  select(northing, westing, cover) |> 
  mutate(cover = as.numeric(cover),
         presence = cover > 0) |> 
  st_as_sf(coords = c("westing", "northing"), crs = 4326) |> 
  st_transform(3005)

# Convert the data frame to a spatial points data frame
# coordinates(field_points) <- ~Westing+Northing

# Set the CRS if known
# crs(field_points) <- CRS("+init=epsg:4326") # WGS84

# If needed, transform the points to the BC Albers projection (EPSG: 3005)
# field_points <- spTransform(field_points, CRS("+init=epsg:3005"))

# Extract raster data for these points
# field_points@data <- extract(clim_variables_albers, field_points)

raster_extract <- extract(c(dem_albers, ndvi, climate_layers), vect(field_points))

# For instance, using a linear model:
lm_result <- lm(variable ~ ., data = field_points@data)
summary(lm_result)

# Use the linear model to create a map prediction
raster_model <- predict(c(dem_albers, ndvi, climate_layers), lm_result, 
                        filename = "cheatgrass_model.tif", overwrite = TRUE)
