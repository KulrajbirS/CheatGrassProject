# Required libraries
library(sf)
library(rgdal)
library(raster)
library(ClimateNA)
library(sp)

# Load the LiDAR data
lidar_data <- readLAS("LiDAR Points.las")

# Transform points to WGS84
lidar_data_wgs84 <- st_transform(lidar_data, 4326)

# Create a DEM for the WGS84 projection and write it to a .asc file
dem <- raster::raster(lidar_data_wgs84)
writeRaster(dem, filename="dem.asc", format="ascii", NAflag=-9999)

# Use ClimateBC to generate climate variables
clim_variables <- ClimateNA::climateNA("dem.asc")
clim_variables <- as.data.frame(lapply(clim_variables, function(x) as.numeric(as.character(x))))

# Transform the rasters to BC Albers projection
clim_variables_albers <- projectRaster(clim_variables, crs = "+init=epsg:3005")

# Load the NDVI and multispectral bands
ndvi <- raster("Multispectral_ndvi.tif")

# Create a DEM from the BC Albers LiDAR, using the same cell size and extents as the NDVI and MS data
dem_albers <- raster::raster(lidar_data)
res(dem_albers) <- res(ndvi)
extent(dem_albers) <- extent(ndvi)

# Resample all layers to fit within the bounds of the NDVI and multispectral bands
dem_albers <- resample(dem_albers, ndvi)
clim_variables_albers <- resample(clim_variables_albers, ndvi)

# Load the field point data
field_points <- read.csv("Survey.csv", stringsAsFactors=FALSE)

# Convert the data frame to a spatial points data frame
coordinates(field_points) <- ~Westing+Northing

# Set the CRS if known
crs(field_points) <- CRS("+init=epsg:4326") # WGS84

# If needed, transform the points to the BC Albers projection (EPSG: 3005)
field_points <- spTransform(field_points, CRS("+init=epsg:3005"))

# Extract raster data for these points
field_points@data <- extract(clim_variables_albers, field_points)

# Now you can proceed with creating your models

# For instance, for a linear model:
lm_result <- lm(variable ~ ., data = field_points@data)
summary(lm_result)
