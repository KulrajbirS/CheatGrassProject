library(Rsagacmd)
library(terra)


saga_path <- "C:/SAGA-GIS/saga-9.1.0_x64/saga_cmd.exe"
saga <- saga_gis(saga_bin = saga_path, raster_backend = "terra")

f <- system.file("ex/elev.tif", package = "terra")
dem <- rast(f)

slope <- saga$ta_morphometry$slope_aspect_curvature(elevation = dem, unit_slope = "degree", slope = "slope.sgrd", method =  .all_outputs = FALSE)

twi <- saga$ta_preprocessor$sink_drainage_route_detection(elevation = dem, sinkroute = "sinkroute.sgrd", threshold = 1, thrsheight = 100) |> 
  saga$ta_preprocessor$sink_removal(dem = dem, sinkroute = _) |> 
  saga$ta_hydrology$flow_accumulation_recursive(elevation = _, sinkroute = "sinkroute.sgrd", flow = "flow.sgrd", .all_outputs = FALSE) |> 
  saga$ta_hydrology$topographic_wetness_index_twi(slope = slope, area = _, conv = 1, method = 1)

writeRaster(twi, "twi.tif")
