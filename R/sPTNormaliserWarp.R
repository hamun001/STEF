

sPTNormaliserWarp <- function(rasterBrick, windowwidth,spatiaNormPercentile,output_name,cpus) {
  library("raster")
  library("rgdal")
  require("doParallel")
  windowwidth <- windowwidth
  spatiaNormPercentile <- spatiaNormPercentile
  output_name <- output_name
  cpus <- cpus
  sfQuickInit(cpus=cpus)
  rad <- rasterEngine(inraster=rasterBrick, fun=sPTNormaliser,window_dims=c(windowwidth,windowwidth),
                      args=list(spatiaNormPercentile =spatiaNormPercentile))
  writeRaster(rad, filename=output_name,datatype ="FLT4S", overwrite=TRUE)
  sfQuickStop()
}

