
spatioTemporalMetricsExtractorWarp <- function(rasterBrick, windowwidth,my_dates,spatiaNormPercentile,threshold,startMonitoryear,output_name,cpus,tryCatchError=T, density = F) {
  library("raster")
  library("rgdal")
  require("doParallel")
  #rada <- rasterBrick
  my_dates <-my_dates
  mYear <- startMonitoryear
  windowwidth <- windowwidth
  density <- density
  threshold <- threshold
  spatiaNormPercentile <- spatiaNormPercentile
  output_name <- output_name
  cpus <- cpus
  tryCatchError <- tryCatchError
  sfQuickInit(cpus=cpus)
  rad <- rasterEngine(inraster=rasterBrick, fun=spatioTemporalMetricsExtractortryCatch,window_dims=c(windowwidth,windowwidth),
                      args=list(mYear = mYear,density = density,my_dates =my_dates,threshold = threshold,spatiaNormPercentile =spatiaNormPercentile, windowwidth =windowwidth,tryCatchError=tryCatchError))
  writeRaster(rad, filename=output_name,datatype ="FLT4S", overwrite=TRUE)
  sfQuickStop()
}
