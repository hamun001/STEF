


bfastMonitorSpaNormWarp <- function(rasterBrick, my_dates, startMonitoryear,windowwidth,spatiaNormPercentile,output_name,
                             minumum_observations,magThreshold,type ="OLS-MOSUM",plot = F,cpus =1) {
  library("raster")
  library("rgdal")
  require("doParallel")
  myear <- startMonitoryear
  my_dates <- my_dates
  windowwidth <- windowwidth
  minumum_observations <-minumum_observations
  magThreshold <- magThreshold
  type <- type
  plot <- plot
  cpus <- cpus
  spatiaNormPercentile <- spatiaNormPercentile
  output_name <- output_name
  sfQuickInit(cpus=cpus)

  rad <- rasterEngine(inraster=rasterBrick,fun=ybfastmonitorNorm ,window_dims=c(windowwidth,windowwidth),
                      args=list(myear = myear,plot = plot,
                                history = c("all"),my_dates =my_dates,spatiaNormPercentile =spatiaNormPercentile ,
                                minumum_observations  = minumum_observations,
                                magThreshold = magThreshold,type =type))
  writeRaster(rad, filename=output_name,datatype ="FLT4S", overwrite=TRUE)
 sfQuickStop()
}

