
bfastMonitorWarp <- function(rasterBrick, my_dates, startMonitoryear,
                             minumum_observations,magThreshold,type ="OLS-MOSUM",plot = F,cpus =1) {
  library("raster")
  library("rgdal")
  require("doParallel")
  rada <- brick(rasterBrick)
  myear <- startMonitoryear
  my_dates <- my_dates
  minumum_observations <-minumum_observations
  magThreshold <- magThreshold
  type <- type
  plot <- plot
  cpus <- cpus
  sfQuickInit(cpus=cpus)

  rad <- rasterEngine(inraster=rada,fun=ybfastmonitor ,window_dims=c(3,3),
                      args=list(myear = myear,plot = plot,
                                history = c("all"),my_dates =my_dates,
                                minumum_observations  = minumum_observations,
                                magThreshold = magThreshold,type =type))
  sfQuickStop()
  rad
}

