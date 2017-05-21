
#'@title Detects forest cover disturbances as structural change in the time series of satellite images using \code{\link[bfast]{bfastmonitor}}  (Verbesselt et al., 2012)
#'
#'@description This function detects forest cover disturbances as structure change at sub-annual scales in the time series of satellite images. The function is a wrapper around \code{\link[bfast]{bfastmonitor}} function,
#'but this implementation allows for spatial normalisation of the time series to reduce seasonality before apply \code{\link[bfast]{bfastmonitor}} in sequential manner (see Hamunyela et al., 2016).
#'Spatial normalisation is performed within this function.The function is intended for near real-time deforestation monitoring using satellite image time series.
#'NOTE: This function is meant to be called within the \code{\link[spatial.tools]{rasterEngine}} function. See the usage example below.
#'
#'@param inraster    Input raster stack. Note that this input is literally the a local space-time data cube around the pixel being processes, whose dimensions are defined by the user.
#'@param my_dates    A vector containing the acquisition date for each image/layer in the image stack.
#'@param myear      The  start datetime for the minotoring period. This parameter has to be defined either a decimal year (e.g. 2015.356)
#'if the start of the monitoring period is not at the begining of the year or as a integer year (e.g. 2015) if the monitoring period start at the begining of the year
#'@param spatiaNormPercentile  The upper percentile ( e.g. 95 percentile ) for determining the value to use
#' to spatially normalise the values in the local data cube. Spatial normalisation is done to reduce seasonal variation in the data cube
#' and also to reduce inter-sensor differences in the time series (data cube) when using multi-sensor time series.
#'@param history  Specification of the start of the stable history period.See \code{\link[bfast]{bfastmonitor}} function. Here, the default is "all", which means all observations in the history period would be used for stable period.
#'@param plotT  Logical. If TRUE, the time series for each pixel would be plotted during the monitoring.
#'@param minumum_observations The minimum number of valid observations that must be in the pixel time series for such pixel to be analysed for forest cover disturbances. Set to 0 if there is no restriction on the minimum observations.
#'@param magThreshold   Magnitude of change threshold for accepting or rejecting the detected change is forest cover distrubance. It can only be a negative value. Set to 0 if  each negative change detected must be accepted as true forest cover change.
#'@param type     Character specifying the type of monitoring process e.g. "OLS-MOSUM" or "OLS-CUSUM". See \code{\link[bfast]{bfastmonitor}} function for more details.
#'@export
#'@import rgdal
#'@import raster
#'@import doParallel
#'@import zoo
#'@importFrom lubridate decimal_date
#'@seealso \code{\link[bfast]{bfastmonitor}}
#'@references 1. Hamunyela, E., Verbesselt, J., Herold, M. (2016) Using spatial context to improve early detection of deforestationfrom Landsat time series. Remote Sensing of Environment,172, 126–138.
#'   \url{http://dx.doi.org/10.1016/j.rse.2015.11.006}
#'
#' 2. Verbesselt J, Zeileis A, Herold M (2012). Near real-time disturbance detection using satellite image time series. Remote Sensing Of Environment, 123, 98–108.
#'  \url{http://dx.doi.org/10.1016/j.rse.2012.02.022}
#'@return Returns a vector containing (1) the timing of a forest cover disturbance (the datetime of the image in which the disturbance is detected), (2)
#'  and the magnitude of change (the difference between the predicted  and observed value)
#'
#'@examples
#'\dontrun{
#'
#' #create a raster stack
#' ra <- raster(ncols=360, nrows=180)
#' ra[] <- rnorm(ncell(ra))
#' for (i in 1: 86){
#' ro <- raster(ncols=360, nrows=180)
#' ro[] <- rnorm(ncell(ro))
#' ra <- stack(ra, ro)
#' }
#'  #generate fake  date stamps for layers
#' imagedate <- decimal_date(seq(as.Date("2010-1-1"), as.Date("2016-12-30"), by = "month"))
#'
#' #single pixel processing example:
#' rad <- rasterEngine(inraster=rasterBrick, fun=ybfastmonitorNorm,window_dims=c(15,15),
#'                    args=list(myear = 2014,,my_dates =imagedate,spatiaNormPercentile =95,plotT = F,history = c("all"), minumum_observations = 5, magThreshold =0,type = "OLS-MOSUM"))
#'
#' #paralell processing example:
#' ## register the cores
#' sfQuickInit(cpus=5)
#' rad <- rasterEngine(inraster=ra, fun=ybfastmonitorNorm,window_dims=c(15,15),
#'                    args=list(myear = 2014,,my_dates =imagedate,spatiaNormPercentile =95,plotT = F,history = c("all"), minumum_observations = 5, magThreshold =0,type = "OLS-MOSUM"))
#' writeRaster(rad, filename="test.tif",datatype ="FLT4S", overwrite=TRUE)
#' # Unregister the cores
#' sfQuickStop()
#' }

#'@author Eliakim Hamunyela
#'@details To be completed.
#'
#'


#BFASTMONITOR SET_UP
ybfastmonitorNorm <- function(inraster,myear = myear,plotT = F,history = c("all"),my_dates =my_dates,spatiaNormPercentile,
                          minumum_observations = minumum_observations,magThreshold = magThreshold,type =type, spatialNormalise=spatialNormalise) {
  library("raster")
  library("rgdal")
  library("Hmisc")
  library("signal")
  library("plyr")
  library("lubridate")
  library("zoo")
  library("bfast")

  a <- brick(inraster)
  dFrv <- t(getValues(a))

  spNormaliser <- function(x ){
    y <- round(as.numeric(quantile(x, c(spatiaNormPercentile/100), na.rm =T)),digits = 2)
    x <- as.numeric(x)/y
    return(x)
  }

  if (spatialNormalise == T){
  dFrv <- apply(dFrv, 1, spNormaliser)
  dFrv <- t(dFrv)
  }
  proCell<- (as.numeric(dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)]))
  zv <- subset(as.numeric(proCell), !is.na(as.numeric(proCell)))
  magnitudex <- 0
  if (length (zv) > minumum_observations){
    dtat <- as.data.frame(proCell)
    dtat$dates <- my_dates
    dtat$my_date <- decimal_date(dtat$dates)
    dtat <- subset(dtat, !is.na(dtat$proCell))
    dtat$proCell<- removedips(dtat$proCell )
    historyPeriod2 <- subset(dtat, round (decimal_date(dtat$dates), digits = 5) < myear)
      tim1 <- c(round (decimal_date(dtat$dates), digits =0))
      tim2 <- unique(tim1)
      timex <- subset (tim2, tim2 >= myear)
      coun <- 1
      breakpointx  <- NA
      while (is.na(breakpointx) & coun < length(timex)){
        bpts <- bfastts(dtat$proCell,  dtat$dates , type = c("irregular"))
        stmon <- timex[coun]
        bfm <- bfastmonitor(data = bpts, formula = response ~ 1, start=c(stmon, 1),plot = plotT, history = history,type =type)
        breakpointx <- round(bfm$breakpoint, digit = 5)
        if (!is.na(breakpointx)){
          pred <- mean(as.numeric((fitted(bfm$model))))
          breakpointx <- round(bfm$breakpoint, digit = 5)
          tbreak <- breakpointx - 0.004
          xbreak <- breakpointx + 0.004
          cdte <- subset(dtat, round(dtat$my_date , digits = 2) >= round( tbreak , digits =2) & round(dtat$my_date , digits = 3) <= round( xbreak  , digits =3) )
          observedx <- cdte$proCell
          magnitudex <-  round(observedx - pred,digits =5)
          if (magnitudex < 0 & magnitudex < magThreshold){ #
            breakpointx <- breakpointx
            magnitudex <- magnitudex
          }else {
            breakpointx <- NA
            magnitudex <- NA
          }

        }else {
          breakpointx <- NA
          magnitudex <- NA
        }
        coun <- coun +1
      }
  }else {
    breakpointx <- NA
    magnitudex <- NA
  }

  output <- c(as.numeric(breakpointx),as.numeric(magnitudex))
  return(output)
}
