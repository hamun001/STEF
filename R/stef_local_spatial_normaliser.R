

#'@title Normalise satellite images spatially  to reduce seasonal variations during forest cover change detection
#'
#'@description This function normalise satellite images spatially (see,  Hamunyela et al., 2016)
#' NOTE: This function must always be called within the \code{\link[spatial.tools]{rasterEngine}} function. See the usage example
#'@param inraster    Input raster stack or single image (e.g. NDVI). The images in the image stack can have different value range
#'@param spatiaNormPercentile  Numeric. The upper percentile ( e.g. 95th percentile ) for determining the value to use
#' for spatial normalisation. The default is 95
#'@export
#'@import rgdal
#'@import raster
#'@import doParallel
#'@importFrom lubridate decimal_date
#'@seealso \code{\link[SpaceTimeChangeDetection]{rasterEngine}}
#'@references 1. Hamunyela, E., Verbesselt, J., Herold, M. (2016) Using spatial context to improve early detection of deforestationfrom Landsat time series. Remote Sensing of Environment,172, 126–138.
#'   \url{http://dx.doi.org/10.1016/j.rse.2015.11.006}
#'
#' 2. Hamunyela, E., Verbesselt, J.,de Bruin, S., Herold, M. (2016). Monitoring Deforestation at Sub-Annual Scales as Extreme Events in Landsat Data Cubes. Remote Sensing, 8(8), 651.
#'  \url{http://dx.doi.org/10.3390/rs8080651}
#'@return Returns spatially normalised image or image stack
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
#' rad <- rasterEngine(inraster=rasterBrick, fun=stef_local_spatial_normaliser,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(spatiaNormPercentile =95))
#'
#' #paralell processing example:
#' ## register the cores
#' sfQuickInit(cpus=5)
#' rad <- rasterEngine(inraster=rasterBrick, fun=stef_local_spatial_normaliser,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(spatiaNormPercentile =95))
#' writeRaster(rad, filename="test.tif",datatype ="FLT4S", overwrite=TRUE)
#' # Unregister the cores
#' sfQuickStop()
#' }

#'@author Eliakim Hamunyela
#'@details To be completed.
#'
#'
#'
#'
stef_local_spatial_normaliser <- function(inraster,spatiaNormPercentile =95,...) {
  library("raster")
  library("rgdal")
  require("doParallel")
  a <- brick(inraster)
  dFrv <- as.data.frame(t(getValues(a)))
  proC<- (as.numeric(dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)]))
  a01 <- subset(proC, !is.na(proC))
  if (length (a01) > 0){
    spNormaliser <- function(x ){
      y <- round(as.numeric(quantile(x, c(spatiaNormPercentile/100), na.rm =T)),digits = 2)
      x <- as.numeric(x)/y
      return(x)
    }

    dFrv <- apply(dFrv, 1, spNormaliser)
    dFrv <- t(dFrv)
    proCell<- (as.numeric(dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)]))
  }else{ proCell <- proC}
  return(proCell)
}

