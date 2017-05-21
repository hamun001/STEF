#'@title Spatially normalise satellite images to reduce seasonal variations during forest cover change detection
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
#' ##@return Returns spatially normalised image or image stack
#'@examples
#'\dontrun{
#'
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
#' rad <- rasterEngine(inraster=rasterBrick, fun=sPTNormaliser,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(spatiaNormPercentile =95))
#'
#' #paralell processing example:
#' ## register the cores
#' sfQuickInit(cpus=5)
#' rad <- rasterEngine(inraster=rasterBrick, fun=sPTNormaliser,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(spatiaNormPercentile =95))
#' writeRaster(rad, filename="test.tif",datatype ="FLT4S", overwrite=TRUE)
#' # Unregister the cores
#' sfQuickStop()
#' }

#'@author Eliakim Hamunyela
#'@details To be completed.
##@section Slots:{\code{slot1} {:} {Matrix of class \code{"numeric"}, containing data from slot1}}
## \describe{
##   \itemize{
##   \item{\code{layer 1} {-} {Matrix of class \code{"numeric"}, containing data from slot1}}
##   \item{\code{layer 2} {-} {Matrix of class \code{"numeric"}, containing data from slot1}}
##   \item{\code{layer 3} {-} {Matrix of class \code{"numeric"}, containing data from slot1}}
##  }
##  {\code{slot2}:{Object of class \code{"character"}, containing data that needs to go in slot2.}}
##   \itemize{
##  \item{\code{layer 1} {-} {Matrix of class \code{"numeric"}, containing data from slot1}}
##  }
##}
#'

LandsatDataMaskCloudsShadowWaterSnow <- function(fmasklayerList, vi = vi){
  require(raster)
  require(rgdal)
  #Check for index scences where FMASK layer is available
  fmasklayerList <- fmasklayerList
  fl_ndvi <- c()
  for (i in 1:length(fmasklayerList)){
    tx <- substr(fmasklayerList [i], 1,21)
    suffixo <- substr(fmasklayerList[i], (nchar(fmasklayerList[i])- 3), (nchar(fmasklayerList[i])))
    dtz <- paste( tx ,"_sr_",vi, suffixo, sep ="")
    fl_ndvi <- c(fl_ndvi, dtz)
  }
  # List of FMASK and NDVI Layers
  suffixo <- suffixo
  mb <- fmasklayerList
  nv <- fl_ndvi
  for (i in 1: length(mb)){
    dC <- as.character(as.vector(substr(nv[i], 1, 29))) #34
    xdr <- raster (nv[i])
    smr <- raster(mb[i])
    xa <- paste (dC,"_masked_water",suffixo,sep="")
    xdrMas1 <- mask(xdr, smr, filename =xa, maskvalue = 1,overwrite=TRUE)
    xa1 <- paste (dC,"_masked_shadow",suffixo,sep="")
    xdrMas2 <- mask(xdrMas1, smr,filename =xa1, maskvalue = 2,overwrite=TRUE)
    xa2 <- paste (dC,"_masked_snow",suffixo,sep="")
    xdrMas3 <- mask(xdrMas2, smr, filename =xa2, maskvalue = 3,overwrite=TRUE)
    xa3 <- paste (dC,"_masked_cloud",suffixo,sep="")
    xdrMas4 <- mask(xdrMas3, smr,filename =xa3, maskvalue = 4,overwrite=TRUE)
    xn <- paste (dC,"_masked_all",suffixo,sep="")
    writeRaster(xdrMas4, filename=xn, overwrite=TRUE)
    #delete unneccessary layers
    xaoo <- paste ("*masked_cloud",suffixo,sep="")
    xao1 <- paste ("*masked_snow",suffixo,sep="")
    xao2 <- paste ("*masked_shadow",suffixo,sep="")
    xao3 <- paste ("*masked_water",suffixo,sep="")
    do.call(file.remove,list(list.files(path, pattern=glob2rx(xaoo), recursive=F)))
    do.call(file.remove,list(list.files(path, pattern=glob2rx(xao1), recursive=F)))
    do.call(file.remove,list(list.files(path, pattern=glob2rx(xao2), recursive=F)))
    do.call(file.remove,list(list.files(path, pattern=glob2rx(xao3), recursive=F)))
  }
}
