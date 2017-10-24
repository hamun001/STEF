

#'@title Global spatial normalisation for SAR time series 
#'
#'@description Performs global spatial normalisation of satellite data to reduce seasonal variations and intersensor dfferences. For local spatial normalisation, see \code{\link[STEF]{stef_local_spatial_normaliser}} 
#'@param inraster    Input raster stack or single raster. 
#'@param isStack     Logical. Set to TRUE if the input raster is a raster stack.
#'@param xpercentile The upper  percentile to use for normalisation. Default is 0.95, which represent 95th percentile.
#'@param output_filename  The name of the output normalised raster or raster stack. The name must contain the file format (e.g. raster_normalised.tif)
#'@import spatial.tools
#'@export
#'@import rgdal
#'@import raster
#'@import doParallel
#'@importFrom lubridate decimal_date
#'@references 1. Hamunyela, E., Verbesselt, J., Herold, M. (2016) Using spatial context to improve early detection of deforestationfrom Landsat time series. Remote Sensing of Environment,172, 126â€“138.
#'   \url{http://dx.doi.org/10.1016/j.rse.2015.11.006}
#'
#' 2. Hamunyela, E., Verbesselt, J.,de Bruin, S., Herold, M. (2016). Monitoring Deforestation at Sub-Annual Scales as Extreme Events in Landsat Data Cubes. Remote Sensing, 8(8), 651.
#'  \url{http://dx.doi.org/10.3390/rs8080651}
#'
#'   3. Hamunyela, E., Reiche, J., Verbesselt, J., & Herold, M. (2017). Using space-time features to improve detection of forest disturbances from Landsat time series. Remote Sensing, 9(6), 515.
#'    \url{http://dx.doi:10.3390/rs9060515}
#
#'
#'@return 
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
#' 
#' # example:
#'
#' 
#' ## synthetic aperture radar (SAR)
#' stef_global_spatial_SAR_normaliser(ra, isStack = T, xpercentile = 0.95,output_filename ="ra_global_normalised.tif")
#' }
#'
#'@author Eliakim Hamunyela
#'@details To be completed.
#'
#'
#'


stef_global_spatial_SAR_normaliser <- function(inraster,  isStack = T, xpercentile =0.95, output_filename= NULL){
  if (xpercentile > 1){
    print("xpercentile  can't be greater than 1. defaulting to  0.95")
    xpercentile <- 0.95
  }
  
  if (!is.null(output_filename)){
    
    if (isStack){
      #inraster <- brick(inraster)
      inRast  <- subset(inraster, 1, drop = F)
      inRastInverse <- 1/abs(inRast)
      xPer <- quantile(inRastInverse,  probs = xpercentile, type=7,names = FALSE)
      norInRast <- inRastInverse/xPer
      
      for (i in 2: nlayers(inraster)){
        inRast2 <- subset(inraster, i, drop =F)
        inRastInverse2 <- 1/abs(inRast2)
        xPer2 <- quantile(inRastInverse2,  probs =xpercentile, type=7,names = FALSE)
        norInRast2 <- inRastInverse2/xPer2
        norInRast <- stack(norInRast,norInRast2)
      }
      
    }else{
      inRastInverse <- 1/abs(inraster)
      xPer <- quantile(inRastInverse,  probs = xpercentile, type=7,names = FALSE)
      norInRast <- inRastInverse/xPer
    }
    
    writeRaster(norInRast, filename=output_filename,datatype ="FLT4S", overwrite=TRUE) 
    
  }else{ 
    print("please provide the output filename with file format")
    stop
    
  }
  
} 

