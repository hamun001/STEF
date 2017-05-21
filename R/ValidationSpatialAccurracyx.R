

#'@title Performances spatial accuracy assessment for change detection results, returning overall, producer's and user's accuracy.
#'
#'@description This function performances spatial accuracy assessment of change detection results by comparing the date of change for validation data and the date of change from change detection algorithm. See the usage example
#'@param x    A data frame  containing the date for change as per validation data and the date of change from Algorithm.The date for validation data and that of change results from algorithm must be in different columns, and the date must be specified as decimal date e.g. 2015.467.
#'At the location where no change occured and detected, the sample point must be allocated one common unique value.
#'@param NoDataValue    Numeric or character. A unique value allocated to sample points where no change occured or detected. It can be NA, for example
#'@param EarlyDateIsCommission      Logical. If TRUE, change detected earlier than the date from the validation data is regarded as error of commission
#'@param TotalSamplesize  Numeric. The total number of sample points in the validation data
#'@param snumberOfsamplefromChange  Numeric. Number of samples from change stratum
#'@param columWithChangeDateFromValidationData  Numeric. The column number in "x (the input data frame)" that contains the date of change as per validation data. The dates must be decimal year (2015.645).
#'@param columWithChangeDateFromAlgorithm Numeric. The column number in  "x (the input data frame)" that contains the date of change as per change detection algorithm outpu. The date stamps must be decimal year (2015.645)
#'@export
#'@import rgdal
#'@import raster
#'@import doParallel
#'@importFrom lubridate decimal_date
#'@return Returns a overall, producer's and user's accuracy
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
#' rad <- rasterEngine(inraster=rasterBrick, fun=spatioTemporalMetricsExtractortryCatch,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(mYear = 2014,density = F,my_dates =imagedate,threshold = 0.01,spatiaNormPercentile =95, windowwidth =15,tryCatchError=T))
#'
#' #paralell processing example:
#' ## register the cores
#' sfQuickInit(cpus=5)
#' rad <- rasterEngine(inraster=ra, fun=spatioTemporalMetricsExtractortryCatch,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(mYear = 2014,density = F,my_dates =imagedate,threshold = 0.01,spatiaNormPercentile =95, windowwidth =15,tryCatchError=T))
#' writeRaster(rad, filename="test.tif",datatype ="FLT4S", overwrite=TRUE)
#' # Unregister the cores
#' sfQuickStop()
#' }

#'@author Eliakim Hamunyela
#'@details To be completed.
#'
#'

spatialAccurayAssessment <- function(x,NoDataValue = NULL,  EarlyDateIsCommission = NULL, TotalSamplesize = NULL,
                                     snumberOfsamplefromChange = NULL,columWithChangeDateFromValidationData =NULL,columWithChangeDateFromAlgorithm =NULL ){
 require(lubridate)

  if (is.null(TotalSamplesize) |is.null(snumberOfsamplefromChange)){
    stop("TotalSamplesize and snumberOfsamplefromChange must not be null")

  }

  x1 <- x[,columWithChangeDateFromValidationData]
  y1 <- x[,columWithChangeDateFromAlgorithm]

  if (!is.na(NoDataValue)){
  x1 <- replace (x1,  x1 ==NoDataValue, NA)
  y1 <- replace (y1,  y1 ==NoDataValue, NA)
  }
  dtaxa <- as.data.frame(x1)
  dtaxa$y1 <- y1

  xo <- c("Overall Acc", "producer's Acc", "User's Acc")

    va0 <- c()
    for (j in 1: nrow(dtaxa)) {
      xt1 <- dtaxa$x1[j]
      xt2 <- dtaxa$y1[j]
      print(c(xt1,xt2 ))

      #check if it is true positive or true negative
      if (is.na(xt1) & is.na(xt2)){
        ba0 <- "ok"
        va0 <- c(va0, ba0)
        tl0 <- F
      }else {tl0 <- T}

      tr0 <- xt2 + 0.002
      if (tl0 == T & !is.na (xt1) & !is.na(xt2) & xt1 <= tr0 |tl0 == T & !is.na (xt1) & !is.na(xt2) & xt1 > tr0 & EarlyDateIsCommission == F){
          ba0 <- "ok"
          va0 <- c(va0, ba0)
          et3 <- F
      }else {et3 <- T}

      # check if it is an omission error
      if (tl0 == T & et3 == T &!is.na (xt1) & is.na (xt2) ){
        ba0 <- "om"
        va0 <- c(va0, ba0)
        be1 <- F
      }else { be1 <- T}

      #Check if it is a commission error

      if (tl0 == T & et3 == T & be1 ==T & !is.na (xt1) & !is.na (xt2) |tl0 == T & et3 == T & be1 ==T & !is.na (xt1) & !is.na (xt2) & xt1 > xt2 & EarlyDateIsCommission == T ){
        ba0 <- "co"
        va0 <- c(va0, ba0)
        ze0 <- F
      }else {ze0 <- T}

      if (tl0 == T & et3 == T & ze0  == T & is.na (xt1) & !is.na (xt2)){
        ba0 <- "co"
        va0 <- c(va0, ba0)
      }
    }
    # calculate overall accuracy and producer's accuracy, and user's accuracy

    pPA <- (snumberOfsamplefromChange -length(subset  (va0 , va0  == "om")))/snumberOfsamplefromChange  * 100
    uSA <- (snumberOfsamplefromChange -length(subset  (va0 , va0  == "om")))/ (snumberOfsamplefromChange -length(subset  (va0 , va0  == "om")) + length(subset  (va0 , va0  == "co")))*100
    oKA <- length(subset  (va0 , va0  == "ok"))/TotalSamplesize * 100
    dox <- c(oKA, pPA,uSA)
    xo <- rbind(xo, dox)


  return(xo)
}









