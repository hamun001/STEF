

#'@title Fast Alerting: Detects forest disturbances as extreme events in space-time data cubes of satellite images. An extreme event is flaged when there is an anomaly at a pixel. 
#'
#'@description This function mean for fast alerting. It detects extreme events (observations) at sub-annual scales in the space-time data cubes of satellite images.
#' Extreme observations are identified at pixel level, and the extremeness is defined based on percentile (e.g 5th percentile). The observations in the data cube are normalised spatially to reduce seasonality, a procedure carried out within this function, if specified, otherwise it is assumed that data have been normlised already. Once an extreme event is detected,the function extracts several space-time feature from local data cube.
#' The function is intended for near real-time deforestation monitoring using satellite image time series.It is scalable to multiple cores (See usage example).
#' NOTE: This function must always be called within the \code{\link[spatial.tools]{rasterEngine}} function. See the usage example
#'@param inraster    Input raster stack. Note that this input is literally the a local spac-time data cube whose dimensions are defined by the user.
#'@param my_dates    A vector containing the acquisition date for each image/layer in the image stack.
#'@param mYear       The  start datetime for the minotoring period. This parameter has to be defined either as a decimal year (e.g. 2015.356)
#'if the start of the monitoring period is not at the begining of the year or as a integer year (e.g. 2015) if the monitoring period start at the begining of the year.
#'@param spatiaNormPercentile   The upper percentile ( e.g. 95 percentile ) for determining the value to use
#' to spatially normalise the values in the local data cube. Spatial normalisation is done to reduce seasonal variation in the data cube
#' and also to reduce inter-sensor differences in the time series (data cube) when using multi-sensor time series.
#'@param threshold  A percentile/ extremeness for determining if the observation being monitored (typically newly acquired observation) is an extreme. Typically the value should be below 5th percentile
#'@param densityPlot  Logical. If TRUE, for every pixel in the entire raster stack, the observations in the local data cube would be plotted as density function when detecting deforestation.
#'@param windowwidth  A numeric value specifying the spatial extent (deminsions) of the data cube. The data cube always has a square spatial extent, hence one numeric value as input, and
#'the windowwidth can only be defined as an odd number,larger than 1
#'@param tryCatchError   Logical. If TRUE, pixels were the algorithm faces an error during  processing would be ignored, and NA values would be return for each of such pixel. If FALSE, the process would stop.
#'@param sPatioNormalixse Logical. If TRUE, pixel values at each timestep will be normalised spatially using  spatiaNormPercentile value, and normalisation wil be local (within the data cube)
#'@param rf_modelx Model object based on training dataset. Set to NULL if no training has been done 
#'@param changeProbability Numeric. The probability threshold for accepting an anomaly as forest disturbances. It is user-defined based on the training data.
#'@export
#'@import spatial.tools
#'@import rgdal
#'@import raster
#'@import doParallel
#'@import randomFores
#'@importFrom lubridate decimal_date
#'@seealso \code{\link[spatial.tools]{rasterEngine}}
#'@seealso \code{\link[STEF]{rasterEngine}}
#'@references 1. Hamunyela, E., Verbesselt, J., Herold, M. (2016) Using spatial context to improve early detection of deforestationfrom Landsat time series. Remote Sensing of Environment,172, 126-138.
#'   \url{http://dx.doi.org/10.1016/j.rse.2015.11.006}
#'
#' 2. Hamunyela, E., Verbesselt, J.,de Bruin, S., Herold, M. (2016). Monitoring Deforestation at Sub-Annual Scales as Extreme Events in Landsat Data Cubes. Remote Sensing, 8(8), 651.
#'  \url{http://dx.doi.org/10.3390/rs8080651}
#'
#'   3. Hamunyela, E., Reiche, J., Verbesselt, J., & Herold, M. (2017). Using space-time features to improve detection of forest disturbances from Landsat time series. Remote Sensing, 9
#
#'
#'@return it either returns two values (date of disturbance and  probability of disturbance) per pixel if rf_modelx model object and changeProbability are provided, otherwise it returns 18 values per pixel, namely: (1) Numeric (e.g. 2017.353) - The date of the image in which the potential forest disturbnace is detected. (2) Number of valid observations in the local data cube over the reference period.
#' (3) Continuous value  - Standard deviation of the observations in the reference period of a local data cube. (4) Continuous value  - The threshold for identifying negative anomalies, computed from the local data cube over the reference period, corresponding to the specified percentile(5th percentile). (5) Magnitude of change (the difference between the actual threshold and the negative anomalies,
#' (5) Discrete-Numeric value - Number of consecutive negative anomalies (is either 2 for two consecutive negative anomalies  or 1 if only the last observation in the time series is a negative anomaly).
#' (6) Discrete-Numeric value - Number of 8-connected neigbours for the focal pixel which also experienced negative anomalies at T1. (7) Binary (yes =1, no =0) value - Indicates whether any of the 8-connected neighbours for the focal pixel also experienced negative anomalies at T1. (8) Discrete-Numeric value - Number of 8-connected neighbours of the focal pixel which are also experiencing negative anomalies at T2.
#' (9) Binary (yes =1, no =0) value - Indicates whether any of the 8-connected neighbours of the focal pixel is also experiencing negative anomalies at T2. (10) Discrete-Numeric value - Number of pixels in the local data cube which also experienced negative anomalies at T1.
#' (11) Discrete-Numeric value - Number of pixels in local data cube with Negative anomalies at T2. (12) Binary (yes =1, no =0) value - Indicates whether any of the 8-connected neighbours for the focal pixel is already non-forest in the reference period.
#' (13) Discrete-Numeric value - Number of the 8-connected neighbours for the focal pixel which are already non-forest in the reference period. 14) Discrete-Numeric value - Number of pixels within the local data cube which have been masked as non-forest in the reference period.
#' (16)  Continuous value  - Cumulative sum of residuals for spatial variability at T1. (17) Continuous value  -  Cumulative sum of resdiuals for spatial variability at T2. (18) Continuous value  - Temporal linear trend in spatial variability.
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
#' rad <- rasterEngine(inraster=rasterBrick, fun=stef_monitor,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(mYear = 2014,density = F,my_dates =imagedate,threshold = 0.05,spatiaNormPercentile =95, windowwidth =15,tryCatchError=T,changeProbability=0.5,rf_modelx=rf_modelx))
#'
#' #paralell processing example:
#' ## register the cores
#' sfQuickInit(cpus=5)
#' rad <- rasterEngine(inraster=ra, fun=stef_monitor,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(mYear = 2014,density = F,my_dates =imagedate,threshold = 0.05,spatiaNormPercentile =95, windowwidth =15,tryCatchError=T, sPatioNormalixse =T,changeProbability=0.5,rf_modelx=rf_modelx))
#' writeRaster(rad, filename="test.tif",datatype ="FLT4S", overwrite=TRUE)
#' # Unregister the cores
#' sfQuickStop()
#' }

#'@author Eliakim Hamunyela
#'@details To be completed.
#'
#'
#'

stef_monitor_fastAlertsProb <- function(inraster,my_dates, mYear, spatiaNormPercentile,threshold, densityPlot,windowwidth,tryCatchError,sPatioNormalixse,changeProbability,rf_modelx,...) {
  library("raster")
  library("rgdal")
  library("spatial.tools")
  require("doParallel")
  require("randomForest")
  
  spatioTemporalMetricsExtractorxcprob <- function(dFrvx,my_dates, mYear,spatiaNormPercentile,threshold, densityPlot,windowwidth,sPatioNormalixse,changeProbability,rf_modelx){
    dFrv <- dFrvx
    windowwidth <- windowwidth
    qdates <- my_dates
    proC<- (as.numeric(dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)]))
    a01 <- as.data.frame(proC)
    a01$date <- qdates
    
    a01 <- subset(a01, a01$date > mYear & !is.na(a01$proC) & !is.nan(a01$proC) & !is.infinite(a01$proC))
    xp <- NA
    CH <- NA
    dav <- as.numeric(c(NA, NA))
    if (length (a01$date) > 0){
      #rownames(dFrv) <- qdates
      dFrv$deci_date <- qdates
      de_date <- dFrv$deci_date
      dFrv$deci_date <- NA
      spNormaliser <- function(x ){
        y <- round(as.numeric(quantile(x, c(spatiaNormPercentile/100), na.rm =T)),digits = 2)
        x <- as.numeric(x)/y
        return(x)
      }
      removedips <- function (x) {
        #x <- na. approx (x, rule = 2)
        y <- as.numeric (x)
        leng <- length (x) - 2
        for(i in 1: leng )
        {
          ## moving window - check distance
          b <- i + 2
          c <- b - 1
          if(any(is.na(x[b]) ,is.na(x[c]) ,is.na(x[i])))
            next
          mida <- x[c] - x[i]
          midc <- x[c] - x[b]
          1
          # Find 20 percent
          threshold1 <- ( -1/100) * x[i]
          threshold2 <- ( -1/100) * x[b]
          # check threshold
          
          if( mida < 0 & midc < 0 & mida < threshold1 | midc < threshold2 ) {
            y[c] <- (x[b] + x[i]) / 2}
        }
        return (y)
      }
      dFrv <- as.data.frame(dFrv)
      dFrv <- dFrv[,1:(ncol(dFrv) -1)]
      
      if (sPatioNormalixse ==T){
        
        dFrv <- apply(dFrv, 1, spNormaliser)
        dFrv <- as.data.frame(dFrv)
        dFrv <- t(dFrv)
      }
      dFrv <- as.data.frame(dFrv)
      xdFrvx <- dFrv
      proCell<- (as.numeric(dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)]))
      
      #pull out the 8-connected neighbours of focal pixel
      nxs <- ncol(dFrv)-trunc(ncol(dFrv)/2)
      n2 <- as.numeric(dFrv[,nxs-1])
      n3 <- as.numeric(dFrv[,nxs+1])
      n4 <- as.numeric(dFrv[,nxs-windowwidth])
      n5 <- as.numeric(dFrv[,nxs+windowwidth])
      n6 <- as.numeric(dFrv[,nxs+windowwidth+1])
      n7 <- as.numeric(dFrv[,nxs+windowwidth-1])
      n8 <- as.numeric(dFrv[,nxs-windowwidth+1])
      n9 <- as.numeric(dFrv[,nxs-windowwidth-1])
      neig <- as.data.frame(cbind(n2,n3,n4,n5,n6,n7,n8,n9))
      #
      
      # check if one or more 8-connected neigbours are non-forest (were masked)
      n2n <- length(subset(n2, !is.na(n2)))
      n3n <- length(subset(n3, !is.na(n3)))
      n4n <- length(subset(n4, !is.na(n4)))
      n5n <- length(subset(n5, !is.na(n5)))
      n6n <- length(subset(n6, !is.na(n6)))
      n7n <- length(subset(n7, !is.na(n7)))
      n8n <- length(subset(n8, !is.na(n8)))
      n9n <- length(subset(n9, !is.na(n9)))
      don <- c(n2n, n3n, n4n, n5n,n6n, n7n, n8n,n9n)
      nonforest <- length(subset(don, don == 0))
      if (nonforest != 0){
        nonforestNeig <- 1
        NOofnonforestNeig <- nonforest
      }else{
        nonforestNeig <- 0
        NOofnonforestNeig <- 0
      }
      xdata <- as.data.frame(as.numeric(proCell))
      xdata$deci_date <- de_date
      colnames(xdata) <- c("x", "deci_date")
      atx <- subset(xdata, !is.na(xdata$x) & !is.nan(xdata$x))
      if (length(atx$x ) > 0){
        ata <- subset(xdata, !is.na(xdata$x))
        cdata <- subset(xdata, is.na(xdata$x))
        xdata <- rbind(ata, cdata)
        xdata <- xdata[order(xdata$deci_date),]
        dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)] <- removedips(xdata$x)
        dFrv$deci_date <- xdata$deci_date
        dFrv <-subset (dFrv, round(dFrv$deci_date ,digits = 4) <  mYear)
        dFrv$deci_date <- NA
        
        dtaa <- as.numeric((unlist(dFrv)))
        rm(dFrv)
        dtaax <- subset(dtaa, !is.na(dtaa) & !is.nan(dtaa) & dtaa > 0)
        rm(dtaa)
        qt <- round(as.numeric(quantile(dtaax, c(threshold), na.rm =T)),digits = 4)
        vq <- round(sd(dtaax, na.rm =T), digits = 4)
        
        # number of extreme values at each time step
        extremesIncube <- function(y,qt ){
          yo <- as.numeric(y)
          exm <- length(subset(yo, !is.na(yo) & yo < qt))
          return(exm)
        }
        ExtreCube <- as.numeric(apply(xdFrvx, 1, extremesIncube, qt))
        
        # number of non-forest pixels in the cube prior to monitoring
        mmean <- function(x){
          ox <- as.numeric(x)
          #if(length(subset(ox, !is.na(ox)) & !is.nan(ox) ) > 0){
          sa <- mean(subset(ox, !is.na(ox)))
          #}else{sa <- NA}
          return(sa)
        }
        PnFo <- as.numeric(apply(xdFrvx, 2, mmean))
        noOfNonForPixelsInCube <- length(subset(PnFo,   is.na(PnFo)))
        
        #calculate variability at each time step
        sdIncube <- function(ya){
          yax <- as.numeric(ya)
          #if(length(subset(yax, !is.na(yax) & !is.nan(yax) & yax > 0)) > 1 ){
          sdx <- sd( subset(yax, !is.na(yax) & !is.nan(yax) & yax > 0))
          #}else{sdx <- NA}
          return(sdx)
        }
        PnSD <- as.numeric(apply(xdFrvx, 1, sdIncube))
        
        # cumullative sum of the variability residuals
        Pndif <- PnSD - mean(PnSD, na.rm = T)
        PnSDo <- replace(Pndif, is.na(Pndif), 0)
        sdPn <- cumsum(PnSDo)
        
        #calculate spatio-temporal cv 
        cvIncube <- function(yx){
          yaxo <- as.numeric(yx)
          #if(length(subset(yax, !is.na(yax) & !is.nan(yax) & yax > 0)) > 1 ){
          cvx <- cv( subset(yaxo, !is.na(yaxo) & !is.nan(yaxo) & yaxo > 0))
          #}else{sdx <- NA}
          return(cvx)
        }
        PnCV <- as.numeric(apply(xdFrvx, 1, cvIncube))
        
        #calculate pixel-time series CUMSUM
        pixelMean <- median(proCell, na.rm =T)
        pixelRs <- proCell - pixelMean
        pixelRsx <- replace(pixelRs, is.na(pixelRs), 0)
        pixelCumsum <- cumsum(pixelRsx)
        rm(xdFrvx)
        
        # check the number of 8-connected neigbours whose observations are also extremes
        pacths <- function(x,qta){
          y1 <- subset(x, !is.na(x) & x < qta)
          b2 <- length(y1)
          return(b2)
        }
        dpatchx <- apply(neig, 1, pacths, qt)
        xdata$y <- as.numeric( dpatchx)
        xdata$extreme <-ExtreCube
        xdata$sdcum <- sdPn
        xdata$sdPnSD <- PnSD
        xdata$pixelCumsum <- pixelCumsum
        xdata$PnCV <- PnCV
        
        # linear trend in the variability
        if (length(subset(PnSD, !is.na(PnSD) & !is.nan(PnSD))) > 2){
          lmModelSlope <-   as.numeric((coef(lm(formula = PnSD ~ deci_date, data = xdata))))[2]
        }else{lmModelSlope <- NA}
        vqs <- length(dtaax)
        xdatap <- subset (xdata, !is.na(xdata$x) & !is.nan(xdata$x) & xdata$x < qt)
        xdatax <- subset (xdatap$deci_date , xdatap$deci_date >=  mYear)
        # Check for  extremes
        if (length(xdatap$x) > 0 & length(xdatax) !=0){
          xdata <- subset (xdata, !is.na(xdata$x) & !is.nan(xdata$x ))
          xdat1 <- subset (xdata , xdata$deci_date <  mYear & xdata$x > qt )
          if (length(xdat1$x) > 0){
            xdato <- subset (xdata$deci_date , round(xdata$deci_date, digits = 4) < mYear)
            xdata <- xdata[(length(xdato)-2):length(xdata$x),]
            xdataT <- xdata[1:2,]
            xdata <- xdata[3:length(xdata$x),]
            if (length(xdata$x) > 1){
              
              # xdata <- subset ( xdata, xdata$deci_date > mYear)
              bba <- xdata$x
              bbax <- replace(bba,  bba > qt, 0)
              xbbax <- replace(bbax,   bbax > 0, 1)
              seq_sum <- function (x){
                
                if (length(x) > 1){
                  zy <- x[1]
                  for (i in 2: length(x)){
                    
                    e <- i -1
                    fe <- x[e] + x[i]
                    zy <- c(zy, fe)
                  }
                  
                }else{
                  
                  zy <- x 
                }
                
                return(zy)
              }
              aax1 <-  seq_sum(xbbax)
              indexx <- seq (1, length(aax1))
              xdata$con <- aax1
              xdata$indexx <- indexx
              x1datx <- subset ( xdata, xdata$con > 0)
              if(length(x1datx$con) > 0){
                coxo <- length (x1datx$indexx)
                coxn <- 1
                chan <- NA
                while(is.na(chan) & coxn <= coxo ){
                  countx <- x1datx$indexx[coxn]
                if(countx == 1){
                  sxz1x <- length(xdataT$y)
                  prePatch <- xdataT$y[sxz1x]
                  preNExtremes <- xdataT$extreme[sxz1x]
                  presdcum <- xdataT$sdcum[sxz1x]
                  prePnCV <- xdata$PnCV[sxz1x]
                  prepixelCumsum <-xdata$pixelCumsum[sxz1x]
                }else{
                  
                  sxz1x <- countx -1
                  prePatch <- xdata$y[sxz1x]
                  preNExtremes <- xdata$extreme[sxz1x]
                  presdcum <- xdata$sdcum[sxz1x]
                  prePnCV <- xdata$PnCV[sxz1x]
                  prepixelCumsum <-xdata$pixelCumsum[sxz1x]
                }
                currentv <- round(xdata$deci_date[countx], digits = 4)
                xp <- round(xdata$x[countx],digits = 4)
                pxPatch <- xdata$y[countx]
                postNExtremes <- xdata$extreme[countx]
                postsdcum <- xdata$sdcum[countx]
                popixelCumsum <-xdata$pixelCumsum[countx]
                poPnCV <- xdata$PnCV[countx]
                sdTrend <- lmModelSlope
                if (prePatch !=0){
                  NboursStep1 <- 1
                  prePatch <-  prePatch
                }else{NboursStep1 <- 0}
                
                if (pxPatch != 0){
                  pxNboursStep2 <- 1
                  pxPatch <-  pxPatch
                }else{pxNboursStep2 <- 0}
                vt8 <- xp -qt
                CH <- 2
                davx <- as.numeric(c(currentv,vqs,vq,qt,vt8,CH,prePatch,NboursStep1,pxPatch, pxNboursStep2,preNExtremes,postNExtremes,
                                    nonforestNeig,NOofnonforestNeig,noOfNonForPixelsInCube,presdcum,postsdcum,sdTrend,prepixelCumsum,popixelCumsum,
                                    prePnCV,poPnCV))
                
                newdata <- data.frame( matrix(ncol = 22, nrow = 1))
                colnames( newdata)<-c ("currentv","vqs","vq","qt","vt8", "CH",
                                                           "prePatch","NboursStep1","pxPatch", "pxNboursStep2","preNExtremes","postNExtremes",
                                                           "nonforestNeig","NOofnonforestNeig","noOfNonForPixelsInCube","presdcum","postsdcum","sdTrend",
                                       "prepixelCumsum","popixelCumsum","prePnCV","poPnCV")
                newdata[1,] <-  davx
                newdata <- newdata[,2:22]
                                    
                # predicted_class <- predict(rf_modelx,newdata)
                probThreshold <- predict(rf_modelx, newdata, type = "prob")
                probThreshold <- as.numeric(round(probThreshold, digits = 4))
                probThreshold <- probThreshold[2]
                if (probThreshold > changeProbability){
                  dav <- c(currentv, probThreshold)
                  chan <- 0
                }else{
                  dav <- as.numeric(c(NA, NA))
                  chan <- NA
                  coxn <-coxn + 1
                }
                }
              }else{
                dav <- as.numeric(c(NA,NA))
              } 

            }else{
              currentv <- NA
              dav <- as.numeric(c(NA, NA))
            }
          }else{
            currentv <- NA
            dav <- as.numeric(c(NA, NA))
          }
        }else{
          currentv <- NA
          dav <- as.numeric(c(NA, NA))
        }
        if (densityPlot == T){
          plot(density(dtaax))
          abline(v = xp, col = "blue", lty = 2)
          abline(v = qt, col = "red", lty = 2)
        }
      }else{
        dav <- as.numeric(c(NA, NA))
      }
    }else{dav <- as.numeric(c(NA, NA))}
    return (dav)
  }
  
  
  a <- brick(inraster)
  dFrvx <- as.data.frame(t(getValues(a)))
  my_dates <- my_dates
  mYear <-mYear
  spatiaNormPercentile <- spatiaNormPercentile
  threshold <- threshold
  density <- density
  windowwidth <- windowwidth
  tryCatchError <- tryCatchError
  changeProbability <- changeProbability
  rf_modelx <- rf_modelx
  
  if(tryCatchError){
    result <-  tryCatch (spatioTemporalMetricsExtractorxcprob(dFrvx,my_dates, mYear,spatiaNormPercentile,threshold, density,windowwidth,sPatioNormalixse,changeProbability,rf_modelx),
                         error=function(e) as.numeric(c(NA, NA)),
                         warning=function(w) as.numeric(c(NA, NA)))
    if (inherits(result, c("error","warning"))){
      result <- as.numeric(c(NA, NA))
    } else {
      result <- result
    }
  }else{
    
    result <- spatioTemporalMetricsExtractorxcprob(dFrvx,my_dates, mYear,spatiaNormPercentile,threshold, density,windowwidth,sPatioNormalixse,changeProbability,rf_modelx)
  }
  return(result)
}

