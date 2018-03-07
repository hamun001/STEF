
#'@title Detects forest regrowth as non-extremes in space-time data cubes of satellite images after forest disturbance
#'
#'@description This function detects forest regrowth as non-extremes in the space-time data cubes of satellite images.
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
#'@export
#'@import spatial.tools
#'@import rgdal
#'@import raster
#'@import doParallel
#'@import segmented
#'@import strucchange
#'@importFrom lubridate decimal_date
#'@seealso \code{\link[spatial.tools]{rasterEngine}}
#'@seealso \code{\link[STEF]{rasterEngine}}
#'@references 1. Hamunyela, E., Verbesselt, J., Herold, M. (2016) Using spatial context to improve early detection of deforestationfrom Landsat time series. Remote Sensing of Environment,172, 126â€“138.
#'   \url{http://dx.doi.org/10.1016/j.rse.2015.11.006}
#'
#' 2. Hamunyela, E., Verbesselt, J.,de Bruin, S., Herold, M. (2016). Monitoring Deforestation at Sub-Annual Scales as Extreme Events in Landsat Data Cubes. Remote Sensing, 8(8), 651.
#'  \url{http://dx.doi.org/10.3390/rs8080651}
#'
#'   3. Hamunyela, E., Reiche, J., Verbesselt, J., & Herold, M. (2017). Using space-time features to improve detection of forest disturbances from Landsat time series. Remote Sensing, 9
#
#'
#'@return Per pixel, it returns 18 values. (1) Numeric (e.g. 2017.353) - The date of the image in which the potential forest disturbnace is detected. (2) Number of valid observations in the local data cube over the reference period.
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
#'                    args=list(mYear = 2014,density = F,my_dates =imagedate,threshold = 0.01,spatiaNormPercentile =95, windowwidth =15,tryCatchError=T))
#'
#' #paralell processing example:
#' ## register the cores
#' sfQuickInit(cpus=5)
#' rad <- rasterEngine(inraster=ra, fun=stef_monitor,window_dims=c(windowwidth=15,windowwidth =15),
#'                    args=list(mYear = 2014,density = F,my_dates =imagedate,threshold = 0.01,spatiaNormPercentile =95, windowwidth =15,tryCatchError=T, sPatioNormalixse =T))
#' writeRaster(rad, filename="test.tif",datatype ="FLT4S", overwrite=TRUE)
#' # Unregister the cores
#' sfQuickStop()
#' }

#'@author Eliakim Hamunyela
#'@details To be completed.
#'
#'
#'

stef_regrowth_monitor2 <- function(inraster,my_dates, mYear, spatiaNormPercentile,threshold, densityPlot,windowwidth,tryCatchError,sPatioNormalixse,minimumRegrowthTime,min_time_length_afterDisturbance,...) {
  library("raster")
  library("rgdal")
  library("spatial.tools")
  require("doParallel")
  
  spatioTemporalMetricsExtractorxc_regrowth_monitor2 <- function(dFrvx,my_dates, mYear,spatiaNormPercentile,threshold, densityPlot,windowwidth,sPatioNormalixse, minimumRegrowthTime, min_time_length_afterDisturbance){
    dFrv <- dFrvx
    windowwidth <- windowwidth
    qdates <- my_dates
    proC<- (as.numeric(dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)]))
    a01 <- as.data.frame(proC)
    a01$date <- qdates
    
    a01 <- subset(a01, a01$date > mYear & !is.na(a01$proC) & !is.nan(a01$proC) & !is.infinite(a01$proC))
    xp <- NA
    CH <- NA
    som <- rep(NA, 19)
    samo <- c(rep(NA, 19))
    dav <- as.numeric(samo)
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
        dtaax <- subset(dtaa, !is.na(dtaa) & !is.nan(dtaa))
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
          sdx <- sd( subset(yax, !is.na(yax) & !is.nan(yax)))
          #}else{sdx <- NA}
          return(sdx)
        }
        
        PnSD <- as.numeric(apply(xdFrvx, 1, sdIncube))
        
        # cumullative sum of the variability residuals
        PnSDxo <- as.data.frame(PnSD)
        PnSDxo$date <- xdata$deci_date
        PnSDxo1 <- subset(PnSDxo, PnSDxo$date < mYear)
        
        Pndif <- PnSD - mean(PnSDxo1$PnSD, na.rm = T)
        PnSDo <- replace(Pndif, is.na(Pndif), 0)
        sdPn <- cumsum(PnSDo)
        
        #calculate spatio-temporal cv 
        cvIncube <- function(yx){
          yaxo <- as.numeric(yx)
          #if(length(subset(yax, !is.na(yax) & !is.nan(yax) & yax > 0)) > 1 ){
          cvx <- cv( subset(yaxo, !is.na(yaxo) & !is.nan(yaxo)))
          #}else{sdx <- NA}
          return(cvx)
        }
        PnCV <- as.numeric(apply(xdFrvx, 1, cvIncube))
        
        #calculate pixel-time series CUMSUM
        timRx <- subset (xdata, xdata$deci_date <  mYear)
        proTel <- timRx$x
        pixelMean <- median(proTel, na.rm =T)
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
        xdata$vquotent <- xdata$sdcum / xdata$pixelCumsum
        
        
        # linear trend in the variability
        if (length(subset(PnSD, !is.na(PnSD) & !is.nan(PnSD))) > 2){
          lmModelSlope <-   as.numeric((coef(lm(formula = PnSD ~ deci_date, data = xdata))))[2]
        }else{lmModelSlope <- NA}
        vqs <- length(dtaax)
        xdatap <- subset (xdata, !is.na(xdata$x) & !is.nan(xdata$x) & xdata$x < qt)
        xdatax <- subset (xdatap$deci_date , xdatap$deci_date >=  mYear)
        
        if (length(xdatap$x) > 0 & length(xdatax) !=0){
          xdatamx <- subset (xdata, !is.na(xdata$x) & !is.nan(xdata$x ))
          xdat1 <- subset (  xdatamx ,   xdatamx$deci_date <  mYear &   xdatamx$x > qt )
          if (length(xdat1$x) > 0){
            xdato <- subset (xdatamx$deci_date , round(xdatamx$deci_date, digits = 4) < mYear)
            xdatamx <-   xdatamx[length(xdato):length(xdatamx$x),]
            if (length(xdatamx$x) > 1){
              countx <- 1
              chan <- NA
              # Check for consecutive extremes
              while ( is.na(chan) & countx < (length(xdatamx$x) - 2)){
                sxz1 <- countx+1
                if (xdatamx$x[countx] < qt & xdatamx$x[sxz1] < qt){
                  currentv <- round(xdatamx$deci_date[sxz1], digits = 4)
                  date_of_last_observa <- xdatamx$deci_date[length(xdatamx$deci_date)]
                  
                  diff_distarbaDate_LastDate <- date_of_last_observa - currentv
                  
                  if (diff_distarbaDate_LastDate > min_time_length_afterDisturbance){
                    #extract values after the detected disturbance 7 check for brreakpoints
                    after_disturba <- subset(xdata, xdata$deci_date > currentv)
                    #bp_ts <- breakpoints(after_disturba$sdcum ~ after_disturba$deci_date)
                    
                    # fit a segmented regression
                    #out.glm<-glm(sdcum~deci_date,data=after_disturba)
                    #o_seg<-segmented(out.glm,seg.Z=~deci_date,psi= c(after_disturba$deci_date[bp_ts$breakpoints]),control=seg.control(display=FALSE))
                    
                    #break_dates <- after_disturba$deci_date[bp_ts$breakpoints]
                    #date_for_last_break <- tail(break_dates, 1)
                    
                    #observ_last_seg <- subset(subset(xdatamx, xdatamx$deci_date >= date_for_last_break))
                    
                    bba <- after_disturba$x
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
                    
                    print("I passed here")
                    print(aax1[length(aax1)])
                    print(aax1[length(aax1) -1] )
                    print("What")
                    if (aax1[length(aax1)] != 2 | aax1[length(aax1) -1] != 2){
                    
                      
                    print("ALL Just fine")  
                    after_disturba$indexz <-seq(1, length(after_disturba$x),  by= 1)
                    
                    after_disturba$choco <- aax1
                    
                    mz_con <- subset(after_disturba,  after_disturba$choco ==2)
                    print("Yeah")
                    if (length(mz_con$indexz) != 0){
                    print("Ok")
                    anded <- mz_con$indexz[length(mz_con$indexz) ]
                    
                    non_ext <-subset(after_disturba,  after_disturba$indexz > anded)
                    print("Goooo")
                    }else{
                      
                    non_ext <- after_disturba
                      
                    }
                    
                    time_length_after_regrowth_start <- non_ext$deci_date[length(non_ext$deci_date)]  - non_ext$deci_date[1] 
                    print("I was here")
                    if (time_length_after_regrowth_start >= minimumRegrowthTime){
                      print("Hello world")
                      lastExtremes <- after_disturba$extreme[length(after_disturba$extreme)]
                      extremes_at_disturbance <- xdata$extreme[sxz1]
                      extremes_at_start_regrowth <- non_ext$extreme[1]
                      
                      prsdcum_at_disturbance <- xdata$sdcum[sxz1]
                      prsdcum_last <- after_disturba$sdcum[length(after_disturba$sdcum)]
                      prsdcum_at_start_regrowth <- non_ext$sdcum[1]
                      print("shall see")
                      
                      pixelCumsum_at_disturbance <- xdata$pixelCumsum[sxz1]
                      pixelCumsum_last <- after_disturba$pixelCumsum[length(after_disturba$pixelCumsum)]
                      pixelCumsum_at_start_regrowth <- non_ext$pixelCumsum[1]
                      print("I am home")
                      prPnCV_at_disturbance <- xdata$PnCV[sxz1]
                      prPnCV_last <- after_disturba$PnCV[length(after_disturba$PnCV)]
                      prPnCV_at_start_regrowth <- non_ext$prPnCV[1]
                      
                      patch_at_disturbance <- xdata$y[sxz1]
                      patch_last <- after_disturba$y[length(after_disturba$y)]
                      patch_at_start_regrowth <- non_ext$patch[1]
                      print("Welcome")
                      sdTrend <- lmModelSlope
                                          
                      dav <- c(currentv,diff_distarbaDate_LastDate,time_length_after_regrowth_start, sdTrend,patch_at_disturbance,patch_at_start_regrowth,patch_last, 
                               extremes_at_disturbance,extremes_at_start_regrowth,lastExtremes,prsdcum_at_disturbance,prsdcum_at_start_regrowth ,prsdcum_last,
                               pixelCumsum_at_disturbance,pixelCumsum_at_start_regrowth,pixelCumsum_last,prPnCV_at_disturbance,prPnCV_at_start_regrowth,prPnCV_last)                     
                      
                      print(dav)
                      chan <- 0
                      
                    }else{
                      dav <- as.numeric(samo)
                      chan <- 0
                      
                    }
                    
                    # start from here
#                     
#                     sz_con <- length(subset(aax1, aax1 ==2))
#                     if (sz_con == 0){
#                       
#                       #pull out the slope for spatioCumsum
#                       
#                       spatio_seg <- slope(o_seg)
#                       last_segmet_slope_spatio <- spatio_seg$deci_date[nrow(spatio_seg$deci_date)]
#                       
#                       # fit a segmented regression to pixel_cumsum
#                       
#                       bp_ts_pixel <- breakpoints(after_disturba$pixelCumsum ~ after_disturba$deci_date)
#                       
#                       out.glm_pixel<-glm(pixelCumsum~deci_date,data=after_disturba)
#                       o_seg_pixel<-segmented(out.glm_pixel,seg.Z=~deci_date,psi= c(after_disturba$deci_date[bp_ts_pixel$breakpoints]),control=seg.control(display=FALSE))
#                       
#                       sar_seg <- slope(o_seg_pixel)
#                       last_segmet_slope_pixel <- sar_seg$deci_date[nrow(sar_seg$deci_date)]
#                       
#                       # fit a segmented regression to spatioCumsum/pixel_cumsum
#                       bp_ts_vquotent <- breakpoints(after_disturba$vquotent ~ after_disturba$deci_date)
#                       
#                       out.glm_vquotent<-glm(vquotent~deci_date,data=after_disturba)
#                       o_seg_vquotent<-segmented(out.glm_vquotent,seg.Z=~deci_date,psi= c(after_disturba$deci_date[bp_ts_vquotent$breakpoints]),control=seg.control(display=FALSE))
#                       
#                       vquotent_seg <- slope(o_seg_vquotent)
#                       last_segmet_slope_vquotent <- vquotent_seg$deci_date[nrow(vquotent_seg$deci_date)]
#                       secondLast_segmet_slope_vquotent <- vquotent_seg$deci_date[(nrow(vquotent_seg$deci_date) -1)]
#                       
#                       # pull out other metrics
#                       
#                       lastExtremes <- after_disturba$extreme[length(after_disturba$extreme)]
#                       extremes_at_disturbance <- xdata$extreme[sxz1]
#                       prsdcum_at_disturbance <- xdata$sdcum[sxz1]
#                       prsdcum_last <- after_disturba$sdcum[length(after_disturba$sdcum)]
#                       pixelCumsum_at_disturbance <- xdata$pixelCumsum[sxz1]
#                       pixelCumsum_last <- after_disturba$pixelCumsum[length(after_disturba$pixelCumsum)]
#                       prPnCV_at_disturbance <- xdata$PnCV[sxz1]
#                       prPnCV_last <- after_disturba$PnCV[length(after_disturba$PnCV)]
#                       patch_at_disturbance <- xdata$y[sxz1]
#                       patch_last <- after_disturba$y[length(after_disturba$y)]
#                       sdTrend <- lmModelSlope
#                       
#                       dav <- c(currentv,diff_distarbaDate_LastDate,sdTrend,patch_at_disturbance,patch_last, extremes_at_disturbance,lastExtremes,prsdcum_at_disturbance,
#                                prsdcum_last, pixelCumsum_at_disturbance,pixelCumsum_last,prPnCV_at_disturbance,prPnCV_last,
#                                last_segmet_slope_spatio,last_segmet_slope_pixel,secondLast_segmet_slope_vquotent,last_segmet_slope_vquotent)
#                       chan <- 0
#                       
#                       
#                     }else{
#                       
#                       dav <- as.numeric(samo)
#                       chan <- 0
#                     }
#                     
#                     
                  }
                  }else{
                    
                    dav <- as.numeric(samo)
                    chan <- 0
                  }
                  
                }else {
                  dav <- as.numeric(samo)
                  countx <- countx + 1
                }
              }
              
            }else{
              dav <- as.numeric(samo)
            }
          }else{
            dav <- as.numeric(samo)
          }
        }else{
          dav <- as.numeric(samo)
        }
        if (densityPlot == T){
          plot(density(dtaax))
          abline(v = xp, col = "blue", lty = 2)
          abline(v = qt, col = "red", lty = 2)
        }
      }else{
        dav <- as.numeric(samo)
      }
    }else{dav <- as.numeric(samo)}
    #print(dav)
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
  minimumRegrowthTime <- minimumRegrowthTime
  min_time_length_afterDisturbance <- min_time_length_afterDisturbance
  samo <- c(rep(NA, 19))
  if(tryCatchError){
    result <-  tryCatch (spatioTemporalMetricsExtractorxc_regrowth_monitor2(dFrvx,my_dates, mYear,spatiaNormPercentile,threshold, density,windowwidth,sPatioNormalixse,minimumRegrowthTime,min_time_length_afterDisturbance ),
                         error=function(e) as.numeric(samo),
                         warning=function(w) as.numeric(samo))
    if (inherits(result, c("error","warning"))){
      result <- as.numeric(samo)
    } else {
      result <- result
    }
  }else{
    
    result <- spatioTemporalMetricsExtractorxc_regrowth_monitor2(dFrvx,my_dates, mYear,spatiaNormPercentile,threshold, density,windowwidth,sPatioNormalixse,minimumRegrowthTime,min_time_length_afterDisturbance)
  }
  return(result)
}

