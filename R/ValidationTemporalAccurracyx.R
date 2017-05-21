
#****************************************************************************
#TEMPORAL DELAY

#RETURN THE MEDIAN NUMBER OF OBSERVATIONS BEFORE CHANGED DETECTED
#*****************************************************************************

temporalDeley <- function(vdata,dtes, TmSeriesOfSamplePixels,  returnMedianDelay = T , nOfindices =1, colmWith =NULL){
  require(lubridate)
  nOfindices  <- nOfindices
  dtes <- dtes
  dt90 <- vdata
  tda <- TmSeriesOfSamplePixels
  xme <- c()
  tmo <- as.data.frame(matrix(nrow =1, ncol =2))
  xx0 <-  colmWith
  for (i in 1:nOfindices){
    xx0 <- xx0+2
    t4xx <- c()
    t5x <- c()
    for (i in 1: nrow(tda)){
      tc <- as.numeric(tda[i,])
      tc <- tc[1:length(tc)]
      tc2 <- as.data.frame(tc)
      tc2$date <- round (dtes, digits = 4)
      tx5 <- dt90[i, ]
      a <- tx5$ChangeDate
      if (!is.na(a)){
        b90 <- as.Date(as.character(a), format = "%Y%j")
        b80 <- decimal_date(b90)
        tx12 <- round(b80,digits =4)-0.002
        tx13 <- round (tx5[,xx0],digits = 4)+0.002
        tx12 <- tx12-0.002
        tx13 <- tx13+0.002
        tz2 <- subset (tc2, (tc2$date >= tx12 & tc2$date <= tx13))
        x1x <-  length(subset (tz2$tc, !is.na(tz2$tc)))
        x2x <- subset (tz2$date, !is.na(tz2$tc))
        if (x1x != 0){
          x1x1 <-  x1x - 1
        }else{
          x1x1 <- NA
          x1x2 <- NA
        }
      }else {
        x1x1 <- NA
      }
      t4xx <- c(t4xx, x1x1)
    }
    t5z <- replace (t4xx, t4xx == -1, 0 )
    tmo <- cbind(tmo,t5z)
    x <- subset(t5z, !is.na(t5z))
    xme <- c(xme, median(x))
  }

  if (returnMedianDelay){
    xmd <-  xme
  }else{ xmd <-  t4xx }

  return(xmd)

}
