


#BFASTMONITOR SET_UP
ybfastmonitor <- function(inraster,myear = myear,plot = F,history = c("all"),my_dates =my_dates,
                          minumum_observations = minumum_observations,magThreshold = magThreshold,type =type) {
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
    if (length(historyPeriod2$dates) > 5){
      tim1 <- c(round (decimal_date(dtat$dates), digits =0))
      tim2 <- unique(tim1)
      timex <- subset (tim2, tim2 >= myear)
      coun <- 1
      breakpointx  <- NA
      while (is.na(breakpointx) & coun < length(timex)){
        bpts <- bfastts(dtat$proCell,  dtat$dates , type = c("irregular"))
        stmon <- timex[coun]
        bfm <- bfastmonitor(data = bpts, formula = response ~ harmon, start=c(stmon, 1),plot = plot, history = history,type =type)
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
      magnitudex <- -88888
    }

  }else {
    breakpointx <- NA
    magnitudex <- -99999
  }

  output <- c(as.numeric(breakpointx),as.numeric(magnitudex))
  print(output)
  return(output)
}






