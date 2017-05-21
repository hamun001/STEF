
spatioTemporalMetricsExtractor <- function(inraster,proC, mYear, a01,spatiaNormPercentile,threshold, density,windowwidth,...){

  require ("raster")
  require("rgdal")
  require("lubridate")
  require("zoo")
  require("bfast")
  require("doParallel")


  print("Tested here")
  windowwidth <- windowwidth
  qdates <- my_dates
  a <- brick(inraster)
  dFrv <- as.data.frame(t(getValues(a)))
  is.na(dFrv) <- do.call(cbind,lapply(dFrv, is.nan))
  is.na(dFrv) <- do.call(cbind,lapply(dFrv, is.infinite))
  proC<- (as.numeric(dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)]))
  a01 <- as.data.frame(proC)
  a01$date <- qdates
  a01 <- subset(a01, a01$date > mYear & !is.na(a01$proC) & !is.nan(a01$proC) & !is.infinite(a01$proC))
  xp <- NA
  CH <- NA
  dav <- as.numeric(c(NA, NA,NA,NA,NA, NA, NA,NA, NA,NA,NA,NA))

  if (length (a01$date) > 0){
    rownames(dFrv) <- qdates
    dFrv$deci_date <- qdates
    de_date <- dFrv$deci_date
    dFrv$deci_date <- NA
    spNormaliser <- function(x ){
      y <- round(as.numeric(quantile(x, c(spatiaNormPercentile/100), na.rm =T)),digits = 2)
      #z <- median(subset (x, x >= y), na.rm =T)
      x <- as.numeric(x)/y
      return(x)
    }
    print ("here")
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

    dFrv <- apply(dFrv, 1, spNormaliser)
    dFrv <- t(dFrv)
    proCell<- (as.numeric(dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)]))

    #pull out the neighbours of focal pixel
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

    # check if one or more neigbours are non-forest (masked)
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
    print(paste("boo:_:", nonforest, sep=""))
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
    print(paste("shooo:_:",length(atx$x ), sep=""))
    if (length(atx$x ) > 0){
      ata <- subset(xdata, !is.na(xdata$x))
      cdata <- subset(xdata, is.na(xdata$x))
      xdata <- rbind(ata, cdata)
      xdata <- xdata[order(xdata$deci_date),]
      print("kim")
      dFrv[,ncol(dFrv)-trunc(ncol(dFrv)/2)] <- removedips(xdata$x)
      print("Eliakim")
      dFrv$deci_date <- xdata$deci_date
      dFrv <-subset (dFrv, round(dFrv$deci_date ,digits = 4) <  mYear)
      dFrv$deci_date <- NA

      dtaa <- as.numeric((unlist(dFrv)))
      rm(dFrv)
      dtaax <- subset(dtaa, !is.na(dtaa) & !is.nan(dtaa) & dtaa > 0)
      rm(dtaa)
      qt <- round(as.numeric(quantile(dtaax, c(threshold), na.rm =T)),digits = 4)
      vq <- round(sd(dtaax, na.rm =T), digits = 4)
      #qt <- qt - abs(optimalmagnitude)
      pacths <- function(x,qta){
        y1 <- subset(x, !is.na(x) & x < qta)
        b2 <- length(y1)
        return(b2)
      }
      dpatchx <- apply(neig, 1, pacths, qt)
      xdata$y <- as.numeric( dpatchx)
      vqs <- length(dtaax)
      xdatap <- subset (xdata, !is.na(xdata$x) & !is.nan(xdata$x) & xdata$x < qt)
      #xdata <- subset (xdata, !is.na(xdata$x) & xdata$x < qt)
      #xdatap <- subset (xdata, !is.na(xdata$x) & xdata$x < qt )
      xdatax <- subset (xdatap$deci_date , xdatap$deci_date >=  mYear)
      print(paste("pooo:_:",length(xdatap$x), sep=""))
      if (length(xdatap$x) > 0 & length(xdatax) !=0){
        xdata <- subset (xdata, !is.na(xdata$x) & !is.nan(xdata$x ))
        xdat1 <- subset (xdata , xdata$deci_date <  mYear & xdata$x > qt )
        print(paste("poooMom:_:",length(xdat1$x), sep=""))
        if (length(xdat1$x) > 0){
          xdato <- subset (xdata$deci_date , round(xdata$deci_date, digits = 4) < mYear)
          xdata <- xdata[length(xdato):length(xdata$x),]
          print(paste("pooono:_:",length(xdata$x), sep=""))
          if (length(xdata$x) > 1){
            countx <- 1
            chan <- NA
            while ( is.na(chan) & countx < (length(xdata$x) - 2)){
              sxz1 <- countx+1
              #sxz2 <- countx+2
              #if (xdata$x[countx] < qt & xdata$y[countx] >= patchsize & xdata$x[sxz1] < qt & xdata$x[sxz2] < qt ) {# | xdata$x[countx] < qt & xdata$y[sxz1] >= patchsize & xdata$x[sxz1] < qt & xdata$x[sxz2] < qt ){
              print(paste("nooono:_:",xdata$x[countx],xdata$x[sxz1], sep=""))
              if (xdata$x[countx] < qt & xdata$x[sxz1] < qt){
                currentv <- round(xdata$deci_date[sxz1], digits = 4)
                xp <- round(xdata$x[sxz1],digits = 4)
                prePatch <- xdata$y[countx]
                pxPatch <- xdata$y[sxz1]
                print(paste("peee:_:",prePatch, sep=""))
                if (prePatch !=0){
                  NboursStep1 <- 1
                  prePatch <-  prePatch
                }else{NboursStep1 <- 0}
                print(paste("beee:_:",pxPatch, sep=""))

                if (pxPatch != 0){
                  pxNboursStep2 <- 1
                  pxPatch <-  pxPatch
                }else{pxNboursStep2 <- 0}
                vt8 <- xp -qt
                CH <- 2
                dav <- as.numeric(c(currentv,vqs,vq,qt,vt8,CH,prePatch,NboursStep1,pxPatch, pxNboursStep2,nonforestNeig,NOofnonforestNeig))
                chan <- 0
              }else {
                currentv <- NA
                dav <- as.numeric(c(currentv,NA,NA,NA,NA, NA, NA,NA, NA,NA,NA,NA))
                countx <- countx + 1
              }
            }

            print(paste("keeoo:_:",dav[1], sep=""))
          if (is.na(dav[1])){
            print(paste("meeoo:_:",xdata$x[length(xdata$x)], sep=""))
            if(xdata$x[length(xdata$x)] < qt) {
              currentv <- round(xdata$deci_date[length(xdata$x)], digits = 4)
              xp <- round(xdata$x[length(xdata$x)],digits = 4)
              prePatch <-xdata$y[(length(xdata$x)-1)]
              pxPatch <- xdata$y[length(xdata$x)]
              print(paste("kuuoo:_:",prePatch, sep=""))
              if (prePatch != 0){
                NboursStep1 <- 1
                prePatch <-  prePatch
              }else{NboursStep1 <- 0}
              print(paste("kiioo:_:",pxPatch, sep=""))

              if (pxPatch != 0){
                pxNboursStep2 <- 1
                pxPatch <-  pxPatch
              }else{pxNboursStep2 <- 0}
              vt8 <- xp -qt
              CH <- 1
              dav <- as.numeric(c(currentv,vqs,vq,qt,vt8, CH,prePatch,NboursStep1,pxPatch, pxNboursStep2,nonforestNeig,NOofnonforestNeig))
            }else{
              dav <- as.numeric(c(NA, NA,NA,NA,NA, NA, NA,NA, NA,NA,NA,NA))
            }
          }
          }else{
            currentv <- NA
            dav <- as.numeric(c(NA, NA,NA,NA,NA, NA, NA,NA, NA,NA,NA,NA))
          }
        }else{
          currentv <- NA
          dav <- as.numeric(c(NA, NA,NA,NA,NA, NA, NA,NA, NA,NA,NA,NA))
        }
      }else{
        currentv <- NA
        dav <- as.numeric(c(NA, NA,NA,NA,NA, NA, NA,NA, NA,NA,NA,NA))
      }
      if (density == T){
        plot(density(dtaax))
        abline(v = xp, col = "blue", lty = 2)
        abline(v = qt, col = "red", lty = 2)
      }
    }else{
      dav <- as.numeric(c(NA, NA,NA,NA,NA, NA, NA,NA, NA,NA,NA,NA))
    }
  }else{dav <- as.numeric(c(NA, NA,NA,NA,NA, NA, NA,NA, NA,NA,NA,NA))}
  # return change date, number of obs, std, qt, Diff,
  #No. of consercative temporal extremes (two or one),
  #No of neigbour changing at t1
  #print(dav)
  return (dav)
  }


