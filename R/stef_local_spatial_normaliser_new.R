stef_local_spatial_normaliser_new <- function(TS, window_size=c(9, 9), spatiaNormPercentile=95, CPUs, ...){
  
  # local spatial normalisation
  parallel_local_spatial_normalisation <- function(inraster, spatiaNormPercentile=95, ...){
    
    # local spatial neighbourhood window as vector
    dFrv <- as.vector(inraster)
    
    # pull out central ts value of the spatial window
    proC <- dFrv[length(dFrv) - trunc(length(dFrv) / 2)]
    
    # check if value
    if(!is.na(proC)){
      
      # apply spatial normalisation over all ts values in the spatial neighbourhood window
      dFrv <- spNormaliser(dFrv)
      
      # pull out the normalised ts values of the grid cell that is located right in the middle of the spatial neighbourhood window 
      proCell <- dFrv[length(dFrv) - trunc(length(dFrv) / 2)]
      
    }else{proCell <- proC}
    
    # return normalised central ts value
    return(proCell)
    
  }
  
  # Sampling & normalisation function
  spNormaliser <- function(x){
    
    # sample the user defined percentile over all values in the spatial neighbourhood window 
    y <- round(as.numeric(quantile(x, c(spatiaNormPercentile / 100), na.rm=T)), digits=2)
    
    # normalise the values in the spatial neighbourhood window  
    x <- as.numeric(x) / y
    return(x)
    
  }
  
  # register the cores
  sfQuickInit(cpus=CPUs)
  
  for(i in 1:nlayers(TS)){
    
    # run local spatial normalisation in parallel mode
    t <- system.time(outraster <- rasterEngine(inraster=TS[[i]], fun=parallel_local_spatial_normalisation, window_dims=c(window_size[1], window_size[2]), args=list(spatiaNormPercentile=95)))
    
    cat(as.numeric(t[3]) / 60, "min elapsed\n")
    
    # re-build normalised raster stack
    if(i == 1){
      
      TSn <- stack(outraster)
      
    }else{
      
      TSn <- addLayer(TSn, outraster)
      
    }
    
    cat("layer", i, "normalised & re-stacked\n")
    
  }
  
  # unregister the cores
  sfQuickStop()
  return(TSn)
  
}