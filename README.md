# -------------------------STEF--------------------------------

## Welcome to STEF (Space-Time Extremes and Features) package! 
 
STEF is a R package that provides functionality for space-time forest change monitoring using observations from multiple satellites. STEF exploits both spatial and temporal information in the time series of satellite observations to identify forest disturbnaces. Potential forest disturbances are identified as "extreme events" in local data cubes of satellite observations. Once potential forest disturbance is idenified at pixel level, STEF extracts several spatio-temporal features. With training data, the user can then use the machine learning algorithm (e.g. Random forest) to discriminate real forest disturbances from false detections using the spatio-temporal features. STEF is developed for near real-time forest change monitoring, but it can also be used to map historcal forest disturbances. 



### Installing STEF from GitHub

```{r, eval=F, echo=T}

## First install devtools package

install.packages("devtools")


## Then install STEF package 

install_github("hamun001/STEF")

```

## STEF's functionality

In addition to forest change detection, STEF provides the functionality for reducing seasonality and inter-sensor differences in satellite image time series. Seasonality and inter-sensor differences are reduced through spatial normaliastion.  STEF provides two spatial normalisation approaches: Local and global. With local normalisation, a value of each pixel in the image is normalised using information derived from its neighbourhood. The size of the neigbourhood is defined by the user. The global normalisation uses information derived over entire image mosaic or single scence. STEF also provides functionality of monitoring forest change, and the functionality for calculating spatial accuracy is included (thanks to Dr Nandika Tsendbazar).  Examples on how to use STEF's functionality are shown below.

### Using STEF to reduce seasonal variations and inter-sensor differences in image time series 

#### 1. Local spatial normalisation 

Local spatial normalisation is recommended for areas where deciduous and evergreen forests co-exist. If you just want to use STEF to reduce seasonality and inter-sensor differences in the time series, and use other another change detection algorithm to detect forest disturbances, you should use *stef_local_spatial_normaliser* function. This function is applied over entire image using *rasterEngine* function  from *spatial.tools* package. If you want to use local spatial normalisation and detect forest changes using STEF, there is no need to do the normalisation separately; *stef_monitor* function which does change detection performs local spatial normalisation internally; you just specify that local normalisation should be done.  

```{r, eval=F, echo=T} 

## load the required packages (make sure these packages are installed)

require("spatial.tools")
require("raster")
require("rgdal")
require("doParallel")
require("STEF")

## read a raster stack

ras <- brick()


## sequential processing example:

rad <- rasterEngine(inraster=rasterBrick, fun=stef_local_spatial_normaliser,window_dims=c(windowwidth=15,windowwidth =15),
                   args=list(spatiaNormPercentile =95))
                   
                   
## paralell processing example:

## register the cores

sfQuickInit(cpus=5)
rad <- rasterEngine(inraster=rasterBrick, fun=stef_local_spatial_normaliser,window_dims=c(windowwidth=15,windowwidth =15),
                   args=list(spatiaNormPercentile =95))
                   
## Unregister the cores

sfQuickStop()

```

#### 2. Global spatial normalisation 

Global spatial normalisation is recommended for areas where  deciduous and evergreen forests do not co-exist. It is much fast than local normalisation, and it should be applied separately even when you want to use STEF for change detection. *stef_global_spatial_normaliser* function which does global normalisation can be applied to a single image or a stack of images; 

```{r, eval=F, echo=T}

## load the required packages 

require("raster")
require("rgdal")
require("doParallel")
require("spatial.tools")
require("STEF")

## read a raster stack

ras <- brick()

## apply global spatial normalisation

stef_global_spatial_normaliser(ras, isStack = T, xpercentile = 0.95,output_filename ="ra_global_normalised.tif")

```
### Monitoring forest changes using STEF 

STEF detects forest disturbances using *stef_monitor* function, which is called using *rasterEngine* function from *spatial.tools* package. The change detection can be scaled to many cores or can be done sequentiallly. 


```{r, eval=F, echo=T}

## load the required packages

require("raster")
require("rgdal")
require("doParallel")
require("spatial.tools")
require("STEF")

## sequential processing example:

rad <- rasterEngine(inraster=rasterBrick, fun=stef_monitor,window_dims=c(windowwidth=15,windowwidth =15),
                    args=list(mYear = 2014,density = F,my_dates =imagedate,threshold = 0.01,spatiaNormPercentile =95, windowwidth=15,tryCatchError=T))

## paralell processing example:

## register the cores

sfQuickInit(cpus=5)

rad <- rasterEngine(inraster=ra, fun=stef_monitor,window_dims=c(windowwidth=15,windowwidth =15),
        args=list(mYear = 2014,density = F,my_dates =imagedate,threshold = 0.01,spatiaNormPercentile =95, windowwidth=15,tryCatchError=T,sPatioNormalixse =T))
        
## unregister the cores

sfQuickStop()


```

## References 

1. Hamunyela, E., Verbesselt, J., Herold, M. (2016). Using spatial context to improve early detection of deforestationfrom Landsat time series. Remote Sensing of Environment,172, 126-138.\url{http://dx.doi.org/10.1016/j.rse.2015.11.006}

2. Hamunyela, E., Verbesselt, J.,de Bruin, S., Herold, M. (2016). Monitoring Deforestation at Sub-Annual Scales as Extreme Events in Landsat Data Cubes. Remote Sensing, 8(8), 651.\url{http://dx.doi.org/10.3390/rs8080651}

3. Hamunyela, E., Reiche, J., Verbesselt, J., & Herold, M. (2017). Using space-time features to improve detection of forest disturbances from Landsat time series. Remote Sensing, 9(6), 515.\url{http://dx.doi:10.3390/rs9060515}

4. Hamunyela, E., Reiche, J., Verbesselt, J., Tsendbazar, N.E., & Herold, M. Combining Sentinel-2 and Landsat time series for small-scale forest change monitoring. In Prep.

5. Reiche, J., Hamunyela, E., Verbesselt, J., Hoekman, D., & Herold, M. Improving near-real time deforestation monitoring in tropical dry forests by combining dense Sentinel-1 time series with Landsat and ALOS-2 PALSAR-2. Remote Sensing of Environment. In review.
