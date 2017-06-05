# -------------------------STEF--------------------------------

## Welcome to STEF (Space-Time Extremes and Features) package! 
 
STEF is a R package that provides functionality for space-time forest change monitoring using observations from multiple satellites. STEF exploits both spatial and temporal information in the time series of satellite observations to identify forest disturbnaces. Potential forest disturbances are identified as "extreme events" in local data cubes of satellite observations. Once potential forest disturbance is idenified at the pixel level, STEF extracts several spatio-temporal features. Using  Machine learning algorithm (e.g. Random forest),spatio-temporal features are the used calculate the probability of forest disturbances. A probability threshold is then used to discriminate forest disturbances from false detections. STEF is developed for near real-time forest change monitoring, but it can also be used to map historcal forest disturbances.STEF is developed specifically to facilitate robust:

1. Near real-time forest change detection in dry forest where strong seasonality in photosynthesis exist.
2. Multi-sensor data combination (E.g. Landsat, Sentinels and RapidEye)
3. Detection of small-scale and low-magnitude forest changes

# -------------------------How to use STEF--------------------------------


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

STEF detects forest disturbances using *stef_monitor* function, which is called using *rasterEngine* function from *spatial.tools* package. The change detection can be scaled to multiple cores or can be done sequentiallly. 


```{r, eval=F, echo=T}

## load the required packages

require("raster")
require("rgdal")
require("doParallel")
require("spatial.tools")
library(utils)
require("STEF")
require("randomForest")

# read the raster stack (NDVI raster stack,  30m resolution). It contains 183 NDVI layers, for 2013 through 2016, acquired by Landsat 7, 8 and Sentinel-2A sensors. Note that this test data set is pre-processed already, and is ready for analysis. Also note that the data set is already normalised spatially using the global spatial normalisation.

re_test <- brick("data/L7L8S2_30m_ndvi_kafa_Global_subset.tif")

# extract a subset from the test data set

re_testx <- crop(re_test, c(810000, 820000, 820000,830000))

# read the image acquistion dates. I save these in STEF/data directory

my_dates <- readRDS(file = "data/my_dates.rds")

#Plot the first layer in the raster stack

plot (raster(re_test, 1))

## We start the monitoring in 2016. Note this may take long, depending on your computer 


## sequential processing example:

rad <- rasterEngine(inraster=re_testx, fun=stef_monitor,window_dims=c(windowwidth=15,windowwidth =15),
                    args=list(mYear = 2016,density = F,my_dates =my_dates,threshold = 0.05,spatiaNormPercentile =95, windowwidth=9,tryCatchError=F))


## paralell processing example:

## register the cores

sfQuickInit(cpus=16)

rad <- rasterEngine(inraster=re_testx, fun=stef_monitor,window_dims=c(windowwidth=9,windowwidth =9),
        args=list(mYear = 2016,density = F,my_dates =my_dates,threshold = 0.05,spatiaNormPercentile =95, windowwidth=9,tryCatchError=T,sPatioNormalixse =F))
        
## unregister the cores

sfQuickStop()

# use random forest model to calculate the probability of forest disturbance 

## read the training data (this data set was used by Hamunyela et al (2017))

training_data <- readRDS(file = "data/trainingData.rds") # RC ==Real change; FC == False change

# set the seed and training the random forest model
set.seed(100)

rf_model <- randomForest(label ~ ., data=training_data,  ntree = 501, importance=TRUE,
                         proximity=TRUE, probability = T)


# Apply the trained random forest model to the entire image 

## Exclude the first layer (date of change)
radx <-subset(rad, c(2:nlayers(rad)))

## rename the layers to the colnames in the training data
names(radx) <- colnames(training_data)[1:17]

## do the prediction
m3 <- predict(radx,rf_model, type='prob',na.rm=TRUE, index=2)

## plot the computed probability map for forest disturbance
plot (m3)

##Combine the date of change layer with the probability map for forest disturbance 

cDate <-subset(rad, 1)

cMap <- stack(cDate,m3)
names(cMap) <- c("Date_of_forest_disturbance", "Probability_of_forest_disturbance")
plot(cMap)

```

## References 

1. Hamunyela, E., Verbesselt, J., Herold, M. (2016). Using spatial context to improve early detection of deforestationfrom Landsat time series. Remote Sensing of Environment,172, 126-138.\url{http://dx.doi.org/10.1016/j.rse.2015.11.006}

2. Hamunyela, E., Verbesselt, J.,de Bruin, S., Herold, M. (2016). Monitoring Deforestation at Sub-Annual Scales as Extreme Events in Landsat Data Cubes. Remote Sensing, 8(8), 651.\url{http://dx.doi.org/10.3390/rs8080651}

3. Hamunyela, E., Reiche, J., Verbesselt, J., & Herold, M. (2017). Using space-time features to improve detection of forest disturbances from Landsat time series. Remote Sensing, 9(6), 515.\url{http://dx.doi:10.3390/rs9060515}

4. Hamunyela, E., Reiche, J., Verbesselt, J., Tsendbazar, N.E., & Herold, M. Combining Sentinel-2 and Landsat time series for small-scale forest change monitoring. In prep.

5. Reiche, J., Hamunyela, E., Verbesselt, J., Hoekman, D., & Herold, M. Improving near-real time deforestation monitoring in tropical dry forests by combining dense Sentinel-1 time series with Landsat and ALOS-2 PALSAR-2. Remote Sensing of Environment. In review.
