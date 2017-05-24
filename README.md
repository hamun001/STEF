# -------------------------STEF--------------------------------

## Welcome to STEF (Space-Time Extremes and Features) package! 
 
STEF is a R package that provides functionality for space-time forest change monitoring using observations from multiple satellites. STEF exploits both spatial and temporal information in the time series of satellite observations to identify forest disturbnaces. Potential forest disturbances are identified as "extreme events" in local data cubes of satellite observations. Once potential forest disturbance is idenified at pixel level, STEF extracts several spatio-temporal features. With training data, the user can then use the machine learning algorithm (e.g. Random forest) to discriminate real forest disturbances from false detections using the spatio-temporal features. STEF is developed for near real-time forest change monitoring, but it can also be used to map historcal forest disturbances. 



### Installing STEF from GitHub

```{r, eval=F, echo=T} 
#### First install devtools package

install.packages("devtools")


#### Then install STEF package 

install_github("hamun001/STEF")

```

## STEF functionality

STEF provides the functionality for reducing seasonality and inter-sensor differences in satellite image time series. Seasonality and inter-sensor differences are reduced through spatial normaliastion.  STEF provides two spatial normalisation approach: Local and global. With local normalisation, a value of focal pixel is normalised using information derived from its neighbourhood. The size of the neigbourhood is defined by the user. The global normalisation uses information derived over entire image mosaic or single scence. STEF also provides functionality of monitoring forest change, and the functionality for calculating spatial accuracy is included (thanks to Dr Nandika Tsendbazar). Below some examples on how to use STEF's functionality are shown 

### Spatial normalisation to reduce seasonal variations 



```{r, eval=F, echo=T} 

#### read a raster stack
ras <- brick()

#### apply global spatial normalisation

stef_global_spatial_normaliser(ras, isStack = T, xpercentile = 0.95,output_filename ="ra_global_normalised.tif")

```
