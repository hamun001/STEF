# -------------------------STEF--------------------------------

## Welcome to STEF (Space-Time Extremes and Features) package! 
 
STEF is a R package that provides functionality for space-time forest change monitoring using observations from multiple satellites while exploiting both spatial and temporal information in satellite observations to identify forest disturbnaces. STEF identifies potential forest disturbances as "extreme events" in local data cubes of satellite observations, and subsequently extract several spatio-temporal features when a extreme (negative anomalies ) are detected in the pixel-time series. With training data, the spatio-temporal features are then used to discriminate real forest disturbances from false detections. The STEF is developed for near real-time forest change monitoring, but it can also be used to map historcal forest disturbances. 



### Installing STEF from GitHub

```{r, eval=F, echo=T} 
#### First install devtools package

install.packages("devtools")


#### Then install STEF package 

install_github("hamun001/STEF")

```

### Spatial normalisation to reduce seasonal variations 


```{r, eval=F, echo=T} 

#### read a raster stack
ras <- brick()

#### apply global spatial normalisation

stef_global_spatial_normaliser(ras, isStack = T, xpercentile = 0.95,output_filename ="ra_global_normalised.tif")

```
