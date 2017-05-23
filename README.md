# -------------------------STEF--------------------------------

## Welcome to STEF (Space-Time Extremes and Features) package! 
 
STEF is a R package that provides functionality for space-time forest change monitoring using observations from multiple satellites while exploiting both spatial and temporal information in satellite observations to identify forest disturbnaces. The STEF identifies potential forest disturbances as "extreme events" in local data cubes of satellite observations, and subsequently extract several spatio-temporal features once a pixel is flagged as potentially disturbed. With training data, the spatio-temporal features are then used to discriminate real forest disturbances from false detections. The STEF is developed for near real-time forest change monitoring, but it can also be used to map historcal forest disturbances. 



### Installing STEF from GitHub

```{r, eval=F, echo=T} 
#### first install devtools package 

install.packages("devtools")


#### Then install STEF package 

install_github("hamun001/STEF")

```