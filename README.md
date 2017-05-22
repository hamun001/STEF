# stFeatures

## Welcome to stFeatures package. 
 
stFeatures package provides functionality for space-time forest change monitoring from satellite image time series, allowing for concurrent exploittion of spatial and temporal information in satellite observations to identify forest disturbnaces. The algorithm identifies potential forest disturbances as "extreme events" in local data cubes of satellite observations. Once the pixel is flagged as potentially disturbed, several space-time features are extracted from the local data cube. With training data,  the space-time features are used to discriminate real forest disturbances from false detections. The algorithm is developed for near real-time forest change monitoring, but it can also be used to map historcal forest disturbances.