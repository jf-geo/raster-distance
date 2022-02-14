# Raster Distance
A python script for calculating the distance from a point for each cell in an elevation raster. Developed to create rough depth-of-field masks for QGIS 3D exports.  
  
Made available under a [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) license.

## Background
I use QGIS' 3D functionality to create mockup 3D visualisations before creating production versions in Blender. QGIS currently lacks the ability to produce depth-of-field effects with it's 3D viewing capabilities, so I approximate that using blurred layers and masks based on the distance between the raster and the camera position. This is the script I use to produce the 'distance from the camera' raster.

## Functionality
#### Current functionality:  
- Convert QGIS 3D coordinates to real world coordinates, so long as you know the false origin of the 3D view.  
- Calculate the distances between a real world point and each cell in an elevation raster.  
  
#### Potential future functionality:  
- Calculate the false origin of the 3D view based on base layers being used.  
- Process rasters strip by strip instead of loading everything into memory.  
- Calculate distance based on an orthographic view.  
- Calculate actual depth-of-field blurring based on camera parameters.  
- Wrap functionality into a QGIS processing algorithm.  
- Create plugin (maybe).  
  
#### Known bugs/issues:  
- Script produces an error when trying to transfer a NoDataValue from the input raster to the output raster.  
- Script loads the entire input raster band into memory and doesn't check whether or not that will cause a MemoryError.  
  
## Requirements  
Requires numpy and gdal python packages, which will be included by default if using a QGIS 3+ python environment.  
  
## Example  
To do.  
  
