import os
import geospatialroutine.Routine as RT

"""
NAME:VectorDataLoad
AUTHOR: Saif Aati (saif@caltehc.edu)
PURPOSE:
    Loading into variables the displacement data from a correlation map (raster format) or text file.
    In order to draw a vector field.
INPUT:
    - ewFile: raster file contains the displacement in the x-direction
    - nsFile: raster file contains the displacement in the y-direction
    - ewBand: band number of the EW raster (default 1)
    - nsBand : band number of NS raster (default 1) 
    - txt: Set to 1 if the data are read from a txt file (0 default is raster file)
    - mapProj: Variable receiving the projection into which are expressed the
    coordinates.If data are from a text file, the projection is asked to be added by the user
    
OUTPUT:
    - ewVal: 2D-Array contains displacement in x-direction
    - nsVal : 2D-Array contains displacement in y-direction
    - xCoord : 2D-array contains x coordinates 
    - yCoord : 2D-Array contains y coordinates 
    - verifVal: This allows to check if the variables returned are good. 
    If a problem occurred during the process, this variable is set to 1.

"""


def VectorDataLoad_FromRaster(ewFile, nsFile, ewBand=1, nsBand=2, verifVal=0):
    ewFileInfo = RT.GetRasterInfo(inputRaster=ewFile)
    nsFileInfo = RT.GetRasterInfo(inputRaster=nsFile)
    ewVal = RT.ImageAsArray(ewFileInfo,ewBand)
    nsVal = RT.ImageAsArray(nsFileInfo,nsBand)
    if ewFileInfo.get("Dimension") == nsFileInfo.get("Dimension"):
        xCoord, yCoord = RT.XCoord_YCoord_as_Grid(rasterInfo=ewFileInfo)

    else:
        verifVal = 1
        print("The two rasters array must have the same size and same projection map")
        return

    return ewVal, nsVal, xCoord, yCoord


def VectorDataLoad_FromTxt(ewFile, nsFile, mapProj):
    "the fill must be of form: xCoord, yCoord, xValue, yValue ...."
    return


if __name__ == '__main__':
    path = "C:\\temp"
    ewPath = os.path.join(path, "EW_Cum.tif")
    nsPath = os.path.join(path, "NS_Cum.tif")
    VectorDataLoad_FromRaster(ewFile=ewPath, nsFile=nsPath)
