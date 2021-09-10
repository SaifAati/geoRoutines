import os

import numpy as np

import Interpolation
import Main_routine as RT


def DMS2DD(d, m, s):
    dd = np.abs(d) + m / 60 + s / 3600
    if d < 0:
        return -1 * dd
    else:
        return dd


def GroundCoord2Alt(demPath, lon, lat, EPSG=4326, interpType=2):
    demInfo = RT.GetRasterInfo(demPath)

    if EPSG == int(demInfo.get("MapInfo")[2]):
        xPix, yPix = RT.Map2Pixel(imageInfo=demInfo, x=lon, y=lat)
    else:
        print("Convert to the same coordinate system from %s to %s : " % (EPSG, demInfo.get("MapInfo")[2]))
        RT.ConvCoordMap1ToMap2(x=lon, y=lat, targetEPSG=demInfo.get("MapInfo")[2])

    inteMeth = Interpolation.InterpolationMethod(method=interpType)
    print(inteMeth)

    h = Interpolation.InterpolationSA(x=xPix, y=yPix, demInfo=demInfo,
                                      interpolationType=inteMeth)

    return h


# def TiePtsAltitude(demPath, tiePoint, method, interpType):
#     """
#     # Verify if the window around the tie point is onside the DEM image
#     # The array could be outside the DEM image in case the DEM does not cover all the ROI
#     # Set h=0 in case outside, otherwise interpolate the DEM at the location of the Tie point
#     :param tiePoint: Tie point coordinates in map coordinate system of the DEM image : array
#     :param demData : list of the DEM raster information : raster
#     :param method  : int
#             :1= Cosi-Corr method : using window of 7x7 around Floor pixel
#             :2= New interpolation method using non integer-pixel
#     :param interpolationType : int
#             : 0= "nearest"    zero order interpolation
#             : 1= "bilinear"
#             : 2= "bicubic"
#             : 3= "spline"   :not yet implemented
#             : Default = "bilinear"
#     :return: the height of the tie points : List
#     """
#     ## if a DEM is used, retrieve for each tie point the altitude
#     # Verify if the DEM is valid
#
#     interpolationType_ = Interpolation.InterpolationMethod(method=interpType)
#
#     ## Get DEM information
#     demInfo = RT.GetLayerInfo(imagePath=demPath)
#
#     ## Get the number of tie points
#     nbTiepts = np.shape(tiePoint)[0]
#
#     altitude = []
#     ## Loop through all tie points
#     for i in range(nbTiepts):
#         ## Get the pixel of the tie points in the DEM image (from MapToPixel)
#         xTemp, yTemp = RT.Map2Pixel_(imagePath=self.demImgPath, x=tiePoint[i, 0], y=tiePoint[i, 1])
#         print("xTemp,yTemp", xTemp, yTemp)
#         if method == 1:
#             #### Interpolation as implemented by Francois in Cosi-Corr
#             h = self.InterpolationCosiCorr(xTemp=xTemp, yTemp=yTemp, demInfo=demInfo,
#                                            tiePointIndex=i, interpolationType=interpType)
#         if method = =2:
#             #### SA@ interpolation without defineing window around Floor pixel
#             h= self.InterpolationSA(xTemp=xTemp,yTemp=yTemp,demI nfo=demInfo,
#             tiePointIndex = i,  interpolationType = interpolationType_)
#
#             altitude.append(h)
#
#
#     return (altitude)


if __name__ == '__main__':
    path = "D:\\temp\RPCTest\WV"
    rasterPath = os.path.join(path, "20JAN21213736-P1BS-504089615060_01_P005.NTF")
    planetRaster = "D:\\temp\RPCTest\Planet\\20190702_171320_1052_1B_Analytic.tif"
    print(rasterPath)
    # GCPs(rasterPath)
    demRaster = "D:\\temp\RPCTest\WV\\n35_w118_1arc_v3.tif"
    demRasterConv = "D:\\temp\RPCTest\WV\dem_conv.pix"
    # Convert_Ground2Img(rasterPath=rasterPath, lon=, lat=35.5227778,demPath=demRaster)
    srtm90 = "D:\\temp\srtm_13_05.tif"

    lon = DMS2DD(-117, 35, 0.5)  # -117.3711111
    lat = DMS2DD(35, 40, 44.75)  # 35.5227778
    print(lon, lat)
    GroundCoord2Alt(demPath=srtm90, lon=-117.58347222222221, lat=35.67909722222222,interpType=1
                    )

    import rasterio
    dataset = rasterio.open(srtm90)
    print(dataset)
    lon = -117.58347222222221
    lat = 35.67909722222222
    py, px = dataset.index(lon, lat)
    print('Pixel Y, X coords: {}, {}'.format(py, px))
    N=1
    window = rasterio.windows.Window(px , py , N, N)
    clip = dataset.read(window=window)
    print(clip)
    print(window)


