# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2020
import os

import geojson
import numpy as np

from pathlib import Path


def RasterFootprint(rasterPath, z=None, demPath=None, writeFp=True, savingPath=None, debug=False):
    """

    Args:
        rasterPath:
        z:
        demPath:
        writeFp:
        savingPath:

    Returns:

    """
    import geoRoutines.georoutines as geoRT
    if debug:
        print("-------", rasterPath)
    rasterInfo = geoRT.RasterInfo(rasterPath, False)
    w = rasterInfo.rasterWidth
    h = rasterInfo.rasterHeight
    epsg = None
    if rasterInfo.validMapInfo == True:
        epsg = rasterInfo.EPSG_Code

    if epsg != None:
        rasterExtent, _ = rasterInfo.RasterDims()
        # print(rasterExtent)
        if debug:
            print("--- Footprint will be in WGS 84 projection system --- ")
        if epsg != 4326:
            x0Map, y0Map, _ = geoRT.ConvCoordMap1ToMap2(x=rasterExtent[0],
                                                        y=rasterExtent[2],
                                                        targetEPSG=4326, sourceEPSG=epsg)
            xfMap, yfMap, _ = geoRT.ConvCoordMap1ToMap2(x=rasterExtent[1],
                                                        y=rasterExtent[-1],
                                                        targetEPSG=4326, sourceEPSG=epsg)

            lats = x0Map, xfMap, xfMap, x0Map, x0Map
            lons = yfMap, yfMap, y0Map, y0Map, yfMap
        else:
            lats = rasterExtent[2], rasterExtent[2], rasterExtent[3], rasterExtent[3], rasterExtent[2]
            lons = rasterExtent[0], rasterExtent[1], rasterExtent[1], rasterExtent[0], rasterExtent[0]


    elif rasterInfo.rpcs:
        if debug:
            print("--- No geoTransform in input Raster, footprint will be calculated from RFM if exist ")

        try:

            # import RFM.cRFM as RFM
            import geoCosiCorr3D.geoRFM.cRFM as RFM

        except:
            print("Error ! Cant import RFM ")
            return
        rpc = RFM.cRFMModel(RFMFile=rasterPath)
        if z is None:
            if demPath:
                # from geospatialroutine.Temp.TP2GCPs import GroundCoord2Alt
                # z = GroundCoord2Alt(demPath=demPath, lon=rpc.lonOff, lat=rpc.latOff, interpType=2)
                ## TODO
                raise ValueError("Not Implemented")
            elif z == None:
                import warnings
                warnings.warn("No DEM available data, RPC altOff will be used")
                z = rpc.altOff
                if debug:
                    print("z=", z)

        lons, lats, _ = rpc.Img2Ground_RFM(col=[0, 0, w, w, 0],
                                           lin=[0, h, h, 0, 0],
                                           altIni=[z, z, z, z, z],
                                           normalized=False)
        # print(lons, lats)

    else:
        raise ValueError("No geoTransfrom and RPCs in input raster ")

    crs__ = {
        "type": "name",
        "properties": {
            "name": "EPSG:4326"
        }
    }
    footprint = geojson.Feature(geometry=geojson.Polygon([list(zip(lons, lats))]),
                                properties={"name": os.path.basename(rasterPath)}, crs=crs__)

    if writeFp:
        if savingPath:
            return footprint, WriteJson(features=footprint, outputFile=savingPath)

        else:

            return footprint, WriteJson(features=footprint,
                                        outputFile=os.path.join(os.path.dirname(rasterPath),
                                                                Path(rasterPath).stem + "_fp"))
    else:
        return footprint


def WriteJson(features, outputFile, driver="GeoJSON"):
    """

    Args:
        features:
        outputFile:
        driver:

    Returns:

    """

    outputFile = outputFile + VectorDrivers(driver)
    # print(features)
    with open(outputFile, 'w') as outfile:
        geojson.dump(features, outfile, sort_keys=True)
    return outputFile


def VectorDrivers(driver):
    if driver == "GeoJSON":
        return ".geojson"
    if driver == "KML":
        return ".kml"
    if driver == "ESRI Shapefile":
        return ".shp"


def FootprintExtract_Batch(inputdirectory, fileExtension=[".NTF"], fpSavingDirectory=None):
    """

    Args:
        inputdirectory:
        fileExtension:
        fpSavingDirectory:

    Returns:

    """
    from geoRoutines.FilesCommandRoutine import ExtractSubfiles
    filesList = ExtractSubfiles(inputdirectory=inputdirectory, fileExtension=fileExtension)
    print("Total number of files:", len(filesList))

    outputs = []
    errorList = []
    for index, rasterPath in enumerate(filesList):
        basename = Path(rasterPath).stem
        path = os.path.dirname(rasterPath)
        if fpSavingDirectory:
            fpOutput = os.path.join(fpSavingDirectory, basename + "_fp")
        else:
            fpOutput = os.path.join(path, basename + "_fp")

        # print("\n===== Extract Footprint:", index + 1, "/", len(filesList), "=====\n")
        # FootprintExtract(rasterPath=rasterPath, fpOutput=fpOutput)
        if os.path.exists(fpOutput + ".geojson"):
            # print("==> Image footprint exist")
            outputs.append(fpOutput + ".geojson")
        else:
            try:
                _, fpPath = RasterFootprint(rasterPath=rasterPath, writeFp=True, savingPath=fpOutput)
                outputs.append(fpPath)
            except Exception:
                errorList.append(rasterPath)
                pass

    print("errorImages:", len(errorList))
    return outputs


def fp_area(geometry):
    """
    :param fpPath:
    :return: area in km2
    """
    from geoRoutines.Routine_Lyr import GetBoundsGeojson
    x, y = GetBoundsGeojson(geojsonFile=geometry, conv2UTM=True)

    height = np.sqrt((x[0] - x[1]) ** 2 + (y[0] - y[1]) ** 2) / 1000
    width = np.sqrt((x[0] - x[3]) ** 2 + (y[0] - y[3]) ** 2) / 1000
    # print(height,width)
    area = height * width
    # print(area)
    return area


def ComputeFootprint(rsmModel, oProj, save=False, oFolder=None, fileName=None, demInfo=None,
                     rsmCorrectionArray=np.zeros((3, 3))):
    """

    Args:
        rsmModel:
        demInfo:
        rsmCorrectionArray:

    Returns:

    """
    import geoRoutines.georoutines as geoRT
    from geoCosiCorr3D.geoRSM.Pixel2GroundDirectModel import cPix2GroundDirectModel
    iXPixList = [0, rsmModel.nbCols - 1, 0, rsmModel.nbCols - 1]
    iYPixList = [0, 0, rsmModel.nbRows - 1, rsmModel.nbRows - 1]

    # resTemp_ = p_map(Pixel2GroundDirectModel.cPix2GroundDirectModel,
    #                  len(iXPixList) * [rsmModel],
    #                  iXPixList,
    #                  iYPixList,
    #                  len(iXPixList) * [rsmCorrectionArray],
    #                  len(iXPixList) * [demInfo.rasterPath],
    #                  num_cpus=len(iXPixList))
    geoCoordList = []

    for xVal, yVal in zip(iXPixList, iYPixList):
        # print("\n----- xVal:{},yVal:{}".format(xVal, yVal))
        if demInfo != None:
            pix2Ground_obj = cPix2GroundDirectModel(rsmModel=rsmModel,
                                                    xPix=xVal,
                                                    yPix=yVal,
                                                    rsmCorrectionArray=rsmCorrectionArray,
                                                    demFile=demInfo.rasterPath)
        else:
            pix2Ground_obj = cPix2GroundDirectModel(rsmModel=rsmModel,
                                                    xPix=xVal,
                                                    yPix=yVal,
                                                    rsmCorrectionArray=rsmCorrectionArray,
                                                    demFile=None)

        geoCoordList.append(pix2Ground_obj.geoCoords)

    ## Convert foot printcoord to the grid projection system

    geoGround = np.asarray(geoCoordList)
    if oProj == 4326:
        crs__ = {
            "type": "name",
            "properties": {
                "name": "EPSG:4326"
            }
        }

        res = geoGround
        lons, lats = res[:, 0], res[:, 1]
        lons = [res[0, 0], res[1, 0], res[3, 0], res[2, 0], res[0, 0]]
        lats = [res[0, 1], res[1, 1], res[3, 1], res[2, 1], res[0, 1]]
        footprint = geojson.Feature(geometry=geojson.Polygon([list(zip(lons, lats))]),
                                    properties={"name": fileName}, crs=crs__)

        WriteJson(features=footprint, outputFile=os.path.join(oFolder, fileName))

        return geoGround, footprint,os.path.join(oFolder, fileName+".geojson")

    utmGround = geoRT.ConvCoordMap1ToMap2_Batch(X=list(geoGround[:, 1]),
                                                Y=list(geoGround[:, 0]),
                                                Z=list(geoGround[:, -1]),
                                                targetEPSG=oProj)
    topLeftGround = [utmGround[0][0], utmGround[1][0], utmGround[2][0]]
    topRightGround = [utmGround[0][1], utmGround[1][1], utmGround[2][1]]
    bottomLeftGround = [utmGround[0][2], utmGround[1][2], utmGround[2][2]]
    bottomRightGround = [utmGround[0][3], utmGround[1][3], utmGround[2][3]]

    return topLeftGround, topRightGround, bottomLeftGround, bottomRightGround, iXPixList, iYPixList
