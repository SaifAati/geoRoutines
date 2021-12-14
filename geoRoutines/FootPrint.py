# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2020
import os

import geojson
import numpy as np

from pathlib import Path


def RasterFootprint(rasterPath, z=None, demPath=None, writeFp=True, savingPath=None):
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
                print("z=", z)

        lons, lats = rpc.Img2Ground_RFM(col=[0, 0, w, w, 0],
                                        lin=[0, h, h, 0, 0],
                                        alt=[z, z, z, z, z],
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
