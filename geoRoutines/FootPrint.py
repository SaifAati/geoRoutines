# Copyright (C) 2020, SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
import os

import geojson
import numpy as np

import geoRoutines.RoutineWarnings as RW
from osgeo import gdal, ogr
from pathlib import Path

"""
# Say you have an 8-bit image with a NoData value of 0,
#  use gdal_translate to stretch values [1,255] to [1,1], 
and assign 0 to represent the NoData value on output,
#  save result as VRT rather than writing a new raster to disk.
gdal_translate -scale 1 255 1 1 -ot Byte -of vrt -a_nodata 0 input_ortho.tif input_ortho_mask.vrt

# Create a polygon shapefile that outlines the valid (DN==1) regions of the VRT just created
gdal_polygonize.py -8 input_ortho_mask.vrt -f "ESRI Shapefile" mask_footprint.shp mask_footprint DN
"""


def RasterFootprint(rasterPath, z=None, demPath=None, writeFp=True, savingPath=None):
    """

    :param rasterPath:
    :param z:
    :param demPath:
    :param verbose:
    :param writeFp:
    :param savingPath:
    :return:
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
        print("--- Raster not gorefrenced foot-print will be calculated from RPC if exist ")

        try:
            import RFM.cRFM as RFM

        except:
            print("Error ! Cant import RFM ")
            return
        rpc = RFM.cRFMModel(RFMFile=rasterPath)
        if z is None:
            if demPath:
                from geospatialroutine.Temp.TP2GCPs import GroundCoord2Alt
                z = GroundCoord2Alt(demPath=demPath, lon=rpc.lonOff, lat=rpc.latOff, interpType=2)

            elif z == None:
                import warnings
                warnings.warn("no DEM available data, RPC altOff will be used",
                              category=RW.NoDEMWarning)
                z = rpc.altOff
                print("z=", z)

        lons, lats = rpc.Img2Ground_RFM(col=[0, 0, w, w, 0], lin=[0, h, h, 0, 0], alt=[z, z, z, z, z],
                                        normalized=False)
        # print(lons, lats)

    else:
        raise RW.NotGeoreferencedError

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

    :param features:
    :param outputFile:
    :param driver:
    :return:
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


def ExtractSubfiles(inputdirectory, fileExtension=[".NTF"]):
    filesList = []
    for root, dirs, files in os.walk(inputdirectory):
        for name in files:
            if any(name.endswith(ele) for ele in fileExtension):
                file = os.path.join(root, name)
                filesList.append(file)
    return filesList


def FootprintExtract(rasterPath, fpOutput, driver="GeoJSON"):
    """
    Compute the longitude, latitude footprint of an image using its RPC model.

    Args:
        geotiff_path (str): path or url to a GeoTIFF file
        z (float): altitude (in meters above the WGS84 ellipsoid) used to
            convert the image corners pixel coordinates into longitude, latitude

    Returns:
        geojson.Feature object containing the image footprint polygon
    """
    vrtTemp = os.path.join(path, "TempVrt.vrt")
    gdal.Warp(vrtTemp, rasterPath, format="vrt", dstNodata=0, dstAlpha=True)
    # dst_layername = os.path.join(path,"POLYGONIZED_STUFF")
    dst_layername = fpOutput
    drv = ogr.GetDriverByName(driver)
    dst_ds = drv.CreateDataSource(dst_layername + VectorDrivers(driver))
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=None)
    raster = gdal.Open(vrtTemp)
    srcband = raster.GetRasterBand(2)
    gdal.Polygonize(srcband, None, dst_layer, -2, [], callback=None)
    os.remove(vrtTemp)
    return


def FootprintExtract_v2(rasterPath, fpOutput, driver="GeoJSON"):
    # vrtTemp = os.path.join(path, "TempVrt.vrt")
    # gdal.Warp(vrtTemp, rasterPath, format="vrt", dstNodata=0, dstAlpha=True)
    # # dst_layername = os.path.join(path,"POLYGONIZED_STUFF")
    # dst_layername = fpOutput
    # drv = ogr.GetDriverByName(driver)
    # dst_ds = drv.CreateDataSource(dst_layername + VectorDrivers(driver))
    # dst_layer = dst_ds.CreateLayer(dst_layername, srs=None)
    # raster = gdal.Open(vrtTemp)
    # srcband = raster.GetRasterBand(2)
    # gdal.Polygonize(srcband, None, dst_layer, -2, [], callback=None)
    # os.remove(vrtTemp)

    raster = gdal.Open(rasterPath)
    RT.RasterExtent(rasterInfo=RT.GetRasterInfo(rasterPath))
    return


def FootprintExtract_Batch(inputdirectory, fileExtension=[".NTF"], fpSavingDirectory=None):
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

    x, y = RT.GetBoundsGeojson(geojsonFile=geometry, conv2UTM=True)

    height = np.sqrt((x[0] - x[1]) ** 2 + (y[0] - y[1]) ** 2) / 1000
    width = np.sqrt((x[0] - x[3]) ** 2 + (y[0] - y[3]) ** 2) / 1000
    # print(height,width)
    area = height * width
    # print(area)
    return area


if __name__ == '__main__':
    path = "G:\SkySatData\Ridgecrest\s104_20200620T182032Z\\basic_analytic"
    path = "G:\SkySatData\Morenci_Mine_AZ\Morenci_Mine_AZ\L1B\MUL_basic_analytic_converted\Converted"
    path = "G:\Ridgecrest\\0-Raw DATA\Zipped_Data_WV_Chris\Images_unzipped_stored\Pre_Eq"
    PATH = "G:\Ridgecrest\\0-Raw DATA\Zipped_Data_WV_Chris\Images_unzipped_stored\Post_eq"
    path = "G:\Ridgecrest\\0-Raw DATA\Zipped_Data_WV_Chris\Images_unzipped_stored\Post_eq\\2020-01-12_wv1"
    path = "G:\Ridgecrest\Pre_dem_2019"
    path = "/media/cosi/32103DCE103D9A35/2-Data/4-Shisper/0-Shisper_Data/WV_GE/0-Data"
    path = "/media/storage/Saif/3D_approch/Bulk Order NAIP_images_over_Ridgecrest/NAIP_images"
    path = "/media/cosicorr/storage/Saif2/NewWV_data"
    path = "//media/cosicorr/storage/Saif/1-Ridgecrest/6-WV/DG/P1BS_sorted/Before_eq/No_metadata/2016-09-08-WV03"
    path = "/media/cosicorr/storage/Saif/1-Ridgecrest/3D_approch_Ridgecrest/DEMS_correlation/FigureForPaper"
    path = "/media/cosicorr/storage/Saif/1-Ridgecrest/3D_approch_Ridgecrest/0-DEMs/WV_DEM"
    path = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/6p4_3D/Correlation/Set1"
    path = "/home/cosicorr/0-WorkSpace/CorrPaperTS_PS/fREQ_CORR/Freq_corr_B1_62_32_8_4"
    # FootprintExtract_Batch(inputdirectory=path, fileExtension=["tif"],
    #                        fpSavingDirectory=None)
    # orthImg = "G:\Ridgecrest\\2-Ridgecrest_DEMs\Ridgecrest_DEM_Lidar_USGS\Lidar_DEM_50cm.tif"
    # orthPlanet = "D:\\temp\\2019-07-04_181306_Blue_SR.tif"
    # RasterFootprint(rasterPath=orthImg)

    # path = "H:\\2-Data\\4-Shisper\\0-Shisper_Data\Shisper_Region_Planet_data\Sorted-Data\DOVE-C-2019-Shisper\\0-All-2019_PlanetData\\2019-08\\temp"
    # fp1 = os.path.join(path,"20190801_041726_0f49_1B_Analytic_fp.geojson")
    # fp2 = os.path.join(path,"20190801_041727_0f49_1B_Analytic_fp.geojson")
    # fp3 = os.path.join(path,"20190807_052421_1006_1B_Analytic_fp.geojson")
    # overlapObj = Overlapping(imgfp1=fp1,imgfp2=fp2,display=False)
    # rasterInof =RT.GetRasterInfo(path,True)
