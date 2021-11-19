import geojson, geopandas, os
import gdal, ogr, osr
import numpy as np
import matplotlib.pyplot as plt

import geoRoutines.RandomColors as randomcolsRT
from geoRoutines.georoutines import ComputeEpsg,ConvCoordMap1ToMap2


######################################## Read Vectors ##################################################################
def GetVectorInfo(shpInput, displayInfo=True):
    vectorData = ReadVector(shpInput=shpInput)
    dataInfo = {"Type": type(vectorData), "Head": vectorData.head, "CRS": vectorData.crs,
                "EPSG": vectorData.crs.to_epsg()}

    if isinstance(shpInput, ogr.DataSource):
        shp_ds = shpInput
    else:
        shp_ds = ogr.Open(shpInput)

    lyr = shp_ds.GetLayer()
    dataInfo["lyr"] = lyr
    s_srs = lyr.GetSpatialRef()
    dataInfo["srs"] = s_srs
    if dataInfo:
        print(dataInfo)
        # print("Degug ::",vectorData["geometry"].head())
        # print(len(vectorData))
    return dataInfo


def ReadVector(shpInput, show=False):
    # Read file using gpd.read_file()
    data = geopandas.read_file(shpInput)
    if show:
        data.plot()
        plt.show()
    return data


def ReadGeojson(geojsonPath):
    with open(geojsonPath) as f:
        data = geojson.load(f)
    return data


def ReadVector_Batch(shpList):
    fpDataFrame = ReadVector(shpList[0])
    fpDF_List = []
    fpDF_List.append(fpDataFrame)
    crs = fpDF_List[0].crs

    for i in range(1, len(shpList)):
        fpDataFrame = ReadVector(shpList[i])
        fpDF_List.append(fpDataFrame)
        if crs != fpDF_List[i].crs:
            # print("Layers does not have the same CRS")
            fpDF_List[i] = fpDF_List[i].to_crs({'init': fpDF_List[0].crs})

    return fpDF_List


def ReadShapeFile(shpPath):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(shpPath)

    # from Layer
    layer = shp.GetLayer()
    spatialRef = layer.GetSpatialRef()
    # print("spatialRef",spatialRef)
    # from Geometry
    feature = layer.GetNextFeature()
    geom = feature.GetGeometryRef()
    # ring = geom.GetGeometryRef(0)
    point = geom.GetPointCount()
    spatialRef = geom.GetSpatialReference()
    EPSG_Code = spatialRef.GetAttrValue('AUTHORITY', 1)
    shpInfo = {"Shp": shp, "ShpPath": shpPath, "ShpName": os.path.split(shpPath)[-1],
               "MapInfo": [EPSG_Code, spatialRef], "Geom": geom}

    geom = shpInfo.get("Geom")
    points = geom.GetPointCount()
    vertex = []
    for point in range(points):
        vertex.append(geom.GetPoint(point))

    shpInfo["Vertex"] = vertex

    return shpInfo


################################ Write vectors ########################################################################
def WriteJson(features, outputFile, driver="GeoJSON"):
    """

    :param features:
    :param outputFile:
    :param driver:
    :return:
    """

    outputFile = outputFile + VectorDrivers(driver)

    with open(outputFile, 'w') as outfile:
        geojson.dump(features, outfile, sort_keys=True)
    return outputFile


############################### Tools #################################################################################
class Intersection:

    def __init__(self, fp1, fp2, dispaly=False):
        self.fpDF_List = []
        self.intersection = 0
        self.fp1_area = 0
        self.fp2_area = 0
        self.res_inter = None
        self.overlapPerc_fp1 = 0
        self.overlapPerc_fp2 = 0
        self.inter_area = 0
        self.fp1 = fp1
        self.fp2 = fp2
        self.fpList = [self.fp1, self.fp2]
        # print(self.fpDF_List)
        self.__ReadShp()
        # print(self.fpDF_List)
        self.Intersect()
        self.intersection = self.intersection
        self.res_inter = self.res_inter

        if dispaly:
            self.Visualise_overlapping()

    def __ReadShp(self):
        for fp_ in self.fpList:
            fpDataFrame = ReadVector(fp_)
            self.fpDF_List.append(fpDataFrame)

        if self.fpDF_List[0].crs != self.fpDF_List[1].crs:
            print("Layers does not have the same CRS")
            print("*************", self.fpDF_List[0].crs, "*************", self.fpDF_List[1].crs, "**************")
            self.fpDF_List[1] = self.fpDF_List[1].to_crs({'init': self.fpDF_List[0].crs})

    def Intersect(self):
        import geoRoutines.georoutines as geoRT
        self.res_inter = geopandas.overlay(self.fpDF_List[0], self.fpDF_List[1], how='intersection')
        # print(self.res_inter)
        if self.res_inter.index.start == self.res_inter.index.stop:
            self.intersection = 0
        else:
            self.intersection = geojson.Feature(geometry=self.res_inter.loc[0, "geometry"])

        if self.intersection != 0:
            dfCopy = self.res_inter.copy()
            coords = self.intersection['geometry']['coordinates']

            """
            FutureWarning: '+init=<authority>:<code>' syntax is deprecated. '<authority>:<code>' is the preferred initialization method. When making the change, be mindful of axis order changes: https://pyproj4.github.io/pyproj/stable/gotchas.html#axis-order-changes-in-proj-6
            return _prepare_from_string(" ".join(pjargs))
            epsg = "epsg:" + str(RT.ComputeEpsg(lon=coords[0][0][0], lat=coords[0][0][1]))
            dfCopy = dfCopy.to_crs({'init': epsg}) 
            ==> dfCopy = dfCopy.to_crs(RT.ComputeEpsg(lon=coords[0][0][0], lat=coords[0][0][1])) 
            """

            dfCopy = dfCopy.to_crs(epsg=geoRT.ComputeEpsg(lon=coords[0][0][0], lat=coords[0][0][1]))
            self.inter_area = (dfCopy['geometry'].area / 10 ** 6)[0]
            # print(self.inter_area)

            fp1_copy = self.fpDF_List[0].copy().to_crs(epsg=geoRT.ComputeEpsg(lon=coords[0][0][0], lat=coords[0][0][1]))

            self.fp1_area = (fp1_copy["geometry"].area / 10 ** 6)[0]
            # print(self.fp1_area)
            fp2_copy = self.fpDF_List[1].copy().to_crs(epsg=geoRT.ComputeEpsg(lon=coords[0][0][0], lat=coords[0][0][1]))
            self.fp2_area = (fp2_copy["geometry"].area / 10 ** 6)[0]
            # print(self.fp2_area)
            self.overlapPerc_fp1 = (self.inter_area / self.fp1_area) * 100
            self.overlapPerc_fp2 = (self.inter_area / self.fp2_area) * 100

            # print("Intersection area (Km^2): ", self.inter_area, " Overlapping with respect to",
            #       os.path.basename(self.fp1), " (%): ", '%1.1f' % (self.overlapPerc_fp1))
            # print("Intersection area (Km^2): ", self.inter_area, " Overlapping with respect to",
            #       os.path.basename(self.fp1), " (%): ", '%1.1f' % (self.overlapPerc_fp2))
            # else:
            #     self.overlapPerc = 0
            # return self.overlapPerc

    def Visualise_overlapping(self):
        fig, ax = plt.subplots()
        colors = randomcolsRT.GenerateColors(N=len(self.fpDF_List), pastel_factor=0.5)
        if self.intersection != 0:
            self.res_inter.plot(ax=ax, alpha=0.5, cmap='tab10')
        self.fpDF_List[0].plot(ax=ax, facecolor='none', edgecolor=colors[0])
        self.fpDF_List[1].plot(ax=ax, facecolor='none', edgecolor=colors[1])
        plt.show()


def VectorDrivers(driver):
    if driver == "GeoJSON":
        return ".geojson"
    if driver == "KML":
        return ".kml"
    if driver == "ESRI Shapefile":
        return ".shp"


def RasterizeVector(shpInput, refRaster, output, noDataValue=-999999, dtype=gdal.GDT_Float32):
    import geoRoutines.georoutines as geoRT
    rasterInfo = geoRT.RasterInfo(refRaster)
    raster = rasterInfo.raster

    originX = rasterInfo.xOrigin
    originY = rasterInfo.yOrigin
    pixelWidth = rasterInfo.pixelWidth
    pixelHeight = rasterInfo.pixelHeight
    rows, cols = np.shape(rasterInfo.ImageAsArray())
    mb_v = ogr.Open(shpInput)
    mb_l = mb_v.GetLayer()

    ## Prepare destination file
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(output, cols, rows, 1, dtype)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

    ## Set the projection
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outRaster.SetMetadataItem("Author", "SAIF AATI saif@caltech.edu")

    # band.SetNoDataValue(noDataValue)
    # band.FlushCache()
    gdal.RasterizeLayer(outRaster, [1], mb_l)
    outRaster = None

    return output

def ShpDistance(shpDataFrame):
    """

    :param shpDataFrame: the return of  ReadVector
    :return:
    """
    data = shpDataFrame
    if data.empty:
        return 0
    else:
        # data.crs.from_string()
        datacrs = str(data.crs)
        # print(datacrs)
        # print(datacrs.split("epsg:")[1][0:4])
        if datacrs.split("epsg:")[1][0:4] == "4326":
            # print("convertCoordinates")
            lon = data["geometry"][0].centroid.xy[0][0]
            lat = data["geometry"][0].centroid.xy[1][0]
            utmEpsg = "epsg:" + str(ComputeEpsg(lon, lat))
            # print(utmEpsg)
            data = data.to_crs({'init': utmEpsg})
        p1 = data["geometry"][0].centroid
        dcum = 0
        for item_ in data["geometry"]:
            p2 = item_.centroid
            d = p1.distance(p2)
            dcum += d
            p1 = item_.centroid
        # print(dcum)
        return dcum

def Intersection_(fp1, fp2, display=False):
    # fp1Info = GetVectorInfo(fp1)
    # fp2Info = GetVectorInfo(fp2)
    # if fp1Info["EPSG"] != fp2Info["EPSG"]:
    #     print("Layers does not have the same CRS")
    "https://www.earthdatascience.org/workshops/gis-open-source-python/reproject-vector-data-in-python/"
    fpList = [fp1, fp2]
    print(fp1, "\n", fp2)
    fpDF_List = []
    for fp_ in fpList:
        fpDataFrame = geopandas.read_file(fp_)
        # print(fpDataFrame.crs)
        fpDF_List.append(fpDataFrame)
    if fpDF_List[0].crs != fpDF_List[1].crs:
        print("Layers does not have the same CRS")
        fpDF_List[1] = fpDF_List[1].to_crs({'init': fpDF_List[0].crs})

    # intersection = 0
    res_inter = geopandas.overlay(fpDF_List[0], fpDF_List[1], how='intersection')
    # print(res_inter)
    if res_inter.index.start == res_inter.index.stop:
        intersection = 0
    else:
        intersection = geojson.Feature(geometry=res_inter.loc[0, "geometry"])
    if display:

        fig, ax = plt.subplots()
        colors = GenerateColors(N=len(fpDF_List), pastel_factor=0.5)
        if intersection != 0:
            res_inter.plot(ax=ax, alpha=0.5, cmap='tab10')
        fpDF_List[0].plot(ax=ax, facecolor='none', edgecolor=colors[0])
        fpDF_List[1].plot(ax=ax, facecolor='none', edgecolor=colors[1])
        plt.show()

    if intersection != 0:
        dfCopy = res_inter.copy()
        coords = intersection['geometry']['coordinates']
        print(coords)
        epsg = "epsg:" + str(RT.ComputeEpsg(lon=coords[0][0][0], lat=coords[0][0][1]))
        print(epsg)
        dfCopy = dfCopy.to_crs({'init': epsg})
        inter_area = (dfCopy['geometry'].area / 10 ** 6)[0]
        # print(inter_area)
        fp1_copy = fpDF_List[0].copy().to_crs({'init': epsg})
        fp1_area = (fp1_copy["geometry"].area / 10 ** 6)[0]
        # print(fp1_area)
        fp2_copy = fpDF_List[1].copy().to_crs({'init': epsg})
        fp2_area = (fp2_copy["geometry"].area / 10 ** 6)[0]

        # self.overlapPerc = '%1.1f' % ((self.inter_area / fpArea) * 100)
        # print("Intersection area (Km^2): ", self.inter_area, " Overlapping %: ", self.overlapPerc)
    # else:
    #     self.overlapPerc = 0
    # return self.overlapPerc

    return intersection


def GetBoundsGeojson(geojsonFile, conv2UTM=False):
    if isinstance(geojsonFile, str):
        import json
        with open(geojsonFile) as json_file:
            json_data = json.load(json_file)
            poly = json_data['geometry']
            coords = poly['coordinates']

    if isinstance(geojsonFile, dict):
        geometry = geojsonFile
        coords = geometry['geometry']['coordinates']

    lonList = [i for i, j in coords[0]]
    latList = [j for i, j in coords[0]]

    if conv2UTM:
        # print("--- Convert Geojson from WGS84 to UTM")
        easting = []
        northing = []
        for index, lon in enumerate(lonList):
            epsg = ComputeEpsg(lon=lon, lat=latList[index])

            coord = ConvCoordMap1ToMap2(x=latList[index], y=lon, targetEPSG=epsg)

            northing.append(coord[1])
            easting.append(coord[0])

        return easting, northing, epsg
    else:
        epsg = 4326
        return lonList, latList, epsg




######################################## Plotting ######################################################################
def PlotShps(shpList, title=None):
    fpDF_List = ReadVector_Batch(shpList)
    fig, ax = plt.subplots(figsize=(10, 15))
    colors = randomcolsRT.GenerateColors(N=len(fpDF_List), pastel_factor=0.5)

    for index, temp in enumerate(fpDF_List):
        temp.plot(ax=ax, facecolor='none', edgecolor=colors[index], linewidth=3.0)
        if title:
            ax.set_title(title)
        # ax.legend()
    # ax.legend(handles =labels)
    # tempPath = "/home/cosicorr/Desktop/0-WorkSpace/3D-Correlation_project"
    # plt.savefig(os.path.join(tempPath,title+".png"))
    plt.show()


def PlotShps_features(featureList):
    # fpDF_List = ReadVector_Batch(shpList)
    fig, ax = plt.subplots()
    colors = randomcolsRT.GenerateColors(N=len(featureList), pastel_factor=0.5)

    for index, temp in enumerate(featureList):
        temp.plot(ax=ax, facecolor='none', edgecolor=colors[index], linewidth=3.0)

        # ax.legend()
    # ax.legend(handles =labels)
    plt.show()

    return


if __name__ == '__main__':
    # path = "D:\TempCode"
    # fp1 = "D:\TempCode\Main_Rupture_7.1\\7p1_faultTrace_LineString.geojson"
    # fp2 = os.path.join(path, "16JUL31183128-P1BS-504089615020_01_P003_fp.geojson")
    # fp3 = os.path.join(path, "16SEP08185521-P1BS-504089615070_01_P008_fp.geojson")
    # inter = Intersection(fp1=fp2, fp2=fp1, display=False)
    # print(inter == 0)

    shpPath = "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/ShisperMask.geojson"
    oshPath = "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/ShisperMask_rasterized.vrt"
    RasterizeVector(shpInput=shpPath,
                    refRaster="/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/0-AsterDEM_UTM_resampled.vrt",
                    output=oshPath)
