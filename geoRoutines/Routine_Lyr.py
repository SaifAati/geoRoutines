import geojson
import geopandas
import gdal, ogr, osr
import numpy as np
import matplotlib.pyplot as plt

import geoRoutines.RandomColors as randomcolsRT


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


def RasterizeVector(shpInput, refRaster, output, noDataValue=-999999, dtype=gdal.GDT_Float32):
    import geospatialroutine.georoutines as geoRT
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

    # import shapefile
    # sf = shapefile.Reader(shpPath)
    # print(sf.shapeType)
    # shapes = sf.shapes()
    # for shape in shapes:
    #     for vertex in shape.points:
    #         print(vertex)
    return shpInfo


def WriteLayer(lyrOutputPath):
    # # Write point coordinates to Shapefile
    # shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    # if os.path.exists(lyrOutputPath):
    #     shpDriver.DeleteDataSource(lyrOutputPath)
    # outDataSource = shpDriver.CreateDataSource(lyrOutputPath)
    # outLayer = outDataSource.CreateLayer(lyrOutputPath, geom_type=ogr.wkbMultiPoint)
    # featureDefn = outLayer.GetLayerDefn()
    # outFeature = ogr.Feature(featureDefn)
    # outFeature.SetGeometry(multipoint)
    # outLayer.CreateFeature(outFeature)
    # outFeature = None
    return


def shp2array(shp_fn, r_ds=None, res=None, extent=None, t_srs=None, msg=False):
    """Rasterize input shapefile to match existing raster Dataset (or specified res/extent/t_srs)
    """
    if isinstance(shp_fn, ogr.DataSource):
        shp_ds = shp_fn
    else:
        shp_ds = ogr.Open(shp_fn)

    lyr = shp_ds.GetLayer()

    # This returns xmin, ymin, xmax, ymax
    shp_extent = lyr_extent(lyr)

    shp_srs = lyr.GetSpatialRef()
    # print(shp_srs)
    # dst_dt = gdal.GDT_Byte
    ndv = 0
    if r_ds is not None:
        r_extent = ds_extent(r_ds)
        res = get_res(r_ds, square=True)[0]
        if extent is None:
            extent = r_extent
        r_srs = get_ds_srs(r_ds)
        r_geom = ds_geom(r_ds)
        # dst_ns = r_ds.RasterXSize
        # dst_nl = r_ds.RasterYSize
        # Convert raster extent to shp_srs
        cT = osr.CoordinateTransformation(r_srs, shp_srs)
        r_geom_reproj = geom_dup(r_geom)
        r_geom_reproj.Transform(cT)
        r_geom_reproj.AssignSpatialReference(t_srs)
        lyr.SetSpatialFilter(r_geom_reproj)
        # lyr.SetSpatialFilter(ogr.CreateGeometryFromWkt(wkt))
    else:
        if msg:
            print(" -- No reference raster is provided ")
        if res is None:
            sys.exit("Must specify input res")
        if extent is None:
            if msg:
                print("Using input shp extent")
            extent = shp_extent
    if t_srs is None:
        t_srs = r_srs
        if msg:
            print("t_src:", t_srs)
    if not shp_srs.IsSame(t_srs):
        if msg:
            print("Input shp srs: %s" % shp_srs.ExportToProj4())
            print("Specified output srs: %s" % t_srs.ExportToProj4())
        out_ds = lyr_proj(lyr, t_srs)
        outlyr = out_ds.GetLayer()
    else:
        outlyr = lyr
    m_ds = mem_ds(res, extent, srs=t_srs, dtype=gdal.GDT_Byte)
    b = m_ds.GetRasterBand(1)
    geoTransfrom = m_ds.GetGeoTransform()
    b.SetNoDataValue(ndv)
    gdal.RasterizeLayer(m_ds, [1], outlyr, burn_values=[1])
    a = b.ReadAsArray()
    # a = a.astype('Bool')
    return a, geoTransfrom


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


def VectorDrivers(driver):
    if driver == "GeoJSON":
        return ".geojson"
    if driver == "KML":
        return ".kml"
    if driver == "ESRI Shapefile":
        return ".shp"


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
    colors = GenerateColors(N=len(featureList), pastel_factor=0.5)

    for index, temp in enumerate(featureList):
        temp.plot(ax=ax, facecolor='none', edgecolor=colors[index], linewidth=3.0)

        # ax.legend()
    # ax.legend(handles =labels)
    plt.show()

    return


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


##################################################################
## Reference: https://automating-gis-processes.github.io/CSC18/lessons/L2/geopandas-basics.html


def lyr_extent(lyr):
    # Envelope is ul_x, ur_x, lr_y, ll_y (?)
    env = lyr.GetExtent()
    # return xmin, ymin, xmax, ymax
    return [env[0], env[2], env[1], env[3]]


def lyr_proj(lyr, t_srs, preserve_fields=True):
    """Reproject an OGR layer
    """
    # Need to check t_srs
    s_srs = lyr.GetSpatialRef()
    cT = osr.CoordinateTransformation(s_srs, t_srs)

    # Do everything in memory
    drv = ogr.GetDriverByName('Memory')

    # Might want to save clipped, warped shp to disk?
    # create the output layer
    # drv = ogr.GetDriverByName('ESRI Shapefile')
    # out_fn = '/tmp/temp.shp'
    # if os.path.exists(out_fn):
    #    driver.DeleteDataSource(out_fn)
    # out_ds = driver.CreateDataSource(out_fn)

    out_ds = drv.CreateDataSource('out')
    outlyr = out_ds.CreateLayer('out', srs=t_srs, geom_type=lyr.GetGeomType())

    if preserve_fields:
        # add fields
        inLayerDefn = lyr.GetLayerDefn()
        for i in range(0, inLayerDefn.GetFieldCount()):
            fieldDefn = inLayerDefn.GetFieldDefn(i)
            outlyr.CreateField(fieldDefn)
        # get the output layer's feature definition
    outLayerDefn = outlyr.GetLayerDefn()

    # loop through the input features
    inFeature = lyr.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(cT)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        if preserve_fields:
            for i in range(0, outLayerDefn.GetFieldCount()):
                outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outlyr.CreateFeature(outFeature)
        # destroy the features and get the next input feature
        inFeature = lyr.GetNextFeature()
    # NOTE: have to operate on ds here rather than lyr, otherwise segfault
    return out_ds


def gt_corners(gt, nx, ny):
    """Get corner coordinates based on input geotransform and raster dimensions
   """
    ul = [gt[0], gt[3]]
    ll = [gt[0], gt[3] + (gt[5] * ny)]
    ur = [gt[0] + (gt[1] * nx), gt[3]]
    lr = [gt[0] + (gt[1] * nx), gt[3] + (gt[5] * ny)]
    return ul, ll, ur, lr


def get_ds_srs(ds):
    """Get srs object for GDAL Datset
    """
    ds_srs = osr.SpatialReference()
    ds_srs.ImportFromWkt(ds.GetProjectionRef())
    return ds_srs


def corner_extent(ul, ll, ur, lr):
    """Get min/max extent based on corner coord
    """
    xmin = min(ul[0], ll[0], ur[0], lr[0])
    xmax = max(ul[0], ll[0], ur[0], lr[0])
    ymin = min(ul[1], ll[1], ur[1], lr[1])
    ymax = max(ul[1], ll[1], ur[1], lr[1])
    extent = [xmin, ymin, xmax, ymax]
    return extent


def ds_extent(ds, t_srs=None):
    """Return min/max extent of dataset based on corner coordinates

    xmin, ymin, xmax, ymax

    If t_srs is specified, output will be converted to specified srs
    """
    ul, ll, ur, lr = gt_corners(ds.GetGeoTransform(), ds.RasterXSize, ds.RasterYSize)
    ds_srs = get_ds_srs(ds)
    if t_srs is not None and not ds_srs.IsSame(t_srs):
        ct = osr.CoordinateTransformation(ds_srs, t_srs)
        # Check to see if ct creation failed
        # if ct == NULL:
        # Check to see if transform failed
        # if not ct.TransformPoint(extent[0], extent[1]):
        # Need to check that transformed coordinates fall within appropriate bounds
        ul = ct.TransformPoint(*ul)
        ll = ct.TransformPoint(*ll)
        ur = ct.TransformPoint(*ur)
        lr = ct.TransformPoint(*lr)
    extent = corner_extent(ul, ll, ur, lr)
    return extent


def get_res(ds, t_srs=None, square=False):
    """Get GDAL Dataset raster resolution
    """
    gt = ds.GetGeoTransform()
    ds_srs = get_ds_srs(ds)
    # This is Xres, Yres
    res = [gt[1], np.abs(gt[5])]
    if square:
        res = [np.mean(res), np.mean(res)]
    if t_srs is not None and not ds_srs.IsSame(t_srs):
        if True:
            # This diagonal approach is similar to the approach in gdaltransformer.cpp
            # Bad news for large extents near the poles
            # ullr = get_ullr(ds, t_srs)
            # diag = np.sqrt((ullr[0]-ullr[2])**2 + (ullr[1]-ullr[3])**2)
            extent = ds_extent(ds, t_srs)
            diag = np.sqrt((extent[2] - extent[0]) ** 2 + (extent[3] - extent[1]) ** 2)
            res = diag / np.sqrt(ds.RasterXSize ** 2 + ds.RasterYSize ** 2)
            res = [res, res]
        else:
            # Compute from center pixel
            ct = osr.CoordinateTransformation(ds_srs, t_srs)
            pt = get_center(ds)
            # Transform center coordinates
            pt_ct = ct.TransformPoint(*pt)
            # Transform center + single pixel offset coordinates
            pt_ct_plus = ct.TransformPoint(pt[0] + gt[1], pt[1] + gt[5])
            # Compute resolution in new units
            res = [pt_ct_plus[0] - pt_ct[0], np.abs(pt_ct_plus[1] - pt_ct[1])]
    return res


def applyGeoTransform(inX, inY, geoTransform):
    inX = np.asarray(inX)
    inY = np.asarray(inY)
    outX = geoTransform[0] + inX * geoTransform[1] + inY * geoTransform[2]
    outY = geoTransform[3] + inX * geoTransform[4] + inY * geoTransform[5]
    return outX, outY


# Add 0.5 px offset to account for GDAL model - gt 0,0 is UL corner, pixel 0,0 is center
def pixelToMap(pX, pY, geoTransform):
    """Convert pixel coordinates to map coordinates based on geotransform

    Accepts float or NumPy arrays

    GDAL model used here - upper left corner of upper left pixel for mX, mY (and in GeoTransform)
    """
    pX = np.asarray(pX, dtype=float)
    pY = np.asarray(pY, dtype=float)
    pX += 0.5
    pY += 0.5
    mX, mY = applyGeoTransform(pX, pY, geoTransform)
    return mX, mY


def geom_transform(geom, t_srs):
    """Transform a geometry in place
    """
    s_srs = geom.GetSpatialReference()
    if not s_srs.IsSame(t_srs):
        ct = osr.CoordinateTransformation(s_srs, t_srs)
        geom.Transform(ct)
        geom.AssignSpatialReference(t_srs)


def ds_geom(ds, t_srs=None):
    """Return dataset bbox envelope as geom
    """
    gt = ds.GetGeoTransform()
    ds_srs = get_ds_srs(ds)
    if t_srs is None:
        t_srs = ds_srs
    ns = ds.RasterXSize
    nl = ds.RasterYSize
    x = np.array([0, ns, ns, 0, 0], dtype=float)
    y = np.array([0, 0, nl, nl, 0], dtype=float)
    # Note: pixelToMap adds 0.5 to input coords, need to account for this here
    x -= 0.5
    y -= 0.5
    mx, my = pixelToMap(x, y, gt)
    geom_wkt = 'POLYGON(({0}))'.format(', '.join(['{0} {1}'.format(*a) for a in zip(mx, my)]))
    geom = ogr.CreateGeometryFromWkt(geom_wkt)
    geom.AssignSpatialReference(ds_srs)
    if not ds_srs.IsSame(t_srs):
        geom_transform(geom, t_srs)
    return geom


def geom_dup(geom):
    """Create duplicate geometry

    Needed to avoid segfault when passing geom around. See: http://trac.osgeo.org/gdal/wiki/PythonGotchas
    """
    g = ogr.CreateGeometryFromWkt(geom.ExportToWkt())
    g.AssignSpatialReference(geom.GetSpatialReference())
    return g


def mem_ds(res, extent, srs=None, dtype=gdal.GDT_Float32):
    """Create a new GDAL Dataset in memory

    Useful for various applications that require a Dataset
    """
    # These round down to int
    # dst_ns = int((extent[2] - extent[0])/res)
    # dst_nl = int((extent[3] - extent[1])/res)
    # This should pad by 1 pixel, but not if extent and res were calculated together to give whole int
    dst_ns = int((extent[2] - extent[0]) / res + 0.99)
    dst_nl = int((extent[3] - extent[1]) / res + 0.99)
    m_ds = gdal.GetDriverByName('MEM').Create('', dst_ns, dst_nl, 1, dtype)
    m_gt = [extent[0], res, 0, extent[3], 0, -res]
    m_ds.SetGeoTransform(m_gt)
    if srs is not None:
        m_ds.SetProjection(srs.ExportToWkt())
    return m_ds


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
