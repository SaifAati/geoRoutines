import os
import sys
import warnings

import numpy as np

import gdal, gdalconst
import ogr, osr
import rasterio
import matplotlib.pyplot as plt
from pathlib import Path

from scipy.stats import norm
from scipy import stats
# from numba import jit

def GetRasterInfo(inputRaster, printInfo=False):
    """
    :Method: Retrieve all information from Raster image
    :param imagePath: the path of the image : string
    :return: dictionary of metadata
    """
    error = False
    try:
        ## Open the image using GDAL
        raster = gdal.Open(inputRaster)
        geoTrans = raster.GetGeoTransform()
    except:
        print('Problem when loading image with GDAL:', inputRaster, sys.exc_info())
        error = True
        return {"Error": error}
    if error == False:
        xOrigin = geoTrans[0]
        yOrigin = geoTrans[3]
        pixelWidth = geoTrans[1]
        pixelHeight = geoTrans[5]
        nbBand = raster.RasterCount
        rasterWidth = raster.RasterXSize
        rasterHeight = raster.RasterYSize
        rasterDataType = []
        for i in range(nbBand):
            band = raster.GetRasterBand(1)
            rasterDataType.append(gdal.GetDataTypeName(band.DataType))
            del band
        try:
            validMapInfo = True
            projection = raster.GetProjection()
            proj = osr.SpatialReference(wkt=raster.GetProjection())
            EPSG_Code = proj.GetAttrValue('AUTHORITY', 1)
            src = rasterio.open(inputRaster)
            crs = src.crs
            epsg_spare = int(str(crs).split(":")[1])
            if epsg_spare != EPSG_Code:
                EPSG_Code = epsg_spare

            dataInfo = {"Error": error, "Raster": raster, "MapInfo": [validMapInfo, geoTrans, EPSG_Code, projection],
                        "NbBands": nbBand, "DataType": rasterDataType,
                        "OriginCoord": (xOrigin, yOrigin), "Resolution": (pixelWidth, pixelHeight),
                        "Dimension": (rasterWidth, rasterHeight), "ImagePath": inputRaster,
                        "ImageName": os.path.split(inputRaster)[-1]}

            bandInfo = []
            for band_ in range(nbBand):
                bandInfo.append(raster.GetRasterBand(band_ + 1).GetDescription())
            dataInfo["BandInfo"] = bandInfo

            dataInfo["MetaData"] = raster.GetMetadata()

            rpcs = raster.GetMetadata('RPC')
            dataInfo["RPC"] = rpcs

            if EPSG_Code == None and np.allclose(GeoTransfomAsArray(geoTrans), np.eye(3)):
                warnings.warn(
                    "NotGeoreferencedWarning: raster has no geotransform set. The identity matrix will be returned.")

        except:
            validMapInfo = False
            dataInfo = {"Error": error, "Raster": raster, "MapInfo": [validMapInfo], "NbBands": nbBand,
                        "dataType": rasterDataType,
                        "OriginCoord": (xOrigin, yOrigin), "Resolution": (pixelWidth, pixelHeight),
                        "Dimension": (rasterWidth, rasterHeight), "ImagePath": inputRaster,
                        "ImageName": os.path.split(inputRaster)[-1]}
            dataInfo["MetaData"] = raster.GetMetadata()
            rpcs = raster.GetMetadata('RPC')
            dataInfo["RPC"] = rpcs

            bandInfo = []
            for band_ in range(nbBand):
                bandInfo.append(raster.GetRasterBand(band_ + 1).GetDescription())
            dataInfo["BandInfo"] = bandInfo
        if printInfo:
            print(dataInfo)
        return dataInfo


def GeoTransfomAsArray(geoTrans):
    xOrigin = geoTrans[0]
    yOrigin = geoTrans[3]
    pixelWidth = geoTrans[1]
    pixelHeight = geoTrans[5]
    rtnX = geoTrans[2]  #
    rtnY = geoTrans[4]  #
    ## Apply affine transformation
    trans = np.zeros((3, 3))
    trans[0, 0] = pixelWidth
    trans[0, 1] = rtnX
    trans[1, 0] = rtnY
    trans[1, 1] = pixelHeight
    trans[2, 0] = xOrigin
    trans[2, 1] = yOrigin
    trans[2, 2] = 1

    return trans


def ImageAsArray(imageInfo, bandNumber=1):
    """
    :param imagePath:
    :param bandNumber:
    :return: image as array
    """
    raster = imageInfo.get("Raster")
    # Transform image to array
    imageAsArray = np.array(raster.GetRasterBand(bandNumber).ReadAsArray())

    return imageAsArray


def MultiBandsRaster2Array(imageInfo):
    """
    Improve with Xarray
    :param imageInfo:
    :return:
    """

    array = np.empty((imageInfo["NbBands"], imageInfo["Dimension"][1], imageInfo["Dimension"][0]))
    for i in range(imageInfo["NbBands"]):
        array[i] = ImageAsArray(imageInfo=imageInfo, bandNumber=i + 1)

    return array


def Pixel2Map(imageInfo, x, y):
    """
    :Method : convert pixel coordinate to map coordinate,
    :Note: The top Left coordinate of the image with GDAl correspond to (0,0)pix
    :param imagePath: path of the image  : string
    :param x: xPixel coordinate : int or float
    :param y: yPixel coordinate: int or float
    :return: xMap,yMap : tuple  (non integer coordinates)
    """
    xOrigin = imageInfo["OriginCoord"][0]  # x top left coordinate (0.0 pix)
    yOrigin = imageInfo["OriginCoord"][1]  # y top left coordinate (0.0 pix)
    pixelWidth = imageInfo["Resolution"][0]  # xGSD
    pixelHeight = imageInfo["Resolution"][1]  # yGSD
    rtnX = imageInfo["MapInfo"][1][2]  #
    rtnY = imageInfo["MapInfo"][1][4]  #
    ## Apply affine transformation
    mat = np.array([[pixelWidth, rtnX], [rtnY, pixelHeight]])
    trans = np.array([[xOrigin, yOrigin]])
    res = np.dot(mat, np.array([[x, y]]).T) + trans.T
    xMap = res[0].item()
    yMap = res[1].item()
    return (xMap, yMap)


# @jit (nopython=True)
def Pixel2Map_Batch(imageInfo, X, Y):
    """
    :Method : convert pixel coordinate to map coordinate,
    :Note: The top Left coordinate of the image with GDAl correspond to (0,0)pix
    :param imagePath: path of the image  : string
    :param X: list of xPixel coordinates : list of int of float
    :param Y: list of  yPixel coordinates:  list of int or float
    :return: xMap,yMap : tuple  (non integer coordinates)
    """

    # xOrigin = imageInfo["OriginCoord"][0]  # x top left coordinate (0.0 pix)
    # yOrigin = imageInfo["OriginCoord"][1]  # y top left coordinate (0.0 pix)
    # pixelWidth = imageInfo["Resolution"][0]  # xGSD
    # pixelHeight = imageInfo["Resolution"][1]  # yGSD
    # rtnX = imageInfo["MapInfo"][1][2]  #
    # rtnY = imageInfo["MapInfo"][1][4]  #
    # X = np.asarray(X)
    # Y = np.asarray(Y)
    # XY = np.hstack((X, Y))
    # XY_sparse = sparse.csr_matrix(XY).T
    # if X.shape[0] != Y.shape[0]:
    #     print("Erroor: x and Y must have the same size !!!")
    #     return
    # mat1 = sparse.kron(sparse.eye(X.shape[0]), [pixelWidth, rtnX])
    # mat2 = sparse.kron(sparse.eye(Y.shape[0]), [rtnY, pixelHeight])
    # mat = sparse.vstack(blocks=[mat1, mat2])
    # t1 = np.asarray([xOrigin] * X.shape[0])
    # t2 = np.asarray([yOrigin] * X.shape[0])
    # Tran = np.hstack((t1, t2))
    # Tran_sparse = sparse.csr_matrix(Tran).T
    # res_sparse = sparse.csr_matrix.dot(mat, XY_sparse) + Tran_sparse
    # res = sparse.csr_matrix.todense(res_sparse)
    # # print(res.tolist())
    # X_map = res[0:X.shape[0]].tolist()
    # X_map = np.concatenate(X_map, axis=0)
    # Y_map = res[X.shape[0]:X.shape[0] + Y.shape[0]].tolist()
    # # print(Y_map)
    # Y_map = np.concatenate(Y_map, axis=0)
    # print(Y_map,Y_map.shape)
    nb = len(X)
    X_map = []
    Y_map = []
    for x, y in zip(X, Y):
        xMap_, yMap_ = Pixel2Map(imageInfo=imageInfo, x=x, y=y)
        X_map.append(xMap_)
        Y_map.append(yMap_)
    return (X_map, Y_map)


def Pixel2Map_Batch_Parallel(imageInfo, X, Y):
    """
    :Method : convert pixel coordinate to map coordinate,
    :Note: The top Left coordinate of the image with GDAl correspond to (0,0)pix
    :param imagePath: path of the image  : string
    :param X: list of xPixel coordinates : list of int of float
    :param Y: list of  yPixel coordinates:  list of int or float
    :return: xMap,yMap : tuple  (non integer coordinates)
    """

    # xOrigin = imageInfo["OriginCoord"][0]  # x top left coordinate (0.0 pix)
    # yOrigin = imageInfo["OriginCoord"][1]  # y top left coordinate (0.0 pix)
    # pixelWidth = imageInfo["Resolution"][0]  # xGSD
    # pixelHeight = imageInfo["Resolution"][1]  # yGSD
    # rtnX = imageInfo["MapInfo"][1][2]  #
    # rtnY = imageInfo["MapInfo"][1][4]  #
    # X = np.asarray(X)
    # Y = np.asarray(Y)
    # XY = np.hstack((X, Y))
    # XY_sparse = sparse.csr_matrix(XY).T
    # if X.shape[0] != Y.shape[0]:
    #     print("Erroor: x and Y must have the same size !!!")
    #     return
    # mat1 = sparse.kron(sparse.eye(X.shape[0]), [pixelWidth, rtnX])
    # mat2 = sparse.kron(sparse.eye(Y.shape[0]), [rtnY, pixelHeight])
    # mat = sparse.vstack(blocks=[mat1, mat2])
    # t1 = np.asarray([xOrigin] * X.shape[0])
    # t2 = np.asarray([yOrigin] * X.shape[0])
    # Tran = np.hstack((t1, t2))
    # Tran_sparse = sparse.csr_matrix(Tran).T
    # res_sparse = sparse.csr_matrix.dot(mat, XY_sparse) + Tran_sparse
    # res = sparse.csr_matrix.todense(res_sparse)
    # # print(res.tolist())
    # X_map = res[0:X.shape[0]].tolist()
    # X_map = np.concatenate(X_map, axis=0)
    # Y_map = res[X.shape[0]:X.shape[0] + Y.shape[0]].tolist()
    # # print(Y_map)
    # Y_map = np.concatenate(Y_map, axis=0)
    # print(Y_map,Y_map.shape)
    nb = len(X)
    X_map = []
    Y_map = []
    # for i in range(nb):
    #     xMap_, yMap_ = Pixel2Map(imageInfo=imageInfo, x=X[i], y=Y[i])
    #     X_map.append(xMap_)
    #     Y_map.append(yMap_)
    # print(imageInfo)
    # print(X)
    # print(Y)
    imageInfo_ = {"OriginCoord": imageInfo["OriginCoord"], "Resolution": imageInfo["Resolution"],
                  "MapInfo": imageInfo["MapInfo"]}
    from joblib import Parallel, delayed

    def my_fun(i):
        """ We define a simple function here.
        """
        # time.sleep(1)
        xMap_, yMap_ = Pixel2Map(imageInfo=imageInfo_, x=X[i], y=Y[i])
        return xMap_, yMap_

    # ImageInfo = [imageInfo]*len(X)
    res = Parallel(n_jobs=20)(delayed(my_fun)(i) for i in range(nb))
    # print(len(res))
    resArray = np.asarray(res)
    X_Map = list(resArray[:, 0])
    Y_map = list(resArray[:, 1])
    return X_Map, Y_map


def Map2Pixel(imageInfo, x, y):
    """
    :Method : convert coordinate from map space to image space,
        Note: The top Left coordinate of the image with GDAl correspond to (0,0)pix
    :param imagePath: path of the image  : string
    :param x: xMap coordinate : int or float
    :param y: yMap coordinate: int or float
    :return: coordinate in image space : tuple in pix
    """
    xOrigin = imageInfo["OriginCoord"][0]  # x top left coordinate (0.0 pix)
    yOrigin = imageInfo["OriginCoord"][1]  # y top left coordinate (0.0 pix)
    pixelWidth = imageInfo["Resolution"][0]  # xGSD
    pixelHeight = imageInfo["Resolution"][1]  # yGSD
    rtnX = imageInfo["MapInfo"][1][2]
    rtnY = imageInfo["MapInfo"][1][4]

    ## Apply inverse affine transformation
    mat = np.array([[pixelWidth, rtnX], [rtnY, pixelHeight]])
    trans = np.array([[xOrigin, yOrigin]])
    temp = np.array([[x, y]]).T - trans.T
    res = np.dot(np.linalg.inv(mat), temp)
    xPx = res[0].item()
    yPx = res[1].item()

    if xPx < 0 or xPx > imageInfo.get("Dimension")[0]:
        warnings.warn("xPix outside the image dimension", DeprecationWarning, stacklevel=2)
    if yPx < 0 or yPx > imageInfo.get("Dimension")[1]:
        warnings.warn("yPix outside the image dimension", DeprecationWarning, stacklevel=2)

    return (xPx, yPx)


def Map2Pixel_Batch(imageInfo, X, Y):
    """
    :Method : convert coordinate from map space to image space,
        Note: The top Left coordinate of the image with GDAl correspond to (0,0)pix
    :param imagePath: path of the image  : string
    :param X: list of xMap coordinate : list of int or float
    :param Y: list of yMap coordinate: list of int or float
    :return: coordinate in image space : tuple in pix
    """
    # from scipy import sparse
    # from scipy.sparse.linalg import pinv
    xOrigin = imageInfo["OriginCoord"][0]  # x top left coordinate (0.0 pix)
    yOrigin = imageInfo["OriginCoord"][1]  # y top left coordinate (0.0 pix)
    pixelWidth = imageInfo["Resolution"][0]  # xGSD
    pixelHeight = imageInfo["Resolution"][1]  # yGSD
    rtnX = imageInfo["MapInfo"][1][2]
    rtnY = imageInfo["MapInfo"][1][4]

    # ## Apply inverse affine transformation
    # mat = np.array([[pixelWidth, rtnX], [rtnY, pixelHeight]])
    # trans = np.array([[xOrigin, yOrigin]])
    # temp = np.array([[x, y]]).T - trans.T
    # res = np.dot(np.linalg.inv(mat), temp)
    # xPx = res[0].item()
    # yPx = res[1].item()
    #
    # if xPx < 0 or xPx > imageInfo.get("Dimension")[0]:
    #     warnings.warn("xPix outside the image dimension", DeprecationWarning, stacklevel=2)
    # if yPx < 0 or yPx > imageInfo.get("Dimension")[1]:
    #     warnings.warn("yPix outside the image dimension", DeprecationWarning, stacklevel=2)
    ############################################3
    # XMap = np.asarray(X)
    # YMap = np.asarray(Y)
    # XYMap = np.hstack((X, Y))
    # XY_sparse = sparse.csr_matrix(XYMap).T
    # # print(XY_sparse.shape,XMap.shape,YMap.shape)
    # if XMap.shape[0] != YMap.shape[0]:
    #     print("Erroor: x and Y must have the same size !!!")
    #     return
    # mat1 = sparse.kron(sparse.eye(XMap.shape[0]), [pixelWidth, rtnX])
    # mat2 = sparse.kron(sparse.eye(YMap.shape[0]), [rtnY, pixelHeight])
    # mat = sparse.vstack(blocks=[mat1, mat2])
    # # mat = sparse.csr_matrix(mat)
    # mat_invert = np.linalg.inv(mat.toarray())
    # mat_invert = sparse.csr_matrix(mat_invert)
    # # print(mat_invert.shape)
    # t1 = np.asarray([xOrigin] * XMap.shape[0])
    # t2 = np.asarray([yOrigin] * YMap.shape[0])
    # Tran = np.hstack((t1, t2))
    # Tran_sparse = sparse.csr_matrix(Tran).T
    # # print(Tran_sparse.shape)
    # res_sparse = sparse.csr_matrix.dot(mat_invert, XY_sparse - Tran_sparse)
    # res = sparse.csr_matrix.todense(res_sparse)
    # # print(res.tolist())
    # X_pix = res[0:XMap.shape[0]].tolist()
    # X_pix = np.concatenate(X_pix, axis=0)
    #
    # Y_pix = res[XMap.shape[0]:XMap.shape[0] + YMap.shape[0]].tolist()
    # # print(Y_map)
    # Y_pix = np.concatenate(Y_pix, axis=0)
    X_pix = []
    Y_pix = []

    for x, y in zip(X, Y):
        xPix, yPix = Map2Pixel(imageInfo=imageInfo, x=x, y=y)
        X_pix.append(xPix)
        Y_pix.append(yPix)
    return (X_pix, Y_pix)


def Map2Pixel_v2(imageInfo, X, Y):
    """
    :Method : convert coordinate from map space to image space,
        Note: The top Left coordinate of the image with GDAl correspond to (0,0)pix
    :param imagePath: path of the image  : string
    :param x: xMap coordinate : int or float as list 
    :param y: yMap coordinate: int or float as list 
    :return: coordinate in image space : tuple in pix
    """
    import numpy as np
    import scipy.sparse

    xOrigin = imageInfo["OriginCoord"][0]  # x top left coordinate (0.0 pix)
    yOrigin = imageInfo["OriginCoord"][1]  # y top left coordinate (0.0 pix)
    pixelWidth = imageInfo["Resolution"][0]  # xGSD
    pixelHeight = imageInfo["Resolution"][1]  # yGSD
    rtnX = imageInfo["MapInfo"][1][2]
    rtnY = imageInfo["MapInfo"][1][4]

    ## Apply inverse affine transformation
    nbObs = len(X)
    trans = np.array([[xOrigin], [yOrigin]])
    Trans = np.tile(trans, (nbObs, 1))
    L = np.empty((2 * nbObs, 1))
    for i in range(nbObs):
        L[2 * i, 0] = X[i]
        L[2 * i + 1, 0] = Y[i]
    Temp = L - Trans

    mat = np.array([[pixelWidth, rtnX], [rtnY, pixelHeight]])
    matInv = np.linalg.inv(mat)
    AA = scipy.sparse.block_diag((matInv,) * nbObs).toarray()

    Res = np.dot(AA, Temp)
    xPix = np.empty((nbObs, 1))
    yPix = np.empty((nbObs, 1))
    for i in range(nbObs):
        xPix[i] = Res[2 * i, 0]
        yPix[i] = Res[2 * i + 1, 0]
    # xPx = res[0].item()
    # yPx = res[1].item()
    # 
    # if xPx < 0 or xPx > imageInfo.get("Dimension")[0]:
    #     warnings.warn("xPix outside the image dimension", DeprecationWarning, stacklevel=2)
    # if yPx < 0 or yPx > imageInfo.get("Dimension")[1]:
    #     warnings.warn("yPix outside the image dimension", DeprecationWarning, stacklevel=2)

    return xPix, yPix


def XCoord_YCoord_as_Grid(rasterInfo):
    """
    Within this function we create an array which has the same dimension as the original raster.
    In each pixel we put as intensity value the XCoord (respectively the yCoord).
    :param rasterInfo: map information contains the projection and Geotrasform
    (the raster has to be opened with Routine function): dict
    :return:
        xCoord : 2D-Array contains the xCoord
        yCoord: 2D-Array contains the yCoord
    """
    xCoord = np.empty(rasterInfo.get("Dimension"))
    yCoord = np.empty(rasterInfo.get("Dimension"))
    for xPix in range(rasterInfo.get("Dimension")[0]):  # width
        for yPix in range(rasterInfo.get("Dimension")[1]):  # height
            xMap, yMap = Pixel2Map(imageInfo=rasterInfo, x=xPix, y=yPix)
            xCoord[xPix, yPix] = xMap
            yCoord[xPix, yPix] = yMap
    return xCoord, yCoord


def ConvCoordMap1ToMap2(x, y, targetEPSG, z=None, sourceEPSG=4326, display=False):
    """
    :Method: convert point coordinates from source to target system
    :param point x,y : map coordinate (e.g lon,lat)
    :param sourceEPSG: source coordinate system of the point; integer (default geographic coordinate )
    :param targetEPSG: target coordinate system of the point; integer
    :return: point in target coordinate system; list =[xCoord,yCoord,zCoord]

    Note: if the trasfoemation from WGS to UTM, x= lat, y=lon ==> coord =(easting(xMap) ,northing(yMap))
    """

    ## Set the source system

    source = osr.SpatialReference()  # instance from SpatialReference Class
    source.ImportFromEPSG(int(sourceEPSG))  # create a projection system based on EPSG code

    ## Set the target system
    target = osr.SpatialReference()
    target.ImportFromEPSG(int(targetEPSG))

    transform = osr.CoordinateTransformation(source,
                                             target)  # instance of the class Coordinate Transformation
    if z == None:
        coord = transform.TransformPoint(x, y)

    else:
        coord = transform.TransformPoint(x, y, z)
    if display:
        print("{}, {}, EPSG={} ----> Easting:{}  , Northing:{} , EPSG={}".format(x, y, sourceEPSG, coord[0], coord[1],
                                                                                 targetEPSG))
    return coord


def ConvCoordMap1ToMap2_Batch(X, Y, targetEPSG, Z=[], sourceEPSG=4326, display=False):
    """
    :Method: convert point coordinates from source to target system
    :param point x,y : map coordinate (e.g lon,lat)
    :param sourceEPSG: source coordinate system of the point; integer (default geographic coordinate )
    :param targetEPSG: target coordinate system of the point; integer
    :return: point in target coordinate system; list =[xCoord,yCoord,zCoord]

    Note: if the transformation from WGS to UTM, x= lat, y=lon ==> coord =(easting(xMap) ,northing(yMap))
    Note: if the transformation from UTM to WGS 84, x=easting, y=Notthing ==> lat, long
    """

    ## Set the source system
    # print(Z)
    import pyproj

    sourceEPSG_string = "epsg:" + str(sourceEPSG)
    targetEPSG_string = "epsg:" + str(targetEPSG)
    transformer = pyproj.Transformer.from_crs(sourceEPSG_string, targetEPSG_string)
    if len(Z) == 0:
        return transformer.transform(X, Y)

    else:
        return transformer.transform(X, Y, Z)


def ConvertGeo2Cartesian(lon, lat, alt):
    import pyproj
    # ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    # # print(ecef)
    # lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    # # print(lla)
    # # transformer = pyproj.Transformer.from_crs(lla, ecef)
    # x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)
    # print( [x,y,z])
    """https://pyproj4.github.io/pyproj/dev/api/proj.html"""
    transproj = pyproj.Transformer.from_crs("EPSG:4326", {"proj": 'geocent', "ellps": 'WGS84', "datum": 'WGS84'},
                                            always_xy=True)
    xpj, ypj, zpj = transproj.transform(lon, lat, alt, radians=False)
    return [xpj, ypj, zpj]


def ConvertGeo2Cartesian_Batch(Lon, Lat, Alt):
    import pyproj
    # ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    # # print(ecef)
    # lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    # # print(lla)
    # # transformer = pyproj.Transformer.from_crs(lla, ecef)
    # x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)
    # print( [x,y,z])
    """https://pyproj4.github.io/pyproj/dev/api/proj.html"""
    transproj = pyproj.Transformer.from_crs("EPSG:4326", {"proj": 'geocent', "ellps": 'WGS84', "datum": 'WGS84'},
                                            always_xy=True)
    Xpj, Ypj, Zpj = transproj.transform(Lon, Lat, Alt, radians=False)
    return Xpj, Ypj, Zpj


def ConvertCartesian2Geo(x, y, z):
    import pyproj
    # ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    # # print(ecef)
    # lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    # # print(lla)
    # # transformer = pyproj.Transformer.from_crs(lla, ecef)
    # x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)
    # print( [x,y,z])
    """https://pyproj4.github.io/pyproj/dev/api/proj.html"""
    transproj = pyproj.Transformer.from_crs({"proj": 'geocent', "ellps": 'WGS84', "datum": 'WGS84'}, "EPSG:4326",
                                            always_xy=True)
    lon, lat, alt = transproj.transform(x, y, z, radians=False)
    return [lon, lat, alt]


def ConvertCartesian2Geo_Batch(X, Y, Z):
    import pyproj
    # ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    # # print(ecef)
    # lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    # # print(lla)
    # # transformer = pyproj.Transformer.from_crs(lla, ecef)
    # x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)
    # print( [x,y,z])
    """https://pyproj4.github.io/pyproj/dev/api/proj.html"""
    transproj = pyproj.Transformer.from_crs({"proj": 'geocent', "ellps": 'WGS84', "datum": 'WGS84'}, "EPSG:4326",
                                            always_xy=True)
    Lon, Lat, Alt = transproj.transform(X, Y, Z, radians=False)
    return Lon, Lat, Alt


def ImgDims_Map_Pix(imageInfo):
    # imgDimsMap = {"x0Map": imageInfo.get("OriginCoord")[0],
    #               "xfMap": imageInfo.get("OriginCoord")[0] + imageInfo.get("Dimension")[0] *
    #                        imageInfo.get("Resolution")[0],
    #               "y0Map": imageInfo.get("OriginCoord")[1],
    #               "yfMap": imageInfo.get("OriginCoord")[1] - imageInfo.get("Dimension")[1] *
    #                        imageInfo.get("Resolution")[0]}
    import rasterio
    ds = rasterio.open(imageInfo.get("ImagePath"))
    bounds = ds.bounds
    imgDimsMap = {"x0Map": bounds.left,
                  "xfMap": bounds.right,
                  "y0Map": bounds.bottom,
                  "yfMap": bounds.top}
    x0, y0 = Map2Pixel(imageInfo=imageInfo, x=imgDimsMap.get("x0Map"), y=imgDimsMap.get("y0Map"))
    xf, yf = Map2Pixel(imageInfo=imageInfo, x=imgDimsMap.get("xfMap"), y=imgDimsMap.get("yfMap"))
    imgDimsPix = {"x0Pix": x0, "xfPix": xf, "y0Pix": y0, "yfPix": yf}

    return imgDimsMap, imgDimsPix


def WriteRaster(refRasterPath, newRasterPath, Listarrays, numberOfBands,
                descriptions=None,
                newGeoTransform=[],
                options=[],
                resample_alg=gdal.GRA_Bilinear,  # gdal.GRA_Lanczos,
                dtype=gdal.GDT_Float32,
                noData=None,
                metaData=[],
                RPC=None,
                GCPs=None, progress=False):
    """
    Generate raster with multi bands from arrays
    :param refRasterPath: Path of a raster to get projection information
    :param newRasterPath: Name and path of generated raster
    :param Listarrays: list of arrays
    :param numberOfBands: number of bands to store
    :param descriptions:  the description of each band
    :param metaData: list or dictionary of metadata
    :param GCPs :list [list og Gdal GCPS, GCP projection system]
    :param options :JPEG/LZW/PACKBITS/DEFLATE/CCITTRLE/CCITTFAX3/CCITTFAX4/LZMA/ZSTD/LERC/LERC_DEFLATE/LERC_ZSTD/WEBP/NONE
    :param resample_alg: GRA_NearestNeighbour, GRA_Bilinear, GRA_Cubic, GRA_CubicSpline, GRA_Lanczos,	GRA_Average,
                        GRA_Mode,
    :return: New Geotiff raster

    """
    if len(Listarrays) == numberOfBands:
        rasterInfo = GetRasterInfo(inputRaster=refRasterPath)
        raster = gdal.Open(refRasterPath)

        originX, originY = rasterInfo["OriginCoord"]
        pixelWidth, pixelHeight = rasterInfo["Resolution"]
        rows, cols = np.shape(Listarrays[0])

        ## Prepare destination file
        driver = gdal.GetDriverByName('GTiff')
        if options:
            print("Option:", options)
            outRaster = driver.Create(newRasterPath, cols, rows, numberOfBands, dtype, options)
        else:
            outRaster = driver.Create(newRasterPath, cols, rows, numberOfBands, dtype)

        ## Set the transformation
        if not newGeoTransform:
            outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        else:
            print("The new GeoTransformation will be set")
            outRaster.SetGeoTransform((newGeoTransform[0], newGeoTransform[1], 0, newGeoTransform[3], 0,
                                       newGeoTransform[5]))

        ## Set the projection
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outRaster.SetMetadataItem("Author", "SAIF AATI saif@caltech.edu")
        ## Set the metadata
        metaData_ = []
        if metaData:
            if isinstance(metaData, dict):
                for key, value in metaData.items():
                    temp = [key, value]
                    metaData_.append(temp)
            elif isinstance(metaData, list):
                metaData_ = metaData

            for mm in metaData_:
                # print("mm    ",mm)
                if not isinstance(mm[1], dict):
                    if not isinstance(mm[1], str):
                        str_ = str(mm[1])
                        outRaster.SetMetadataItem(mm[0], str_)
                    else:
                        outRaster.SetMetadataItem(mm[0], mm[1])
        if RPC:
            outRaster.SetMetadata(RPC, 'RPC')

        if GCPs:
            outRaster.SetGCPs(GCPs[0], GCPs[1])

        ## Write the data
        for i in range(numberOfBands):
            outband = outRaster.GetRasterBand(i + 1)
            outband.WriteArray(Listarrays[i], resample_alg=resample_alg)
            if noData != None:
                print("No data=", noData)
                outband.SetNoDataValue(noData)
            if descriptions != None:
                outband.SetDescription(descriptions[i])
            if progress:
                print("Writing band number: ", i + 1, " ", i + 1, "/", numberOfBands)

        outband.FlushCache()
    else:
        print("Error number of band does not correspond to number of array")

    return


def WriteRaster_maskValues(inputRaster, maxValue=50, noDataValue=-32767):
    rasterInfo = GetRasterInfo(inputRaster)
    arrayList = []
    for i in range(rasterInfo["NbBands"]):
        array_ = ImageAsArray(rasterInfo, i + 1)
        array_ = np.ma.masked_where(array_ < -maxValue, array_)
        array_ = np.ma.masked_where(array_ > maxValue, array_)
        array_ = array_.filled(fill_value=noDataValue)
        arrayList.append(array_)
    WriteRaster(refRasterPath=inputRaster,
                newRasterPath=os.path.join(os.path.dirname(inputRaster), Path(inputRaster).stem + "_filt.tif"),
                Listarrays=arrayList, numberOfBands=len(arrayList), descriptions=rasterInfo["BandInfo"],
                noData=noDataValue)

    return


# def Array2Raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):
#
#     cols = array.shape[1]
#     rows = array.shape[0]
#     originX = rasterOrigin[0]
#     originY = rasterOrigin[1]
#
#     driver = gdal.GetDriverByName('GTiff')
#     outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
#     outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
#     outband = outRaster.GetRasterBand(1)
#     outband.WriteArray(array)
#     outRasterSRS = osr.SpatialReference()
#     outRasterSRS.ImportFromEPSG(4326)
#     outRaster.SetProjection(outRasterSRS.ExportToWkt())
#     outband.FlushCache()


def DataMemory(dataType):
    gdalTypesLib = [gdal.GDT_Unknown, gdal.GDT_Byte, gdal.GDT_UInt16, gdal.GDT_Int16, gdal.GDT_UInt32, gdal.GDT_Int32,
                    gdal.GDT_Float32, gdal.GDT_Float64, gdal.GDT_CInt16, gdal.GDT_CInt32, gdal.GDT_CFloat32,
                    gdal.GDT_CFloat64]

    dataTypeLib = {"uint8": 1, "int8": 1, "uint16": 2, "int16": 3, "uint32": 4, "int32": 5,
                   "float32": 6, "float64": 7, "complex64": 10, "complex128": 1}
    dataTypeMemo = {"uint8": 8, "int8": 8, "uint16": 16, "int16": 16, "uint32": 32, "int32": 32,
                    "float32": 32, "float64": 64, "complex64": 64, "complex128": 128}
    for index, val in enumerate(dataTypeLib.items()):
        if dataType == val[1]:
            # print(index,val)
            memoType = dataTypeMemo.get(val[0])
    return memoType


def GetDataFromRaster(imageInfo, windDims, bandNumber):
    """
    :method: Extract data from source image on a specific band and window dimension
    :param imageInfo :dictionaire  string
    :param windDims: dimension of the window : List[xMin,XMax,yMin,yMax] in pixel
                    srcWin --- subwindow in pixels to extract: [left_x,width, top_y, height]
    :param bandNumber: band number : int
    :return: the data extracted form the image : array, bool True=outside of image
    """
    imageAsArray = ImageAsArray(imageInfo=imageInfo, bandNumber=bandNumber)
    ## Verify if the window around the tie point is inside the DEM image
    if (windDims[0] >= 0) and (windDims[1] <= imageInfo.get("Dimension")[0]) and \
            (windDims[2] >= 0) and (imageInfo.get("Dimension")[1] >= windDims[3]):
        windowData = imageAsArray[windDims[2]:windDims[3] + 1, windDims[0]:windDims[1] + 1]
        return windowData, False

    else:
        print("ERROR Window out of bounds")
        print(imageInfo.get("Dimension"), windDims)
        return 0, True


def ResizeRaster(refImg, inputImg, outputPath, option=2):
    """
    This module resize the input image to the refimage
    :param refImg: path of the image that has the small size
    :param inputImg: path of the image to be subseted
    :param savingPath: path and the name of the outputimage
    :param option 1 :using routine method where we extract the array
                  2 : using Gdal_translate where we save directly the ouput
    :return:
    """
    refImgInfo = GetRasterInfo(inputRaster=refImg)
    inputImgInfo = GetRasterInfo(inputRaster=inputImg)
    metaData = inputImgInfo["MetaData"]
    metaData["RefRaster_for_resize"] = refImg

    if option == 1:
        windDim_ = [0, refImgInfo.get("Dimension")[0], 0, refImgInfo.get("Dimension")[
            1]]  # srcWin --- subwindow in pixels to extract: [left_x,width, top_y, height]
        newData, verif = GetDataFromRaster(imageInfo=inputImgInfo, windDims=windDim_, bandNumber=1)
        metaData["Resize method"] = " Using option 1; based on routine Writeraster "
        ## Save the new data into raster
        WriteRaster(refRasterPath=inputImg, newRasterPath=outputPath, Listarrays=[newData], numberOfBands=1)
        return 0
    if option == 2:
        xMin, yMin = Pixel2Map(imageInfo=refImgInfo, x=0, y=0)
        xMax, yMax = Pixel2Map(imageInfo=refImgInfo, x=refImgInfo.get("Dimension")[0], y=refImgInfo.get("Dimension")[1])
        windDim = [xMin, yMin, xMax,
                   yMax]  # projWin --- subwindow in projected coordinates to extract: [ulx, uly, lrx, lry]
        gdal.Translate(outputPath, inputImgInfo.get("Raster"), projWin=[windDim[0], windDim[1], windDim[2], windDim[3]])
        metaData["Resize method"] = " Using option 2; based on gdal translate"
        raster = GetRasterInfo(outputPath)["Raster"]
        raster.SetMetadataItem("Author", "SAIF AATI:saif@caltech.edu")
        metaData_ = []
        if metaData:
            if isinstance(metaData, dict):
                for key, value in metaData.items():
                    temp = [key, value]
                    metaData_.append(temp)
            elif isinstance(metaData, list):
                metaData_ = metaData

            for mm in metaData_:
                if not isinstance(mm[1], dict):
                    if not isinstance(mm[1], str):
                        str_ = str(mm[1])
                        raster.SetMetadataItem(mm[0], str_)
                    else:
                        raster.SetMetadataItem(mm[0], mm[1])
        raster.FlushCache()
        return 0


def SubsetImg(inputImg, outputPath, windDim=[], metaData=None, bandList=[1]):
    """
    This module subset the input image, using the dimension windDim
    :param inputImg: path of the image to be subseted
    :param windDim : list of map coordinates projWind =[xMinMap,yMinMap,xMaxMap,yMapMap] + in projected coordinates
    :param savingPath: path and the name of the outputimage

    :return:
    """

    inputImgInfo = GetRasterInfo(inputRaster=inputImg)
    gdal.Translate(outputPath, inputImgInfo.get("Raster"), projWin=[windDim[0], windDim[1], windDim[2], windDim[3]],
                   format="GTiff", outputType=gdal.GDT_Float64, bandList=bandList)
    # inputdataset = gdal.Open(outputPath)
    # inputdataset.SetMetadataItem("Author", "SAIF AATI:saif@caltech.edu")
    #
    # if metaData:
    #     for mm in metaData:
    #         inputdataset.SetMetadataItem(mm[0], mm[1])
    #
    # inputdataset.FlushCache()

    return 0


def NormalizeImage(img):
    """
    A better to normalize the image is to take each value and divide by the largest value experienced by the data type.
    This ensures that images that have a small dynamic range in the image remain small and they are not inadvertently
    normalized so that they become gray.
    we want to have the same intensity scale after converting to np.uint8
    :param img: as array
    :return:
    """

    max = img.max()
    # print("max=", max)
    imgNorm = img.astype(np.float64) / max  # normalize the data to 0 -1
    imgScale = 255 * imgNorm  # scale the image by 255

    imgFinal = imgScale.astype(np.uint8)
    return (imgFinal)


def Interpolation(x, y, mode="quadratic"):
    """

    :param x:
    :param y:
    :param mode:
    :return:
    """
    from scipy.interpolate import interp1d

    if mode == "linear":
        return interp1d(x=x, y=y, kind="linear"), mode
    elif mode == "nearest":
        # return the nearest point along the x-axis
        return interp1d(x=x, y=y, kind="nearest"), mode
    elif mode == "cubic":
        # Spline interpolation 3rd order
        return interp1d(x=x, y=y, kind="cubic"), mode
    elif mode == "quadratic":
        # Spline interpolation 2nd order
        return interp1d(x=x, y=y, kind="quadratic"), mode
    elif mode == "slinear":
        # Spline interpolation 1st order
        return interp1d(x=x, y=y, kind="slinear", assume_sorted=False), mode
    elif mode == "previous":
        # # return the previous point along the x-axis
        return interp1d(x=x, y=y, kind="previous"), mode
    elif mode == "next":
        # # return the next point along the x-axis
        return interp1d(x=x, y=y, kind="next"), mode
    elif mode == "gekko":
        from gekko import GEKKO
        m = GEKKO()
        # m.x = m.Param(value=tNew)
        # m.y = m.Var()
        #
        # m.options.IMODE = 2
        # m.cspline(m.x, m.y, listDatesDays, V_nbPC)
        # m.solve(disp=False)
        # ax2.plot(m.x.value, m.y.value, 'r--', label='cubic spline')


def Create_MultiBandsRaster_form_MultiRasters(rastersList, outputPath, bandNumber=1):
    """

    :param rastersFolder:
    :param outputPath:
    :return:
    """
    arrayList = []
    bandDescription = []
    for img_ in rastersList:
        array = ImageAsArray(imageInfo=GetRasterInfo(img_), bandNumber=bandNumber)
        arrayList.append(array)
        bandDescription.append("Band" + str(bandNumber) + "_" + Path(img_).stem)
    WriteRaster(refRasterPath=rastersList[0], newRasterPath=outputPath, Listarrays=arrayList,
                numberOfBands=len(rastersList), descriptions=bandDescription)
    return


def Create_MultiRasters_from_MultiBandsRaster(inputRaster, output):
    """

    :param inputRaster:
    :param output:
    :return:
    """
    rasterInfo = GetRasterInfo(inputRaster)
    print(rasterInfo["NbBands"])
    for i in range(rasterInfo["NbBands"]):
        array = ImageAsArray(rasterInfo, i + 1)
        WriteRaster(refRasterPath=inputRaster, newRasterPath=os.path.join(output, "Img_Band_" + str(i + 1) + ".tif"),
                    Listarrays=[array], numberOfBands=1)
    return


def Gdal_rasterize(refRaster, shpFile, output):
    ndsm = refRaster
    shp = shpFile
    data = gdal.Open(ndsm, gdalconst.GA_ReadOnly)
    geo_transform = data.GetGeoTransform()
    # source_layer = data.GetLayer()
    x_min = geo_transform[0]
    y_max = geo_transform[3]
    x_max = x_min + geo_transform[1] * data.RasterXSize
    y_min = y_max + geo_transform[5] * data.RasterYSize
    x_res = data.RasterXSize
    y_res = data.RasterYSize

    mb_v = ogr.Open(shp)
    mb_l = mb_v.GetLayer()
    pixel_width = geo_transform[1]
    target_ds = gdal.GetDriverByName('GTiff').Create(output, x_res, y_res, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform((x_min, pixel_width, 0, y_min, 0, pixel_width))
    band = target_ds.GetRasterBand(1)
    NoData_value = -999999
    band.SetNoDataValue(NoData_value)
    band.FlushCache()
    gdal.RasterizeLayer(target_ds, [1], mb_l)

    target_ds = None
    return


def ConvertRaster2WGS84(inputPath, outputPath=None):
    if outputPath == None:
        outputPath = Path(inputPath).stem + "_conv_4326.tif"

    gdal.Warp(outputPath, inputPath, dstSRS='EPSG:4326')
    return outputPath


def ExtractBandFromRaster(inputPath, bandNumber, outputPath=None):
    if not outputPath:
        outputPath = os.path.dirname(inputPath)

    baseName = os.path.basename(inputPath)
    srcDs = GetRasterInfo(inputRaster=inputPath)["Raster"]
    out_ds = gdal.Translate(os.path.join(outputPath, baseName[:-4] + "_band" + str(bandNumber) + '.tif'), srcDs,
                            format='GTiff',
                            bandList=[bandNumber])
    out_ds = None

    return


def BatchExtractBandFromRaster(inputPath, outputPath=None):
    from FilesCommandRoutine import GetFilesBasedOnExtension
    imgsList = GetFilesBasedOnExtension(path=inputPath)
    if not outputPath:
        # ExtractBandFromRaster(inputPath=path, bandNumber=3)
        return
    else:
        for imgPath in imgsList:
            ExtractBandFromRaster(inputPath=imgPath, outputPath=outputPath, bandNumber=3)

    return


def BuildRasterMask(refImg, shpInput, output=None):
    import geospatialroutine.Routine_Lyr as VecRT
    # read raster file
    rasterInfo = GetRasterInfo(refImg)
    rasterDim = rasterInfo.get("Dimension")
    ## Read Shape file
    shpVector = ogr.Open(shpInput)
    shpLayer = shpVector.GetLayer()
    print(shpLayer)
    vectorInfo = VecRT.GetVectorInfo(shpInput=shpInput)

    ## Generate output grid
    # outputGridArray = np.zeros(rasterDim)

    # import Vector_routine as vector
    res = VecRT.RasterizeVector(shpInput=shpInput, refRaster=refImg, output=output)
    # res = None
    print(res)
    return res


def ParseResampMethod(method):
    """
    Parse resampling method
    :param method:
    :return:
    """
    if method == 'near':
        # Note: Nearest respects nodata when downsampling
        return gdal.GRA_NearestNeighbour
    elif method == 'bilinear':
        return gdal.GRA_Bilinear
    elif method == 'cubic':
        return gdal.GRA_Cubic
    elif method == 'cubicspline':
        return gdal.GRA_CubicSpline
    elif method == 'average':
        return gdal.GRA_Average
    elif method == 'lanczos':
        return gdal.GRA_Lanczos
    elif method == 'mode':
        # Note: Mode respects nodata when downsampling, but very slow
        return gdal.GRA_Mode
    elif method == 'sinc':
        print("Not implemented yet, it is the sinc method implmented in cosi-corr")
        sys.exit("Invalid resampling method")
    else:
        algo = None
        sys.exit("Invalid resampling method")


def ComputeEpsg(lon, lat):
    """
    Compute the EPSG code of the UTM zone which contains
    the point with given longitude and latitude

    Args:
        lon (float): longitude of the point
        lat (float): latitude of the point

    Returns:
        int: EPSG code
    """
    # UTM zone number starts from 1 at longitude -180,
    # and increments by 1 every 6 degrees of longitude
    zone = int((lon + 180) // 6 + 1)

    # EPSG = CONST + ZONE where CONST is
    # - 32600 for positive latitudes
    # - 32700 for negative latitudes
    const = 32600 if lat > 0 else 32700
    # print("const=",const,"  zone=",zone)
    return const + zone


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


def Set_crs(epsg=4326):
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(epsg)
    # print("-----", outRasterSRS)

    return outRasterSRS


def Array2Raster(rasterSavingPath, rasterOrigin, pixelWidth, pixelHeight, array, epsg=4326):
    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(rasterSavingPath, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


class Distribution:
    def __init__(self, inputArray, displayValue=True):
        sample = np.ma.masked_invalid(inputArray)
        mask = np.ma.getmask(sample)

        # Remove mask and array to vector
        if isinstance(sample, np.ma.MaskedArray):  # check if the sample was masked using the class numpy.ma.MaskedArray
            sample = sample.compressed()  ## return all the non-masked values as 1-D array
        else:
            if sample.ndim > 1:  # if the dimension of the array more than 1 transform it to 1-D array
                sample = sample.flatten()
        self.sample = sample
        # Estimate initial sigma and RMSE
        (self.mu, self.sigma) = norm.fit(sample)
        self.sigma_ = '%.3f' % (self.sigma)
        temp = np.square(sample)
        temp = np.ma.masked_where(temp<=0,temp)
        self.RMSE = '%.3f' % (np.ma.sqrt(np.ma.mean(temp)))

        self.max = '%.3f' % (np.nanmax(sample))
        self.min = '%.3f' % (np.nanmin(sample))
        self.std = '%.3f' % (np.nanstd(sample))
        self.mean = '%.3f' % (np.nanmean(sample))
        self.median = '%.3f' % (np.nanmedian(sample))
        self.mad = '%.3f' % (stats.median_absolute_deviation(sample))
        self.nmad = '%.3f' % (1.4826 * stats.median_absolute_deviation(sample))
        self.ce68 = stats.norm.interval(0.68, loc=self.mu, scale=self.sigma)
        self.ce90 = stats.norm.interval(0.9, loc=self.mu, scale=self.sigma)
        self.ce95 = stats.norm.interval(0.95, loc=self.mu, scale=self.sigma)
        self.ce99 = stats.norm.interval(0.99, loc=self.mu, scale=self.sigma)


        if displayValue:
            print("mu, sigma", self.mu, self.sigma)
            print("RMSE=", self.RMSE)
            print("max=", self.max, "min=", self.min, "std=", self.std, "mean=", self.mean, "median", self.median,
                  "mad=",
                  self.mad, "nmad=", self.nmad)
            print("CE68=", self.ce68, "CE90=", self.ce90, "CE95=", self.ce95, "CE99=", self.ce99)


def Merge(inputList, output, dtype="Float32"):
    """
    ### Using sub process
    cmd = ["gdal_merge.py"]

    cmd.extend(["-o", output])
    cmd.extend(["-ot "+dtype])
    cmd.extend(inputList)
    print(cmd)
    call = ""
    for cmd_ in cmd:
        call += cmd_ + " "
    print(call)
    os.system(call)
    return
    """

    print("################################")
    print("      Merge using Rasterio      ")
    print("################################")
    src_files_to_mosaic = []
    from rasterio.merge import merge
    for fp in inputList:
        src = rasterio.open(fp)
        src_files_to_mosaic.append(src)
    # Merge function returns a single mosaic array and the transformation info
    mosaic, out_trans = merge(src_files_to_mosaic)
    print(mosaic.shape)
    #### Copy the metadata
    out_meta = src.meta.copy()
    crs = src.crs
    print(crs)
    # Update the metadata
    out_meta.update({"driver": "GTiff", "height": mosaic.shape[1], "width": mosaic.shape[2], "transform": out_trans,
                     "crs": crs})  # "+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs "
    rasterTemp = GetRasterInfo(inputRaster=inputList[0])
    descriptions = rasterTemp["BandInfo"]
    listArrays = []
    for id in range(mosaic.shape[0]):
        with rasterio.open(output, "w", **out_meta) as dest:
            # the .astype(rasterio.int16) forces dtype to int16
            dest.write_band(id + 1, mosaic[id, :, :])
            dest.set_band_description(id + 1, descriptions[id])
        listArrays.append(mosaic[id, :, :])

    WriteRaster(refRasterPath=output, newRasterPath=output, Listarrays=listArrays, numberOfBands=mosaic.shape[0],
                descriptions=descriptions)


if __name__ == '__main__':
    path = "/media/stoage/Saif/3D_approch/NAIP_images"
    import geospatialroutine.FilesCommandRoutine as File_RT

    # imgList = File_RT.GetFilesBasedOnExtension(path="/home/cosicorr/0-WorkSpace/3D-Correlation_project/6p4_3D/Correlation/Set1/3D_correlation_results")
    output = "//home/cosicorr/0-WorkSpace/3D-Correlation_project/6p4_3D/Correlation/Set1/3D_Displacement_Python_Mosaic.tif"
    # Mosaic(inputList=imgList,output=output)
    # rasterInfo = GetRasterInfo(output, printInfo=True)
    # WriteRaster(refRasterPath=output, newRasterPath="/media/storage/Saif/3D_approch/NAIP_RGB_Mosaic.tif",
    #             Listarrays=[ImageAsArray(rasterInfo, bandNumber=1), ImageAsArray(rasterInfo, bandNumber=2),
    #                         ImageAsArray(rasterInfo, bandNumber=3)],numberOfBands=3,dtype=gdal.GDT_Byte)
    #
    # imgList = File_RT.GetFilesBasedOnExtensions(
    #     "//home/cosicorr/0-WorkSpace/3D-Correlation_project/TestPCA_ICA/ENVI/Dz/ICA/Gaussian/ICA_reconst_1_3/PCA/PCA")
    path = "//home/cosicorr/0-WorkSpace/Test_Data_ElMayor/Correlation_Set2_test/Test8"
    imgList = File_RT.GetFilesBasedOnExtensions(path)
    print(imgList)
    # Create_MultiBandsRaster_form_MultiRasters(rastersList=imgList,
    #                                           outputPath=os.path.join(path, "Set2_Test8_3DDA_Dz.tif"), bandNumber=3)

    # Create_MultiRasters_from_MultiBandsRaster(
    #     inputRaster="//home/cosicorr/0-WorkSpace/3D-Correlation_project/TestPCA_ICA/ENVI/Dz/ICA/Gaussian/ICA_reconst_1_3/allDz-6-DZ_ICA_Gaussian_Reconstwith_PC1_3.tif",
    #     output="/home/cosicorr/0-WorkSpace/3D-Correlation_project/TestPCA_ICA/ENVI/Dz/ICA/Gaussian/ICA_reconst_1_3/PCA")

    # 
    path = "/home/cosicorr/0-WorkSpace/Test_Data_ElMayor/Correlation_Set2_test/Test4/Test4-4_3DDA/3DDTiles"
    imgList = File_RT.GetFilesBasedOnExtensions(path)
    Merge(inputList=imgList,
          output="//home/cosicorr/0-WorkSpace/Test_Data_Elayor/Correlation_Set2_test/Test4/Test_4_4_3DD.tif")
