import sys

import numpy as np
from scipy import interpolate
from scipy.interpolate import *
import geospatialroutine.Routine as RT
import warnings

from p_tqdm import p_map

class BicubicInterpolation:

    def __init__(self, rasterInfo, x, y):

        self.imageInfo = rasterInfo
        self.bandNumber = 1
        self.res, self.validation = self.BicubicInterpolation_(x=x, y=y)

    def DFfunction(self, x):
        if np.abs(x) < 1:
            val = (np.abs(x) ** 3) - 2 * (np.abs(x) ** 2) + 1
            return val
        else:
            if (np.abs(x) < 2) and (np.abs(x) >= 1):
                val = -np.abs(x) ** 3 + 5 * (np.abs(x) ** 2) - 8 * np.abs(x) + 4
                return val
            else:
                return 0

    def Afunction(self, i, j, n, dx):

        s1, verif1 = self.Sfunction(i - 1, j + n - 2)
        s2, verif2 = self.Sfunction(i, j + n - 2)
        s3, verif3 = self.Sfunction(i + 1, j + n - 2)
        s4, verif4 = self.Sfunction(i + 2, j + n - 2)
        if verif1 == False and verif2 == False and verif3 == False and verif4 == False:
            a_n = s1 * self.DFfunction(dx + 1) + s2 * self.DFfunction(dx) + \
                  s3 * self.DFfunction(dx - 1) + s4 * self.DFfunction(dx - 2)
            return a_n, False
        else:
            return 0, True

    def Sfunction(self, i, j):

        s, verif = RT.GetDataFromRaster(imageInfo=self.imageInfo, windDims=[i, i, j, j], bandNumber=self.bandNumber)
        # data = np.array([[20, 41, 35, 20], [15, 12, 21, 24], [11, 25, 32, 30], [8, 25, 35, 50]])
        # s = data [j,i]
        if verif == False:
            return s.item(), verif
        else:
            return s, verif

    def BicubicInterpolation_(self, x, y):
        i = int(x - 0.5)
        j = int(y - 0.5)
        dx = (x - 0.5) - i
        dy = (y - 0.5) - j
        a1, verif1 = self.Afunction(i, j, 1, dx)
        a2, verif2 = self.Afunction(i, j, 2, dx)
        a3, verif3 = self.Afunction(i, j, 3, dx)
        a4, verif4 = self.Afunction(i, j, 4, dx)

        if verif1 == False and verif2 == False and verif3 == False and verif4 == False:
            res = a1 * self.DFfunction(dy + 1) + a2 * self.DFfunction(dy) + \
                  a3 * self.DFfunction(dy - 1) + a4 * self.DFfunction(dy - 2)
            return res, False
        else:
            return 0, True



def InterpolationMethod(method):
    """
     Function to convert number into string
    # Switcher is dictionary data type here
    # get() method of dictionary data type returns
    # value of passed argument if it is present
    # in dictionary otherwise second argument will
    # be assigned as default value of passed argument
    :param method:
    :return:
    """
    switcher = {
        0: "nearest",
        1: "linear",
        2: "cubic",
    }

    return switcher.get(method, "default")


def Interpolation1d(x, y, mode="quadratic"):
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
        print("To add")
        # from gekko import GEKKO
        # m = GEKKO()
        # m.x = m.Param(value=tNew)
        # m.y = m.Var()
        #
        # m.options.IMODE = 2
        # m.cspline(m.x, m.y, listDatesDays, V_nbPC)
        # m.solve(disp=False)
        # ax2.plot(m.x.value, m.y.value, 'r--', label='cubic spline')


def InterpolationCosiCorr(xTemp, yTemp, demInfo, tiePointIndex, interpolationType):
    """
    :param xTemp: xTiepoints in demImage space
    :param yTemp: yTiepoints in demImage space
    :param demInfo: list of the dem image
    :param tiePointIndex: Number of the tie point
    param interpolationType : int
            : 0= "nearest"    zero order interpolation
            : 1= "bilinear"
            : 2= "bicubic"
            : 3= "spline"   :not yet implemented
            : Default = "bilinear"
    :return: altitude : float
    """
    ## Define window around this pixel
    ##:param windowDim: window dimension: array[X1, X2, Y1, Y2]
    windowDim = [int(xTemp) - 3, int(xTemp) + 3, int(yTemp) - 3, int(yTemp) + 3]

    ## Verify if the window around the tie point is onside the DEM image
    res, verif = RT.GetDataFromRaster(imageInfo=demInfo, windDims=windowDim, bandNumber=1)
    if verif == False:
        # Get the window data in this window
        windowData = res
        # print("windowData around TiePts=\n", windowData)

        h = Interpolate(data=windowData, points=(xTemp - (int(xTemp) - 3), yTemp - (int(yTemp) - 3)),
                        method=interpolationType)
        print("h=", h, "type  ", InterpolationMethod(interpolationType), '\n')
        return h

    else:
        print("Tie Point number ", tiePointIndex + 1, " is outside DEM coverage - Altitude set to 0.0")
        h = 0
        return h


def InterpolationSA(x, y, demInfo, interpolationType, bandNumber=1):
    """
    Interpolation using non-integer pixel position
    3 methods of interpolation have been implemented:
        1- zero order interpolation (nearest neighbour)
        2- first order interpolation (bilinear interpolation)
        3- second order interpolation (bicubic convolution, Lagrange polynomials)
    The difference between this implementation and the fo;mer used in cosi-corr or by numpy that we take
    retrieve the suitable altitude value from the local neighborhoods using non-integer pixel position
    :param xTemp: xTiepoints in demImage space
    :param yTemp: yTiepoints in demImage space
    :param demInfo: list of the dem image
    :return: altitude : float
    """
    print(" ### SA Interpolation non-integer pixel position, method: ", interpolationType)

    if interpolationType == "nearest":
        #### Zero Order Interpolation  (Nearest)
        xTemp_ = int(np.rint(x))
        yTemp_ = int(np.rint(y))
        point = [xTemp_, xTemp_, yTemp_, yTemp_]

        hNearest, verifNearest = RT.GetDataFromRaster(imageInfo=demInfo, windDims=point, bandNumber=bandNumber)

        if verifNearest == False:
            print("h_Nearest=", hNearest.item())
            return hNearest
        else:
            warnings.warn("Point number is outside DEM coverage - Altitude set to nan")
            h = np.nan
            return h

    if interpolationType == "default" or interpolationType == "linear":
        #### First Order Interpolation  (bilineare Interpolation)
        i_ = int(x - 0.5)
        j_ = int(y - 0.5)
        dx = (x - 0.5) - i_
        dy = (y - 0.5) - j_

        g1, verif1 = RT.GetDataFromRaster(imageInfo=demInfo, windDims=[i_, i_, j_, j_], bandNumber=bandNumber)
        g2, verif2 = RT.GetDataFromRaster(imageInfo=demInfo, windDims=[i_ + 1, i_ + 1, j_, j_], bandNumber=bandNumber)
        g3, verif3 = RT.GetDataFromRaster(imageInfo=demInfo, windDims=[i_, i_, j_ + 1, j_ + 1], bandNumber=bandNumber)
        g4, verif4 = RT.GetDataFromRaster(imageInfo=demInfo, windDims=[i_ + 1, i_ + 1, j_ + 1, j_ + 1],
                                          bandNumber=bandNumber)

        if (verif1 == False and verif2 == False and verif3 == False and verif4 == False):
            hBilinear = g1 + dx * (g2 - g1) + dy * (g3 - g1) + dx * dy * ((g4 - g2) - (g3 - g1))
            print("h_Bilinear= ", hBilinear.item())
            return hBilinear
        else:
            warnings.warn("Point number is outside DEM coverage - Altitude set to nan")
            h = np.nan
            return (h)

    if interpolationType == "cubic":
        #### Second Order Interpolation  (Bicubic Interpolation)
        bic = BicubicInterpolation(rasterInfo=demInfo, x=x, y=y)

        hBicubic, verif = bic.res, bic.validation
        if verif == False:
            print("hBicubic= ", hBicubic)
            return (hBicubic)
        else:
            warnings.warn("Point number is outside DEM coverage - Altitude set to nan")
            h = np.nan
            return (h)


def Interpolate(data, points, method):
    """
    :param data: data information: array
    :param points: Coordinate of the point :list=[x,y]
    :param method: interpolation method : int
                 : 0= " nearest"    zero order interpolation
                 : 1= "linear"
                 : 2= "cubic"
                 : 3= "spline"
                 : Default = "linear"
    :return: interpolated value : float
    """
    row, col = np.shape(data)
    coordList = []
    dataList = []
    # Transform data to list and create correspondent index array coordinate list
    for i in range(row):
        for j in range(col):
            coordList.append((i, j))
            dataList.append((data[i, j]))

    coorList = np.asanyarray(coordList)

    if InterpolationMethod(method=method) == "default":

        methodType = "linear"
        h = interpolate.griddata(coorList, dataList, points, method=methodType)
    else:
        if InterpolationMethod(method=method) != "spline":
            methodType = InterpolationMethod(method=method)
            h = interpolate.griddata(coorList, dataList, points, method=methodType)
            h = h.item()
        else:
            spline = " implement the code"
            methodType = "spline"
            print("Not yet implemented")
            # print(coorList[:,0])
            # h = interpolate.bisplrep(x=coorList[:,0], y=coorList[:,1], z=np.asarray(dataList), s=0)
            # print("h_spline=",h[2].shape)
            # h=np.asarray(h[2])
            # h=np.reshape(h,newshape=(7,7))

    return (h)


################## Interpolation    ###############################################

def Interpolate2D(inArray, x, y, kind='cubic'):
    """
    Procedure for 2D interpolation
 
    :param inArray:
    :param x: list
    :param y: list
    :param kind: RectBivariateSpline,quintic
    :return:
    """
    ## Note: add the possibility to perform sinc iinterpolation :
    ## Hints: https://stackoverflow.com/questions/1851384/resampling-interpolating-matrix
    ## https://www.programcreek.com/python/example/58004/numpy.sinc
    ## https://gist.github.com/gauteh/8dea955ddb1ed009b48e
    ## https://gist.github.com/endolith/1297227
    ## https://www.harrisgeospatial.com/docs/INTERPOLATE.html

    shape = np.shape(inArray)
    lin = shape[0]
    col = shape[1]
    # print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    # print(inArray.shape,x,y)
    if kind == "RectBivariateSpline":
        f = RectBivariateSpline(np.arange(0, lin, 1), np.arange(0, col, 1), inArray, kx=3, ky=3, s=0)
        if len(x) > 1 and len(y) > 1:
            # res = []
            # res = p_map(_2Dinterp, len(x)*[f], x, y,num_cpus=50)
            for x_, y_ in zip(x, y):
                res.append(f(x_, y_).item())

            return res
        else:
            # print(f(x,y))
            return [f(x, y).item()]
    if kind == "cubic":
        f = interp2d(np.arange(0, lin, 1), np.arange(0, col, 1), np.ndarray.flatten(inArray), kind=kind)
        if len(x) > 1 and len(y) > 1:
            res = []
            # res = p_map(_2Dinterp, len(x)*[f], x, y,num_cpus=50)
            for x_, y_ in zip(x, y):
                res.append(f(x_, y_).item())

            return res
        else:
            # print(f(x,y))
            return [f(x, y).item()]
    if kind == "linear" or kind == "nearest":
        #'linear in 1d = binlinear in 2D'
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
        from scipy.interpolate import RegularGridInterpolator
        f = interpolate.RegularGridInterpolator((np.arange(0, lin, 1), np.arange(0, col, 1)), inArray,
                                                method=kind,
                                                bounds_error=False, fill_value= np.nan)
        return f(list(zip(x, y)))


def Interpolate2D_batch(inArray, x_Array, y_Arrary, kind='cubic'):
    nbCpu = 50
    if len(x_Array) > nbCpu:
        print(len(x_Array) / nbCpu)
        sys.exit()

    return


def Nearest_2D(inarray, lin, col):
    """

    :param inarray:
    :param lin:
    :param col:
    :return:
    """
    lin_ = int(np.rint(lin))
    col_ = int(np.rint(col))
    return inarray[lin_, col_]


def Bilinear_2D(inarray, lin, col):
    """
    The bilinear or biquadratic interpolation takes into account the 2 x 2 adjacent grey values of the computed pixels.
    The interpolated gray value is the result of the weighted average of the adjacent grey values in witch the weight
    is given by  the relative coverage of te current pixel.

    Note to take into account of the edges !!!!
    :param inarray:
    :param lin:
    :param col:
    :return:
    """
    i = int(col - 0.5)
    j = int(lin - 0.5)

    dx = col - 0.5 - i
    dy = lin - 0.5 - j

    res = inarray[j, i] + dx * (inarray[j, i + 1] - inarray[j, i]) + dy * (
            inarray[j + 1, i] - inarray[j, i]) + dx * dy * (
                  inarray[j, i] + inarray[j + 1, i + 1] - inarray[j, i + 1] - inarray[j + 1, i])
    return res

def Lagrange_2D(inarray, lin, col):
    def Dfunction(x):
        if np.abs(x) < 1:
            val = (np.abs(x) ** 3) - 2 * (np.abs(x) ** 2) + 1
            return val
        elif (np.abs(x) < 2) and (np.abs(x) >= 1):
            val = -np.abs(x) ** 3 + 5 * (np.abs(x) ** 2) - 8 * np.abs(x) + 4
            return val
        else:
            return 0

    def Afunction(i, j, n, dx):

        s1 = inarray[j + n - 2, i - 1]
        s2 = inarray[j + n - 2, i]
        s3 = inarray[j + n - 2, i + 1]
        s4 = inarray[j + n - 2, i + 2]

        a_n = s1 * Dfunction(dx + 1) + s2 * Dfunction(dx) + s3 * Dfunction(dx - 1) + s4 * Dfunction(dx - 2)
        return a_n

    i = int(col - 0.5)
    j = int(lin - 0.5)
    dx = col - 0.5 - i
    dy = lin - 0.5 - j
    a1 = Afunction(i, j, 1, dx)
    a2 = Afunction(i, j, 2, dx)
    a3 = Afunction(i, j, 3, dx)
    a4 = Afunction(i, j, 4, dx)
    res = a1 * Dfunction(dy + 1) + a2 * Dfunction(dy) + a3 * Dfunction(dy - 1) + a4 * Dfunction(dy - 2)
    return res

def Interpoate1D(X, Y, xCord, kind="linear"):
    # print(X,Y)
    f = interp1d(X, Y, kind=kind, fill_value="extrapolate")
    ynew = f(xCord)  # use interpolation function returned by `inter

    if len(ynew) > 1:
        return ynew
    else:
        return [ynew.item()]

def LinearIterpolation(array, location):

    loc = int(location)
    tmpArr = array[loc]
    if loc + 1 < len(array):
        res = (array[loc + 1] - tmpArr) * (location - loc) + tmpArr
    else:
        res = (array[-1] - tmpArr) * (location - loc) + tmpArr
    return res


## ================================================================================================================= ##

# def Fast2Dinterp(inArray, X, Y):
#     from fast_interp import interp2d
#
#     # https://github.com/CD3/libInterpolate
#     # https://github.com/dbstein/fast_interp
#
#     shape = np.shape(inArray)
#     lin = shape[0]
#     col = shape[1]
#     f = interp2d([0, 0], [lin, col], [1, 1], inArray)  # , k=5, p=[False, True], e=[1, 0])
#
#     def _2Dinterp(f, xVal, yVal):
#
#         return f(xVal, yVal).item()
#
#     if len(X) > 1 and len(Y) > 1:
#         res = []
#         res = p_map(_2Dinterp, len(X) * [f], X, Y, num_cpus=50)
#         # for x_, y_ in zip(x, y):
#         #     res.append(f(x_, y_).item())
#
#         return res
#     else:
#         # print(f(x,y))
#         return [f(X, Y).item()]
