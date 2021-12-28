# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2021
import sys

import numpy as np
from scipy import interpolate
from scipy.interpolate import *

import warnings


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
        ##FIXME
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


def InterpolationMethods(method):
    """
    Function to convert number into string.
    Switcher is dictionary data type here
    get() method of dictionary data type returns
    value of passed argument if it is present
    in dictionary otherwise second argument will
    be assigned as default value of passed argument
    Args:
        method:

    Returns:

    """

    switcher = {
        0: "nearest",
        1: "linear",
        2: "cubic",
    }

    return switcher.get(method, "default")


def Interpolation1d(x, y, mode="quadratic"):
    """
    1D interpolation using Scipy
    Args:
        x:
        y:
        mode: linear,nearest,cubic,quadratic,previous,next,slinear

    Returns:

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
        ##TODO
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

    if InterpolationMethods(method=method) == "default":

        methodType = "linear"
        h = interpolate.griddata(coorList, dataList, points, method=methodType)
    else:
        if InterpolationMethods(method=method) != "spline":
            methodType = InterpolationMethods(method=method)
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
    Args:
        inArray:
        x: list
        y: list
        kind: RectBivariateSpline,quintic

    Returns:
    Notes:
            Add the possibility to perform sinc interpolation :
                Hints:  https://stackoverflow.com/questions/1851384/resampling-interpolating-matrix
                        https://www.programcreek.com/python/example/58004/numpy.sinc
                        https://gist.github.com/gauteh/8dea955ddb1ed009b48e
                        https://gist.github.com/endolith/1297227
                        https://www.harrisgeospatial.com/docs/INTERPOLATE.html
    """

    shape = np.shape(inArray)
    lin = shape[0]
    col = shape[1]
    # print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    # print(inArray.shape,x,y)

    if kind == "RectBivariateSpline":
        f = RectBivariateSpline(np.arange(0, lin, 1), np.arange(0, col, 1), inArray, kx=3, ky=3, s=0)
        if len(x) > 1 and len(y) > 1:
            res = []
            # res = p_map(_2Dinterp, len(x)*[f], x, y,num_cpus=50)
            for x_, y_ in zip(x, y):
                # print(f(x_, y_).item())
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
        # 'linear in 1d = binlinear in 2D'
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html

        f = interpolate.RegularGridInterpolator(points=(np.arange(0, lin, 1), np.arange(0, col, 1)),
                                                values=inArray,
                                                method=kind,
                                                bounds_error=False,
                                                fill_value=np.nan)

        return f(np.array([x, y]).T)


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
    """

    Args:
        array:
        location:

    Returns:

    """

    loc = int(location)
    tmpArr = array[loc]
    if loc + 1 < len(array):
        res = (array[loc + 1] - tmpArr) * (location - loc) + tmpArr
    else:
        res = (array[-1] - tmpArr) * (location - loc) + tmpArr
    return res
