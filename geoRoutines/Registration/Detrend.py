import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

import geoRoutines.Routine as RT
import geoRoutines.Plotting.Plotting_Routine as RT_Plot
import geoRoutines.Filter.Filtering_Routine as RT_Filter
import geoRoutines.FilesCommandRoutine as FileRT


##### Xarray #####
# import xarray as xr


def detrend_dim(da, dim, deg=1):
    # detrend along a single dimension
    p = da.polyfit(dim=dim, deg=deg)
    fit = xr.polyval(dim, p.polyfit_coefficients)
    return da - fit


def detrend(da, dims, deg=1):
    # detrend along multiple dimensions
    # only valid for linear detrending (deg=1)
    da_detrended = da
    for dim in dims:
        da_detrended = detrend_dim(da_detrended, dim, deg=deg)
    return da_detrended


#####
def IterativeReweighteningLeastSquare(obsArray, jacMat, resDistribution=False, iterations=10, tol=1e-5):
    Z = obsArray
    X = jacMat
    # initiate with standard least square
    A_lsq = np.linalg.lstsq(X, Z, rcond=None)[0]

    # compute absolute value of residuals (fit minus data)
    abs_resid = abs(np.dot(X, A_lsq) - Z)

    # store residuals for output
    iter_mad = np.median(abs_resid)

    # iterate till the fit converges
    A_robust = A_lsq
    for i in range(iterations):
        print(A_robust)
        print('Running iteration', i)

        # compute the scaling factor for the standardization of residuals
        # using the median absolute deviation of the residuals
        # 6.9460 is a tuning constant (4.685/0.6745)
        # equivalent to  http://users.stat.umn.edu/~sandy/courses/8053/handouts/robust.pdf

        if np.median(abs_resid) < tol:
            print('MAR=', np.median(abs_resid))
            print('Convergence reached after', i, 'iteration(s)')

            if i == 0:

                print('This may indicate an issue with the input data as for example constant 0 grids.')
                print('Output will be empty.')
                out_surface = None
                return out_surface

            else:
                break

        abs_res_scale = 6.9460 * np.median(abs_resid)

        # standardize residuals
        w = abs_resid / abs_res_scale

        # compute the robust bisquare weights excluding outliers
        outliers = (w > 1) * 1
        w[outliers.nonzero()] = 0

        good_values = (w != 0) * 1

        # calculate robust weights for 'good' points
        tmp = 1 - np.power(w[good_values.nonzero()], 2)
        w[good_values.nonzero()] = np.power(tmp, 2)

        # get weighted X'es
        XW = np.tile(w, (1, 3)) * X

        a = np.dot(XW.T, X)
        b = np.dot(XW.T, Z)

        # get the least-squares solution to a linear matrix equation
        A_robust = np.linalg.lstsq(a, b)[0]

        # recompute absolute value of residuals (fit minus data)
        abs_resid = abs(np.dot(X, A_robust) - Z)

        iter_mad = np.append(iter_mad, np.median(abs_resid))
        print('MAD =', np.median(abs_resid))

    return A_robust, abs_resid

def Fit2DPlan(inputArray, order=1):
    """
    order1 = val =ax+by+c
    :param inputArray:
    :param order:
    :return:
    """
    ## Generate the output Grid
    shape = inputArray.shape
    xCoord = np.linspace(0, shape[1], shape[1])
    yCoord = np.linspace(0, shape[0], shape[0])
    xx, yy = np.meshgrid(xCoord, yCoord)

    ## Observation Matrix
    obsArray = np.zeros(inputArray.shape)
    obsArray[:] = inputArray

    ## Mask inavalid values and remove this values from the observation matrix
    if not hasattr(inputArray, "mask"):
        inputArray = np.ma.masked_invalid(inputArray)

    if hasattr(inputArray, "mask"):
        mask = inputArray.mask.flatten()
        xx_fl = xx.flatten()[np.invert(mask)]
        yy_fl = yy.flatten()[np.invert(mask)]
        obsArray_fl = obsArray.flatten()[np.invert(mask)]

    else:
        xx_fl = xx.flatten()
        yy_fl = yy.flatten()
        obsArray_fl = obsArray.flatten()

    xx_fl = xx_fl.reshape((len(xx_fl), 1))
    yy_fl = yy_fl.reshape((len(yy_fl), 1))
    obsArray_fl = obsArray_fl.reshape((len(obsArray_fl), 1))

    ### Parameters
    ### functional model f(xi,yj)= axi+byj+c

    ### Construct the Jacobian matrix
    ## jacMat = np.empty(obsArray_fl.shape[0], 3)
    jacMat = np.hstack((xx_fl, yy_fl, np.ones(xx_fl.shape)))
    paramEst, res = IterativeReweighteningLeastSquare(obsArray=obsArray_fl, jacMat=jacMat, resDistribution=False)
    print(paramEst)

    # # Reconstruct unmasked design matrix to get an output values for all input elements (including masked values)

    xx_fl = xx.flatten()
    yy_fl = yy.flatten()
    xx_fl = xx_fl.reshape((len(xx_fl), 1))
    yy_fl = yy_fl.reshape((len(yy_fl), 1))
    A = np.hstack((xx_fl, yy_fl, np.ones(xx_fl.shape)))
    outputGrid = np.reshape(np.dot(A, paramEst), inputArray.shape)

    ## here we wiill take the example of 1st order polynom, we will fit a first order polynom to val=f(x,y)=a1+a2 x+a3 y

    return outputGrid

def Plotting(array1, array2):
    my_cmap = plt.get_cmap('seismic')
    my_cmap = plt.get_cmap('RdBu')
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0, 0])  # row,col
    ax2 = fig.add_subplot(gs[0, 1])  # row,col
    im1 = ax1.imshow(array1, cmap=my_cmap, vmin=-3, vmax=3)
    im2 = ax2.imshow(array2, cmap=my_cmap, vmin=-3, vmax=3)
    RT_Plot.ColorBar(ax=ax2, mapobj=im1)
    plt.show()

    return


def MaskDisplacementMap_(maskFile, map2Mask, band, largeValues=None, visualization=False):
    """

    :param maskFile: mask file path have the same dimension as the map2Mask, composed of 0 and 1. 1 represents the values to be masked
    :param map2Mask: input map path
    :param visualization:
    :return:
    """

    return


def Example(imgPath2Detrend):
    from pathlib import Path
    # path = "//home/cosicorr/0-WorkSpace/3D-Correlation_project/DEMS_correlation/micmac/Med10"
    # imgPath2Detrend = os.path.join(path, "corr_micmac_geoRef_subset.tif")

    img2Detrend_info = RT.GetRasterInfo(imgPath2Detrend, True)
    img2Detrend_array = RT.MultiBandsRaster2Array(imageInfo=img2Detrend_info)

    correctionGrid_EW = Fit2DPlan(inputArray=img2Detrend_array[0])
    img2Detrended_array_corrected_EW = np.subtract(img2Detrend_array[0], correctionGrid_EW)

    # RT_Plot.PlotPlaneFitting(arrayBeforeCorrection=img2Detrend_array[0],
    #                          arrayAfterCorrection=img2Detrended_array_corrected_EW,
    #                          planeArray=correctionGrid_EW, vMin=-3, vMax=3)
    correctionGrid_NS = Fit2DPlan(inputArray=img2Detrend_array[1])
    img2Detrended_array_corrected_NS = np.subtract(img2Detrend_array[1], correctionGrid_NS)

    correctionGrid_dz = Fit2DPlan(inputArray=img2Detrend_array[2])
    img2Detrended_array_corrected_dz = np.subtract(img2Detrend_array[2], correctionGrid_dz)
    RT.WriteRaster(refRasterPath=imgPath2Detrend,
                   newRasterPath=os.path.join(path, Path(imgPath2Detrend).stem + "_detrended.tif"),
                   Listarrays=[img2Detrended_array_corrected_EW, img2Detrended_array_corrected_NS,
                               img2Detrended_array_corrected_dz], numberOfBands=3,
                   descriptions=["E/W detrended", "N/S  detrended", "Dz detrended"])

    return


if __name__ == '__main__':
    # Example()
    # path= "//home/cosicorr/0-WorkSpace/3D-Correlation_project/DEMS_correlation/micmac/med9"
    # imgList= FileRT.ExtractSubfiles(path,[".tif"])
    path = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/Couple_12/3-Otho_PRE_DEM"
    imgList = [
        "/home/cosicorr/0-WorkSpace/3D-Correlation_project/Couple_12/3-Otho_PRE_DEM/3D_displacment_python_Mosaic_parallel.tif"]
    for img_ in imgList:
        Example(img_)
    # path = "E:\OneDrive - California Institute of Technology\\04-Projects\9-3D_approach\Ridgecrest\DetrendingTechnique"
    # img2Detrend = os.path.join(path, "Corr3_3D_NLMF_1p5_21")
    # imgDetrended_withCOSICORR = os.path.join(path, "Corr3_3D_NLMF_1p5_21_detrended")
    # maskPAth = os.path.join(path, "mask_.tif")
    # img2Detrend_info = RT.GetRasterInfo(img2Detrend, printInfo=True)
    # img2Detrended_array = RT.MultiBandsRaster2Array(imageInfo=img2Detrend_info)
    # imgDetrended_withCOSICORR_array = RT.MultiBandsRaster2Array(imageInfo=RT.GetRasterInfo(imgDetrended_withCOSICORR))
    # #
    # # img2Detrended_array_masked = MaskDisplacementMap_(maskFile=maskPAth, map2Mask=img2Detrend,largeValues=None,
    # #                                                   visualization=False,)
    # img2DetrendedArray_masked = RT_Filter.MaskDisplacementMap(maskFile=maskPAth, map2Mask=img2Detrend,
    #                                                              largeValues=None,
    #                                                              visualization=False, bandNumber= 3)
    # # img2DetrendedArray_masked = np.zeros(img2Detrended_array.shape)
    # # for i in range(img2DetrendedArray_masked.shape[0]):
    # #     temp = RT_Filter.MaskDisplacementMap(maskFile=maskPAth, map2Mask=img2Detrend,
    # #                                                                  largeValues=None,
    # #                                                                  visualization=False, bandNumber=i + 1)
    # #
    # #     img2DetrendedArray_masked[i,:,:] = temp[:,:]
    # #     plt.imshow(img2DetrendedArray_masked[i])
    # #     plt.show()
    #
    # # Plotting(array1=img2DetrendedArray_masked[0], array2=imgDetrended_withCOSICORR_array[0])
    # correctionGrid =Fit2DPlan(inputArray=img2DetrendedArray_masked)
    # img2Detrended_array_corrected = np.subtract(img2Detrended_array[2], correctionGrid)
    #
    # # RT_Plot.PlotPlaneFitting(arrayBeforeCorrection=img2Detrended_array[0], arrayAfterCorrection=img2Detrended_array_corrected,
    #                          # planeArray=correctionGrid, vMin=-3, vMax=3)
    # # plt.show()
    # # img2Detrended_array_corrected = Fit2DPlan(inputArray=inputArrayMasked[0])
    # Plotting(array1=img2Detrended_array_corrected,array2 =  imgDetrended_withCOSICORR_array[2])

    # array = np.array ([[1,2,3,5],[1,2,6,8]])
    # print(array)
    # mask = np.array ([[True,True,False,False],[True,True,False,False]])
    # maskBin = mask*1
    # mask_array = maskBin==1
    # print(mask_array)
    # array_Masked = np.ma.masked_array(array,mask =mask_array)
    # print(array_Masked)
