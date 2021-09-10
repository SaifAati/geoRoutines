import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal
from scipy.stats import norm

import geospatialroutine.FilesCommandRoutine as fileRT
import geospatialroutine.Filter.Filtering_Routine as filterRT
import geospatialroutine.Routine as RT
import geospatialroutine.Plotting.Plotting_Routine as pltRT
import geospatialroutine.Correlation_misc as corrmiscRT


def GetStatisitic(inputArray):
    statDict = {}
    val_mean = np.ma.mean(inputArray)
    val_std = np.ma.std(inputArray)
    statDict['mean'] = val_mean
    statDict['std'] = val_std
    statDict['rmse'] = np.sqrt(np.ma.mean(np.square(inputArray)))
    statDict['max'] = np.ma.max(inputArray)
    statDict['min'] = np.ma.min(inputArray)
    if isinstance(inputArray, np.ma.MaskedArray):
        statDict['n_samples'] = inputArray.count()
    else:
        statDict['n_samples'] = inputArray.size

    return statDict


def FitPlaneIter(inputDispArray, tune=4.685,
                 iterations=10, tol=1e-5):
    """
    Fit a plane a 2D array using robust least square (L1-norm)

    :param inputDispArray: 2D array to be fitted, if masked array the mask is not taken into account
    :param tune: Tuckey tuning constant (default 4.685, currently set to common default)
    :param iterations: Maximum number of iterations (default 10)
    :param tol: Stopping criteria. The fiting will stop before iter is reached if the residuals are below
                     this threshold.
    :return: fitted plane as 2D numpy array, MAD of the fits during each iteration
    """
    print("\n--- Start Iterative plane fitting")
    # Create a coordinate matrix
    x_size = inputDispArray.shape[1]
    y_size = inputDispArray.shape[0]
    nx = np.linspace(-1, 1, x_size)
    ny = np.linspace(-1, 1, y_size)
    x, y = np.meshgrid(nx, ny)

    # Construct design matrix
    x_fl = x.flatten()
    y_fl = y.flatten()
    z_ones = np.ones([x.size, 1])

    # Transfer mask
    if hasattr(inputDispArray, 'mask'):
        mask = inputDispArray.mask.flatten()
        x_fl = x_fl[np.invert(mask)]
        y_fl = y_fl[np.invert(mask)]
        z_ones = z_ones[np.invert(mask)]

    X = np.hstack((np.reshape(x_fl, ([len(x_fl), 1])), np.reshape(y_fl, ([len(y_fl), 1])), z_ones))

    # Construct response matrix
    Z = np.zeros(inputDispArray.shape)
    Z[:] = inputDispArray
    Z_fl = Z.flatten()
    # Transfer mask
    if hasattr(inputDispArray, 'mask'):
        Z_fl = Z_fl[np.invert(mask)]
    Z = np.reshape(Z_fl, ([len(Z_fl), 1]))

    # Initiate with standard least square
    A_lsq = np.linalg.lstsq(X, Z, rcond=None)[0]

    # compute absolute value of residuals (fit minus data)
    abs_resid = abs(np.dot(X, A_lsq) - Z)

    # store residuals for output
    iter_mad = np.median(abs_resid)

    # iterate till the fit converges
    for i in range(iterations):
        # print('Running iteration', i)
        # compute the scaling factor for the standardization of residuals
        # using the median absolute deviation of the residuals
        # 6.9460 is a tuning constant (4.685/0.6745)
        # equivalent to  http://users.stat.umn.edu/~sandy/courses/8053/handouts/robust.pdf

        if np.median(abs_resid) < tol:
            # print('MAR median=', np.median(abs_resid))
            print('Convergence reached after', i, 'iteration(s)')
            if i == 0:
                print('This may indicate an issue with the input data as for example constant 0 grids.')
                print('Output will be empty.')
                outputDispArr = None
                return outputDispArr
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
        A_robust = np.linalg.lstsq(a, b, rcond=None)[0]

        # recompute absolute value of residuals (fit minus data)
        abs_resid = abs(np.dot(X, A_robust) - Z)

        iter_mad = np.append(iter_mad, np.median(abs_resid))
        # print('MAD =', np.median(abs_resid))

    # get values of the approximated plane

    # reconstruct unmasked design matrix to get an output values for all input elements (including masked values)
    x_fl = x.flatten()
    y_fl = y.flatten()
    z_ones = np.ones([x.size, 1])
    X = np.hstack((np.reshape(x_fl, ([len(x_fl), 1])), np.reshape(y_fl, ([len(y_fl), 1])), z_ones))

    out_values = np.dot(X, A_robust)

    # # bring it into the format of the input matrix

    # if hasattr(inputDispArray, 'mask'):
    #     outputDispArr = np.zeros(inputDispArray.shape)
    #     outputDispArr[np.invert(inputDispArray.mask)] = out_values.flatten()
    # else:
    outputDispArr = np.reshape(out_values, inputDispArray.shape)
    print("--- End Iterative plane fitting: number of iteration", i)
    return outputDispArr, iter_mad.tolist(),A_robust


def OffsetCorrectionPixelBased_SA(dispMap,
                                  bandNumber,
                                  refPath=None,
                                  snrMap=None,
                                  snrBandNb=None,
                                  snrThreshold=0.333,
                                  maskLargeValues=None,
                                  maskZeros=False,
                                  # Mask large value more than this value and less than -this value
                                  iterations=3,
                                  zoom=1,

                                  showDiagnose=True,
                                  saveDiagnose=False,
                                  diagnoseSavingPath="",
                                  imgIndex=None,

                                  max_displacement=-30,
                                  xlim=[-1, 1], ylim=[0, 3]):
    """
    Approximate the displacement field with a plane and return the corrected grid, a dictionary with the log of
    the corrections and the correction grid that can be used to transform the slave image.
    :param displacement_field: full path to the displacement fields that should be corrected
    :param master_slave: master slave dictionary
    :param ref_path: full path to a reference image, if the displacement field has no georeference the reference
                                of the reference image will be used
    :param snr_filed:  full path to the raster holding the correlation coefficients
    :param snr_threshold:   threshold on the correlation coefficient, pixels below this value will be ignored during the
                                    correction
    :param iterations: number of iterations for iteratively reweighted plane fitting
    :param zoom: optional downsampling, mainly for reducing memory size and processing time
    :param diagnose: if True plotting is activated to illustrate the corrections
    :param maskLargeValues: Bool if we want to mask large values
    :param max_displacement: the values of values to be masked, all displacement values whose absolute magnitude exceeds
                            this value will be masked out for the correction
    :param xlim: x axis limits of the histogram plots
    :param ylim: y axis limits of the histogram plots
    :return: the corrected grid, a dictionary with the log of the corrections and the correction grid
    """

    print("***** OffsetCorrection for image :", dispMap, "*****")
    # Allocate output dictionary
    out_dict = {}
    outDict = {}

    # Displacement raster info
    displacementInfo = RT.GetRasterInfo(inputRaster=dispMap)
    # This allows GDAL to throw Python Exceptions
    gdal.UseExceptions()

    GSD = displacementInfo.get("Resolution")[0]
    bands = displacementInfo.get("NbBands")
    print("GSD:", GSD, ", NbBands:", bands)

    displacementRaster = displacementInfo.get("Raster")
    displBand = displacementRaster.GetRasterBand(bandNumber)
    inputDisp_arrayPix = np.array(displBand.ReadAsArray()) / GSD  ## pixel based
    inputDisp_arrayPix_origin = np.array(inputDisp_arrayPix)

    # Get the  georeference of the displacement field
    if displacementInfo.get("MapInfo")[0] == False:
        print('Input displacement has no georeference. Trying reference image instead...')
        try:
            print("Ref Image:", refPath, " will be used as map reference for the displacement map")
            refRaster = gdal.Open(refPath)
            geotransform = refRaster.GetGeoTransform()
        except:
            raise IOError('Cannot read georeference from ' + refPath)
    else:
        geotransform = displacementInfo.get("MapInfo")[1]

    ## Mask values with SNR
    if snrMap is not None:
        # Import correlation coefficient
        snrBand = displacementRaster.GetRasterBand(snrBandNb)
        inputSNR_array = np.ma.array(snrBand.ReadAsArray())
        inputDisp_arrayPix = Filter.MaskDispWithSNR(snrArray=inputSNR_array, dispArray=inputDisp_arrayPix,
                                                    snrThreshold=snrThreshold, visualization=showDiagnose)
    ## Mask large displacements
    if maskLargeValues:
        print("maskLargeValues")
        inputDisp_arrayPix = np.ma.array(Filter.MaskLargeValues(inputArray=inputDisp_arrayPix * GSD,
                                                                maxDisplacement=maskLargeValues,
                                                                visualization=showDiagnose)) / GSD

    if maskZeros == True:
        inputDisp_arrayPix = Filter.MaskZeros(inputArray=inputDisp_arrayPix)

    inputDisp_array_before_correction = inputDisp_arrayPix  ## store raster before offset correction
    ### Mask nan values
    inputDisp_arrayPix = np.ma.masked_invalid(inputDisp_arrayPix)
    # ## Get statistical infos on the input before offset correction

    outDict["Input"] = {}
    outDict["Input"]["name"] = os.path.basename(dispMap)
    outDict["Input"]["bandNb"] = bandNumber
    outDict["Input"].update(GetStatisitic(inputArray=inputDisp_arrayPix * GSD))

    print("--- Statistic about input displacement map=\n", outDict)

    # inputDisp_arrayGaussianFilterPix, mean_gauss, std_gauss, rmse_gauss = Filter.IterativeGaussianFit(
    #     sample=inputDisp_arrayPix,
    #     showDiagnose=showDiagnose,
    #     saveDiagnose=saveDiagnose,
    #     diagnoseSavePath=diagnoseSavingPath,
    #     imgIndex=imgIndex,
    #     title='Residuals before corrections',
    #     xlim=xlim, ylim=ylim)
    #
    # out_dict['input']['rmse_gauss'] = rmse_gauss
    # out_dict['input']['std_gauss'] = std_gauss
    # out_dict['input']['mean_gauss'] = mean_gauss
    # print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)
    print("\n--- Fit a plane ")
    # # Fit plane
    correctionSurface, iterMad = FitPlaneIter(inputDispArray=inputDisp_arrayPix, iterations=iterations)

    # correction with plane == deramping
    inputDisp_arrayPix -= correctionSurface

    # cmap = 'seismic'
    # # cmap = "rainbow"
    # # cmap = 'gist_earth'
    if showDiagnose or saveDiagnose:
        Plot.PlotPlaneFitting(arrayBeforeCorrection=inputDisp_array_before_correction * GSD,
                              arrayAfterCorrection=inputDisp_arrayPix * GSD, planeArray=correctionSurface * GSD,
                              show=showDiagnose,
                              save=saveDiagnose)

    # Get infos after plane fitting
    outDict['PlaneCorrected'] = {}
    outDict['PlaneCorrected']['iter_mad'] = iterMad
    outDict['PlaneCorrected'].update(GetStatisitic(inputArray=inputDisp_arrayPix * GSD))
    print("--- Statistic about displacement map after Gaussian filter=\n", outDict)

    # inputDisp_arrayGaussianFilterPix, mean_gauss, std_gauss, rmse_gauss = Filter.IterativeGaussianFit(
    #     sample=inputDisp_arrayPix,
    #     showDiagnose=showDiagnose,
    #     title='Residuals after Deramping',
    #     xlim=xlim, ylim=ylim)
    #
    # out_dict['plane_corrected']['rmse_gauss'] = rmse_gauss
    # out_dict['plane_corrected']['std_gauss'] = std_gauss
    # out_dict['plane_corrected']['mean_gauss'] = mean_gauss
    #
    # print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)

    return inputDisp_arrayPix * GSD


def OffsetCorrectionPixelBased(displacement_field,
                               ref_path,
                               bandNumber=1,
                               snr_filed=None,
                               snr_threshold=0.333,
                               iterations=3,
                               zoom=1,
                               showDiagnose=True,
                               saveDiagnose=False,
                               diagnoseSavingPath="",
                               imgIndex=None,
                               maskLargeValues=False,
                               noData=-32767,
                               maskZeros=False,
                               max_displacement=-30,
                               xlim=[-1, 1], ylim=[0, 3]):
    """
    Approximate the displacement field with a plane and return the corrected grid, a dictionary with the log of
    the corrections and the correction grid that can be used to transform the slave image.
    :param displacement_field: full path to the displacement fields that should be corrected
    :param master_slave: master slave dictionary
    :param ref_path: full path to a reference image, if the displacement field has no georeference the reference
                                of the reference image will be used
    :param snr_filed:  full path to the raster holding the correlation coefficients
    :param snr_threshold:   threshold on the correlation coefficient, pixels below this value will be ignored during the
                                    correction
    :param iterations: number of iterations for iteratively reweighted plane fitting
    :param zoom: optional downsampling, mainly for reducing memory size and processing time
    :param diagnose: if True plotting is activated to illustrate the corrections
    :param maskLargeValues: Bool if we want to mask large values
    :param max_displacement: the values of values to be masked, all displacement values whose absolute magnitude exceeds
                            this value will be masked out for the correction
    :param xlim: x axis limits of the histogram plots
    :param ylim: y axis limits of the histogram plots
    :return: the corrected grid, a dictionary with the log of the corrections and the correction grid
    """

    print("***** OffsetCorrection for image :", displacement_field, "*****")
    # Allocate output dictionary
    out_dict = {}

    # Displacement raster info
    displacementInfo = RT.GetRasterInfo(inputRaster=displacement_field)
    # This allows GDAL to throw Python Exceptions
    gdal.UseExceptions()

    GSD = displacementInfo.get("Resolution")[0]
    bands = displacementInfo.get("NbBands")
    print("=== > GSD, bands", GSD, bands)

    displacementRaster = displacementInfo.get("Raster")
    displBand = displacementRaster.GetRasterBand(bandNumber)
    inputDisp_arrayPix = np.array(displBand.ReadAsArray()) / GSD  ## pixel based
    inputDisp_array = np.array(displBand.ReadAsArray())  ## in meter

    # Get the  georeference of the displacement field
    if displacementInfo.get("MapInfo")[0] == False:
        print('Input field has no georeference. Trying reference image instead...')
        try:
            refRaster = gdal.Open(ref_path)
            geotransform = refRaster.GetGeoTransform()
        except:
            raise IOError('Cannot read georeference from ' + ref_path)
    else:
        geotransform = displacementInfo.get("MapInfo")[1]

    ## Mask values with SNR
    if snr_filed is not None:
        # Import correlation coefficient
        snrBand = displacementRaster.GetRasterBand(2)
        inputSNR_array = np.array(snrBand.ReadAsArray())

        inputDisp_arrayPix = filterRT.MaskDispWithSNR(snrField=inputSNR_array, dispArray=inputDisp_arrayPix,
                                                      snrThreshold=snr_threshold)
        inputDisp_array = filterRT.MaskDispWithSNR(snrField=inputSNR_array, dispArray=inputDisp_array,
                                                   snrThreshold=snr_threshold)

    ## Mask large displacements

    inputDisp_arrayPix = np.ma.masked_where(inputDisp_arrayPix == noData, inputDisp_arrayPix)
    print("+++++max", np.ma.max(inputDisp_arrayPix))
    # if maskLargeValues == True:
    #     inputDisp_arrayPix = filterRT.MaskLargeValues(inputArray=inputDisp_arrayPix, maxDisplacement=max_displacement)
    #     inputDisp_array = filterRT.MaskLargeValues(inputArray=inputDisp_array, maxDisplacement=max_displacement)
    # if maskZeros == True:
    #     inputDisp_arrayPix = filterRT.MaskZeros(inputArray=inputDisp_arrayPix)
    #     inputDisp_array = filterRT.MaskZeros(inputArray=inputDisp_array)

    inputDisp_array_before_correction = np.copy(inputDisp_arrayPix)
    # print(inputDisp_array_before_correction)
    ### Mask nan values
    # if isinstance(inputDisp_arrayPix, np.ma.MaskedArray):
    inputDisp_arrayPix = np.ma.masked_invalid(inputDisp_arrayPix)
    # print(inputDisp_array_before_correction)
    ## Get statistical infos on the input
    val_mean = np.ma.mean(inputDisp_arrayPix)
    val_std = np.ma.std(inputDisp_arrayPix)

    val_min = val_mean - val_std
    val_max = val_mean + val_std
    print("===========================================")
    print("val_mean=%f,val_std=%f" % (val_mean, val_std))
    print("val_min=%f,val_max=%f" % (val_min, val_max))
    print("===========================================")
    out_dict['input'] = {}
    out_dict['input']['name'] = displacement_field
    out_dict['input']['mean'] = val_mean
    out_dict['input']['std'] = val_std
    out_dict['input']['rmse'] = np.sqrt(np.ma.mean(np.square(inputDisp_arrayPix)))
    out_dict['input']['max'] = np.ma.max(inputDisp_arrayPix)
    out_dict['input']['min'] = np.ma.min(inputDisp_arrayPix)

    if isinstance(inputDisp_arrayPix, np.ma.MaskedArray):
        out_dict['input']['n_samples'] = inputDisp_arrayPix.count()
    else:
        out_dict['input']['n_samples'] = inputDisp_arrayPix.size

    print("=== > Input data stat:\n", out_dict)
    inputDisp_arrayGaussianFilterPix, mean_gauss, std_gauss, rmse_gauss = filterRT.IterativeGaussianFit(
        sample=inputDisp_arrayPix,
        showDiagnose=showDiagnose,
        saveDiagnose=saveDiagnose,
        diagnoseSavePath=diagnoseSavingPath,
        imgIndex=imgIndex,
        title='Residuals before corrections',
        xlim=xlim, ylim=ylim)

    out_dict['input']['rmse_gauss'] = rmse_gauss
    out_dict['input']['std_gauss'] = std_gauss
    out_dict['input']['mean_gauss'] = mean_gauss
    print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)

    print("\n--- Fit a plan ")

    # Fit plane
    correctionSurface, iter_mad = FitPlaneIter(inputDispArray=inputDisp_arrayPix, iterations=iterations)

    # correction with plane == deramping
    inputDisp_arrayPix -= correctionSurface

    cmap = 'seismic'
    # cmap = "rainbow"
    # cmap = 'gist_earth'
    if showDiagnose or saveDiagnose:
        print("Plotiing ....")
        # plot
        fig, axes = plt.subplots(ncols=3, figsize=(10, 3))
        vis_max = max(abs(val_min), abs(val_max))
        vis_min = -vis_max
        # vis_min =-0.5
        # vis_max = 0.5
        axes[0].imshow(inputDisp_array_before_correction, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[0].set_title('Input Image', fontsize=12)
        axes[0].tick_params(labelsize=12)

        # plot
        axes[1].imshow(correctionSurface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[1].set_title('Fitted plane', fontsize=12)
        axes[1].tick_params(labelsize=12)

        im = axes[2].imshow(inputDisp_arrayPix, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[2].set_title('After plane correction', fontsize=12)
        axes[2].tick_params(labelsize=12)
        fig.subplots_adjust(right=0.91)

        divider = make_axes_locatable(axes[2])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=12)
        axes[0].axis("off")
        axes[1].axis("off")
        axes[2].axis("off")
        if showDiagnose:
            # plt.show()
            plt.draw()
            plt.pause(1)
        if saveDiagnose:
            fig.savefig(diagnoseSavingPath + "img_" + str(imgIndex) + ".png")
        plt.close()

    # Get infos after plane fitting
    out_dict['plane_corrected'] = {}
    out_dict['plane_corrected']['iter_mad'] = iter_mad
    out_dict['plane_corrected']['mean'] = np.ma.mean(inputDisp_arrayPix)
    out_dict['plane_corrected']['std'] = np.ma.std(inputDisp_arrayPix)
    out_dict['plane_corrected']['rmse'] = np.sqrt(np.ma.mean(np.square(inputDisp_arrayPix)))
    out_dict['plane_corrected']['max'] = np.ma.max(inputDisp_arrayPix)
    out_dict['plane_corrected']['min'] = np.ma.min(inputDisp_arrayPix)

    inputDisp_arrayGaussianFilterPix, mean_gauss, std_gauss, rmse_gauss = filterRT.IterativeGaussianFit(
        sample=inputDisp_arrayPix,
        showDiagnose=showDiagnose,
        title='Residuals after Deramping',
        xlim=xlim, ylim=ylim)

    out_dict['plane_corrected']['rmse_gauss'] = rmse_gauss
    out_dict['plane_corrected']['std_gauss'] = std_gauss
    out_dict['plane_corrected']['mean_gauss'] = mean_gauss

    print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)

    return inputDisp_arrayPix * GSD


def OffsetCorrectionGSDBased(displacement_field,
                             master_slave,
                             ref_path,

                             snr_filed=None,
                             snr_threshold=0.333,

                             iterations=3,
                             zoom=1,
                             diagnose=True,
                             maskLargeValues=False,
                             max_displacement=-30,
                             xlim=[-15, 15], ylim=[0, 3]):
    """
    Approximate the displacement field with a plane and return the corrected grid, a dictionary with the log of
    the corrections and the correction grid that can be used to transform the slave image.
    :param displacement_field: full path to the displacement fields that should be corrected
    :param master_slave: master slave dictionary
    :param ref_path: full path to a reference image, if the displacement field has no georeference the reference
                                of the reference image will be used
    :param snr_filed:  full path to the raster holding the correlation coefficients
    :param snr_threshold:   threshold on the correlation coefficient, pixels below this value will be ignored during the
                                    correction
    :param iterations: number of iterations for iteratively reweighted plane fitting
    :param zoom: optional downsampling, mainly for reducing memory size and processing time
    :param diagnose: if True plotting is activated to illustrate the corrections
    :param maskLargeValues: Bool if we want to mask large values
    :param max_displacement: the values of values to be masked, all displacement values whose absolute magnitude exceeds
                            this value will be masked out for the correction
    :param xlim: x axis limits of the histogram plots
    :param ylim: y axis limits of the histogram plots
    :return: the corrected grid, a dictionary with the log of the corrections and the correction grid
    """

    print("***** OffsetCorrection for image :", displacement_field, "*****")
    # Allocate output dictionary
    out_dict = {}

    # Displacement raster info
    displacementInfo = RT.GetRasterInfo(inputRaster=displacement_field)
    # This allows GDAL to throw Python Exceptions
    gdal.UseExceptions()

    GSD = displacementInfo.get("Resolution")[0]
    bands = displacementInfo.get("NbBands")
    print("GSD, bands", GSD, bands)

    displacementRaster = displacementInfo.get("Raster")
    displBand = displacementRaster.GetRasterBand(1)
    # inputDisp_arrayPix = np.array(displBand.ReadAsArray()) / GSD  ## pixel based
    inputDisp_array = np.array(displBand.ReadAsArray())  ## in meter

    # Get the  georeference of the displacement field
    if displacementInfo.get("MapInfo")[0] == False:
        print('Input field has no georeference. Trying reference image instead...')
        try:
            refRaster = gdal.Open(ref_path)
            geotransform = refRaster.GetGeoTransform()
        except:
            raise IOError('Cannot read georeference from ' + ref_path)
    else:
        geotransform = displacementInfo.get("MapInfo")[1]

    ## Mask values with SNR
    if snr_filed is not None:
        # Import correlation coefficient
        snrBand = displacementRaster.GetRasterBand(2)
        inputSNR_array = np.array(snrBand.ReadAsArray())

        inputDisp_array = filterRT.MaskDispWithSNR(snrField=inputSNR_array, dispArray=inputDisp_array,
                                                   snrThreshold=snr_threshold)

    ## Mask large displacements
    if maskLargeValues == True:

        inputDisp_array = filterRT.MaskLargeValues(inputArray=inputDisp_array, maxDisplacement=max_displacement)

        ### Mask nan values
        if isinstance(inputDisp_array, np.ma.MaskedArray):
            inputDisp_array = np.ma.masked_invalid(inputDisp_array)

    ## Get statistical infos on the input
    val_mean = np.ma.mean(inputDisp_array)
    val_std = np.ma.std(inputDisp_array)

    val_min = val_mean - val_std
    val_max = val_mean + val_std
    print("val_min=%f,val_max=%f" % (val_min, val_max))
    out_dict['input'] = {}
    out_dict['input']['name'] = displacement_field
    out_dict['input']['mean'] = val_mean
    out_dict['input']['std'] = val_std
    out_dict['input']['rmse'] = np.sqrt(np.ma.mean(np.square(inputDisp_array)))
    out_dict['input']['max'] = np.ma.max(inputDisp_array)
    out_dict['input']['min'] = np.ma.min(inputDisp_array)

    if isinstance(inputDisp_array, np.ma.MaskedArray):
        out_dict['input']['n_samples'] = inputDisp_array.count()
    else:
        out_dict['input']['n_samples'] = inputDisp_array.size

    print("--- Statistic about input displacement map=\n", out_dict)
    inputDisp_arrayGaussianFilterPix, mean_gauss, std_gauss, rmse_gauss = filterRT.IterativeGaussianFit(
        sample=inputDisp_array,
        diagnose=diagnose,
        title='Residuals before corrections',
        xlim=xlim, ylim=ylim)

    out_dict['input']['rmse_gauss'] = rmse_gauss
    out_dict['input']['std_gauss'] = std_gauss
    out_dict['input']['mean_gauss'] = mean_gauss
    print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)

    print("\n--- Fit a plan ")
    # Fit plane
    correctionSurface, iter_mad = FitPlaneIter(inputDispArray=inputDisp_array, iterations=iterations)

    # correction with plane == deramping
    inputDisp_array -= correctionSurface

    cmap = 'seismic'
    if diagnose:
        print("Plotiing ....")
        # plot
        fig, axes = plt.subplots(ncols=3, figsize=(10, 10))
        vis_max = max(abs(val_min), abs(val_max))
        vis_min = -vis_max
        vis_min = 0
        vis_max = 0.25

        axes[0].imshow(inputDisp_array, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[0].set_title('Input', fontsize=16)
        axes[0].tick_params(labelsize=14)

        # plot
        axes[1].imshow(correctionSurface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[1].set_title('Fitted plane', fontsize=16)
        axes[1].tick_params(labelsize=14)

        im = axes[2].imshow(inputDisp_array, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[2].set_title('After plane correction', fontsize=16)
        axes[2].tick_params(labelsize=14)
        fig.subplots_adjust(right=0.91)

        cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.3])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.ax.tick_params(labelsize=16)
        axes[0].axis("off")
        axes[1].axis("off")
        axes[2].axis("off")
        plt.draw()
        plt.pause(2)
        plt.close()

    # Get infos after plane fitting
    out_dict['plane_corrected'] = {}
    out_dict['plane_corrected']['iter_mad'] = iter_mad
    out_dict['plane_corrected']['mean'] = np.ma.mean(inputDisp_array)
    out_dict['plane_corrected']['std'] = np.ma.std(inputDisp_array)
    out_dict['plane_corrected']['rmse'] = np.sqrt(np.ma.mean(np.square(inputDisp_array)))
    out_dict['plane_corrected']['max'] = np.ma.max(inputDisp_array)
    out_dict['plane_corrected']['min'] = np.ma.min(inputDisp_array)

    inputDisp_arrayGaussianFilterPix, mean_gauss, std_gauss, rmse_gauss = filterRT.IterativeGaussianFit(
        sample=inputDisp_array,
        diagnose=diagnose,
        title='Residuals after Deramping',
        xlim=xlim, ylim=ylim)

    out_dict['plane_corrected']['rmse_gauss'] = rmse_gauss
    out_dict['plane_corrected']['std_gauss'] = std_gauss
    out_dict['plane_corrected']['mean_gauss'] = mean_gauss

    print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)

    return inputDisp_array


def MainCorrOffset_Batch(inputPath, savingPath="", saveGrid=False, gridSavingPath=""):
    imgs = FileRT.FilesInDirectory(path=inputPath)
    workspacePath = Path(inputPath).parents[1]
    if not savingPath:
        folders = FileRT.FilesInDirectory(path=workspacePath)
        print(folders)
        if "OffsetCorrection" in folders:

            tempPath = os.path.join(workspacePath, "OffsetCorrection")
        else:
            tempPath = FileRT.CreateDirectory(directoryPath=workspacePath, folderName="OffsetCorrection")

        savingPath = FileRT.CreateDirectory(directoryPath=tempPath,
                                            folderName=os.path.basename(os.path.normpath(inputPath)) + "_OffsetCorr")

    if saveGrid:
        if not gridSavingPath:
            gridSavingPath = FileRT.CreateDirectory(directoryPath=savingPath, folderName="CorrectionGrids")

    for index, img_ in enumerate(imgs):
        if ".tif" in img_:
            print("\n ------------ image number :", index + 1, "---------------\n")

            imgPath = inputPath + img_

            output = OffsetCorrectionPixelBased(displacement_field=imgPath,
                                                ref_path=imgPath, snr_filed=None, saveDiagnose=saveGrid,
                                                diagnoseSavingPath=gridSavingPath, imgIndex=index + 1,
                                                maskLargeValues=False, maskZeros=True, showDiagnose=False)
            rasterInfo = RT.GetRasterInfo(imgPath)
            metaData = rasterInfo["MetaData"]
            RT.WriteRaster(refRasterPath=imgPath, newRasterPath=savingPath + img_, Listarrays=[output],
                           numberOfBands=1, metaData=metaData)

    return


def animationFunction(path1, path2, animationRecording=False, figTitle="fig", title1="", title2="",
                      animationTitle="", vMin=-4, vMax=4, transparency=False):
    import matplotlib.animation as animation
    import dateutil.parser as dparser

    paths = [path1, path2]
    files = []
    for path_ in paths:
        tmp = []

        fileTmp = FileRT.FilesInDirectory(path_)
        fileTmp.sort()
        for file_ in fileTmp:
            if ".tif" in file_ and ".xml" not in file_:
                tmp.append(file_)
        files.append(tmp)
    print(files[0])
    print(files[1])
    fileList1 = files[0]  # list of raster to be displayed
    dateList = []  # list of date
    satType = []  # list of satellite type
    pathrowList = []
    diffDays = []
    # Fetching Date and satellite type
    for file_ in fileList1:
        # Fetching Date and satellite type
        splitList = file_.split("_")
        dateList.append(splitList[1] + "#" + splitList[4])
        satType.append(splitList[5])
        masterDate = dparser.parse(splitList[1], fuzzy=True)
        slaveDate = dparser.parse(splitList[4], fuzzy=True)
        daysDiff = abs(masterDate - slaveDate)
        diffDays.append(daysDiff.days)
        if "_PR_" in file_:
            pathrowList.append(splitList[7])
        else:
            pathrowList.append("")

    fig, axes = plt.subplots(2, 2)
    fig.canvas.set_window_title(figTitle)

    # subsetImg = "/home/cosi/2-Data/4-Shishper/DownSamplingImg/SL/subsetofsubset"
    subsetImg = "/home/cosi/2-Data/4-Shishper/DownSamplingImg/SL/DownSimpled_SL_V2.tif"
    subsetImgAsRaster = RT.ImageAsArray(subsetImg)
    subsetImgAsRasterNorm = RT.NormalizeImage(subsetImgAsRaster)
    ax1 = axes[0, 0]
    ax2 = axes[0, 1]
    ax3 = axes[1, 0]
    ax4 = axes[1, 1]

    ax1.set_title(title1, fontsize=10)
    ax2.set_title(title2, fontsize=10)
    ax1.imshow(subsetImgAsRasterNorm, cmap="gray")
    ax2.imshow(subsetImgAsRasterNorm, cmap="gray")
    ax1.axis("off")
    ax2.axis("off")
    ax3.get_yaxis().set_visible(False)
    ax4.get_yaxis().set_visible(False)

    # my_cmap = plt.get_cmap('gist_ncar')
    # my_cmap = plt.get_cmap('jet')
    # my_cmap = plt.get_cmap('rainbow')
    # my_cmap = plt.get_cmap('Reds')
    my_cmap = plt.get_cmap('seismic')

    ## Plotting

    # Animation initialization

    imageAsArray1 = np.array(RT.ImageAsArray(path1 + files[0][0]))
    imageAsArray2 = RT.ImageAsArray(path2 + files[1][0])
    if transparency == True:
        imageAsArray1 = np.ma.masked_inside(imageAsArray1, 0.00001, -0.00001)
        imageAsArray2 = np.ma.masked_inside(imageAsArray2, 0.00001, -0.00001)
    vmin = vMin
    vmax = vMax

    im1 = ax1.pcolormesh(imageAsArray1, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)
    im2 = ax2.pcolormesh(imageAsArray2, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)

    txtTitle = dateList[0] + " (" + satType[0] + pathrowList[0] + ")" + "  diff=" + str(diffDays[0])
    # txt = ax1.set_title(100, -40, txtTitle, color="red", fontsize=10)
    # txt = ax1.set_title( txtTitle,fontdict={'fontsize': 20, 'fontweight': 'medium'}, loc= "center", pad= 2.52)
    titlte = fig.suptitle(txtTitle, fontsize=16)

    imageAsArray1_vector = imageAsArray1[~np.isnan(imageAsArray1)].flatten()
    imageAsArray2_vector = imageAsArray2[~np.isnan(imageAsArray2)].flatten()

    (mu, sigma) = norm.fit(imageAsArray1_vector)
    print(mu, sigma)
    RMSE1 = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(imageAsArray1_vector))))
    lower_bin = mu - 3 * sigma
    upper_bin = mu + 3 * sigma
    hist1, bins1 = np.histogram(imageAsArray1_vector, range=[lower_bin, upper_bin], density=False, bins=100)
    bars1 = ax3.bar((bins1[:-1] + bins1[1:]) / 2, hist1, align='center', width=(bins1[1] - bins1[0]))

    (mu, sigma) = norm.fit(imageAsArray2_vector)
    print("2", mu, sigma)
    RMSE2 = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(imageAsArray2_vector))))
    lower_bin = mu - 3 * sigma
    upper_bin = mu + 3 * sigma
    hist2, bins2 = np.histogram(imageAsArray2_vector, range=[lower_bin, upper_bin], density=False, bins=100)
    bars2 = ax4.bar((bins2[:-1] + bins2[1:]) / 2, hist2, align='center', width=(bins2[1] - bins2[0]))

    ax3.set_title('Gaussian fit')

    text1 = ax3.text(0.9, 0.9, 'RMSE=' + str(RMSE1),
                     horizontalalignment='right',
                     verticalalignment='top',
                     transform=ax3.transAxes)

    text2 = ax4.text(0.9, 0.9, 'RMSE=' + str(RMSE2),
                     horizontalalignment='right',
                     verticalalignment='top',
                     transform=ax4.transAxes)

    def updatefig(i):
        ax3.clear()
        ax4.clear()
        if i >= len(files[0]):
            # (ax1,ax2,ax3,ax4).clear()
            print("i>=len(newX)")
            i = len(files[0]) - 1
            print("frame i=", i, "   ", files[0][i])

            imageAsArray1 = np.array(RT.ImageAsArray(path1 + files[0][i]))
            imageAsArray2 = RT.ImageAsArray(path2 + files[1][i])
            if transparency == True:
                imageAsArray1 = np.ma.masked_inside(imageAsArray1, 0.00001, -0.00001)
                imageAsArray2 = np.ma.masked_inside(imageAsArray2, 0.00001, -0.00001)

            im1 = ax1.pcolormesh(imageAsArray1, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)
            im2 = ax2.pcolormesh(imageAsArray2, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)

            txtTitle = dateList[i] + " (" + satType[i] + pathrowList[i] + ")" + "  diff=" + str(diffDays[i])
            # txt = ax1.set_title(100, -40, txtTitle, color="red", fontsize=10)
            # txt = ax1.set_title(txtTitle)
            titlte = fig.suptitle(txtTitle, fontsize=16, color="red")

            imageAsArray1_vector = imageAsArray1[~np.isnan(imageAsArray1)].flatten()
            imageAsArray2_vector = imageAsArray2[~np.isnan(imageAsArray2)].flatten()
            (mu, sigma) = norm.fit(imageAsArray1_vector)
            print("1", mu, sigma)
            RMSE1 = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(imageAsArray1_vector))))
            lower_bin = mu - 3 * sigma
            upper_bin = mu + 3 * sigma
            hist1, bins1 = np.histogram(imageAsArray1_vector, range=[lower_bin, upper_bin], density=False, bins=100)
            bars1 = ax3.bar((bins1[:-1] + bins1[1:]) / 2, hist1, align='center', width=(bins1[1] - bins1[0]))

            (mu, sigma) = norm.fit(imageAsArray2_vector)
            print("2", mu, sigma)
            RMSE2 = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(imageAsArray2_vector))))
            lower_bin = mu - 3 * sigma
            upper_bin = mu + 3 * sigma
            hist2, bins2 = np.histogram(imageAsArray2_vector, range=[lower_bin, upper_bin], density=False, bins=100)
            bars2 = ax4.bar((bins2[:-1] + bins2[1:]) / 2, hist2, align='center', width=(bins2[1] - bins2[0]))

            text1 = ax3.text(0.9, 0.9, 'RMSE=' + str(RMSE1),
                             horizontalalignment='right',
                             verticalalignment='top',
                             transform=ax3.transAxes)

            text2 = ax4.text(0.9, 0.9, 'RMSE=' + str(RMSE2),
                             horizontalalignment='right',
                             verticalalignment='top',
                             transform=ax4.transAxes)
            return im1, im2, titlte, bars1, bars2, text1, text2
        else:

            print("frame i=", i, "   ", files[0][i])
            imageAsArray1 = np.array(RT.ImageAsArray(path1 + files[0][i]))
            imageAsArray2 = RT.ImageAsArray(path2 + files[1][i])
            if transparency == True:
                imageAsArray1 = np.ma.masked_inside(imageAsArray1, 0.00001, -0.00001)
                imageAsArray2 = np.ma.masked_inside(imageAsArray2, 0.00001, -0.00001)

            im1 = ax1.pcolormesh(imageAsArray1, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)
            im2 = ax2.pcolormesh(imageAsArray2, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)

            txtTitle = dateList[i] + " (" + satType[i] + pathrowList[i] + ")" + "  diff=" + str(diffDays[i])

            titlte = fig.suptitle(txtTitle, fontsize=16, color="red")

            imageAsArray1_vector = imageAsArray1[~np.isnan(imageAsArray1)].flatten()
            imageAsArray2_vector = imageAsArray2[~np.isnan(imageAsArray2)].flatten()
            (mu, sigma) = norm.fit(imageAsArray1_vector)
            print("1", mu, sigma)
            RMSE1 = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(imageAsArray1_vector))))
            lower_bin = mu - 3 * sigma
            upper_bin = mu + 3 * sigma
            hist1, bins1 = np.histogram(imageAsArray1_vector, range=[lower_bin, upper_bin], density=False, bins=100)
            bars1 = ax3.bar((bins1[:-1] + bins1[1:]) / 2, hist1, align='center', width=(bins1[1] - bins1[0]))

            (mu, sigma) = norm.fit(imageAsArray2_vector)
            print("2", mu, sigma)
            RMSE2 = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(imageAsArray2_vector))))
            lower_bin = mu - 3 * sigma
            upper_bin = mu + 3 * sigma
            hist2, bins2 = np.histogram(imageAsArray2_vector, range=[lower_bin, upper_bin], density=False, bins=100)
            bars2 = ax4.bar((bins2[:-1] + bins2[1:]) / 2, hist2, align='center', width=(bins2[1] - bins2[0]))
            ax3.set_title('Gaussian fit')

            text1 = ax3.text(0.9, 0.9, 'RMSE=' + str(RMSE1),
                             horizontalalignment='right',
                             verticalalignment='top',
                             transform=ax3.transAxes)

            text2 = ax4.text(0.9, 0.9, 'RMSE=' + str(RMSE2),
                             horizontalalignment='right',
                             verticalalignment='top',
                             transform=ax4.transAxes)
            return im1, im2, titlte, bars1, bars2, text1, text2

        ## define the postion and parameters of the color bar

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax2,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1, 0., 1, 1),
                       bbox_transform=ax2.transAxes,
                       borderpad=0.5)
    cbar = plt.colorbar(im2, pad=0.5, cax=axins)  # , ax=axes.ravel().tolist()
    cbar.set_label(label=' Velocity (m/day)', size=8, weight="bold", style="italic", family='times new roman')
    cbar.ax.tick_params(labelsize=8)
    # cbar = plt.colorbar(im1, ax=axes.ravel().tolist(), label=' Velocity (m/day)')

    ani = animation.FuncAnimation(fig, updatefig, frames=len(files[0]) + 3, interval=100, repeat_delay=3000)

    plt.show()
    return


#######################################################################################################################
def Deramping(displacement_field,
              ref_path,
              bandNumber=1,
              snr_filed=None,
              snr_threshold=0.333,
              iterations=3,
              zoom=1,
              showDiagnose=True,
              saveDiagnose=False,
              diagnoseSavingPath="",
              imgIndex=None,
              maskLargeValues=False,
              noData=-32767,
              maskZeros=False,
              max_displacement=-30,
              xlim=[], ylim=[0, 3]):
    """
        Approximate the displacement field with a plane and return the corrected grid, a dictionary with the log of
    the corrections and the correction grid that can be used to transform the slave image.

    Args:
        displacement_field:
        ref_path:
        bandNumber:
        snr_filed:
        snr_threshold:
        iterations:
        zoom:
        showDiagnose:
        saveDiagnose:
        diagnoseSavingPath:
        imgIndex:
        maskLargeValues:
        noData:
        maskZeros:
        max_displacement:
        xlim:
        ylim:

    Returns:

    """
    """

    :param displacement_field: full path to the displacement fields that should be corrected
    :param master_slave: master slave dictionary
    :param ref_path: full path to a reference image, if the displacement field has no georeference the reference
                                of the reference image will be used
    :param snr_filed:  full path to the raster holding the correlation coefficients
    :param snr_threshold:   threshold on the correlation coefficient, pixels below this value will be ignored during the
                                    correction
    :param iterations: number of iterations for iteratively reweighted plane fitting
    :param zoom: optional downsampling, mainly for reducing memory size and processing time
    :param diagnose: if True plotting is activated to illustrate the corrections
    :param maskLargeValues: Bool if we want to mask large values
    :param max_displacement: the values of values to be masked, all displacement values whose absolute magnitude exceeds
                            this value will be masked out for the correction
    :param xlim: x axis limits of the histogram plots
    :param ylim: y axis limits of the histogram plots
    :return: the corrected grid, a dictionary with the log of the corrections and the correction grid
    """

    print("***** OffsetCorrection for image :", displacement_field, "*****")
    # Allocate output dictionary
    out_dict = {}

    # Displacement raster info
    displacementInfo = RT.GetRasterInfo(inputRaster=displacement_field)
    # This allows GDAL to throw Python Exceptions
    gdal.UseExceptions()

    GSD = displacementInfo.get("Resolution")[0]
    bands = displacementInfo.get("NbBands")
    print("=== > GSD, bands", GSD, bands)

    displacementRaster = displacementInfo.get("Raster")
    displBand = displacementRaster.GetRasterBand(bandNumber)

    inputDisp_array = np.array(displBand.ReadAsArray())  ## in meter

    # Get the  georeference of the displacement field
    if displacementInfo.get("MapInfo")[0] == False:
        print('Input field has no georeference. Trying reference image instead...')
        try:
            refRaster = gdal.Open(ref_path)
            geotransform = refRaster.GetGeoTransform()
        except:
            raise IOError('Cannot read georeference from ' + ref_path)
    else:
        geotransform = displacementInfo.get("MapInfo")[1]

    ## Mask values with SNR
    if snr_filed is not None:
        # Import correlation coefficient
        snrBand = displacementRaster.GetRasterBand(2)
        inputSNR_array = np.array(snrBand.ReadAsArray())

        inputDisp_array = filterRT.MaskDispWithSNR(snrField=inputSNR_array, dispArray=inputDisp_array,
                                                   snrThreshold=snr_threshold)
        inputDisp_array = filterRT.MaskDispWithSNR(snrField=inputSNR_array, dispArray=inputDisp_array,
                                                   snrThreshold=snr_threshold)

    ## Mask large displacements
    inputDisp_array = np.ma.masked_where(inputDisp_array == noData, inputDisp_array)

    if maskLargeValues == True:
        inputDisp_array = filterRT.MaskLargeValues(inputArray=inputDisp_array, maxDisplacement=max_displacement)

    if maskZeros == True:
        inputDisp_array = filterRT.MaskZeros(inputArray=inputDisp_array)


    inputDisp_array_before_correction = np.copy(inputDisp_array)
    # print(inputDisp_array_before_correction)
    ### Mask nan values
    # if isinstance(inputDisp_arrayPix, np.ma.MaskedArray):
    inputDisp_array = np.ma.masked_invalid(inputDisp_array)
    # print(inputDisp_array_before_correction)
    ## Get statistical infos on the input
    val_mean = np.ma.mean(inputDisp_array)
    val_std = np.ma.std(inputDisp_array)

    val_min = val_mean - val_std
    val_max = val_mean + val_std
    print("===========================================")
    print("val_mean=%f,val_std=%f" % (val_mean, val_std))
    print("val_min=%f,val_max=%f" % (val_min, val_max))
    print("===========================================")
    out_dict['input'] = {}
    out_dict['input']['name'] = displacement_field
    out_dict['input']['mean'] = val_mean
    out_dict['input']['std'] = val_std
    out_dict['input']['rmse'] = np.ma.sqrt(np.ma.mean(np.square(inputDisp_array)))
    out_dict['input']['max'] = np.ma.max(inputDisp_array)
    out_dict['input']['min'] = np.ma.min(inputDisp_array)

    if isinstance(inputDisp_array, np.ma.MaskedArray):
        out_dict['input']['n_samples'] = inputDisp_array.count()
    else:
        out_dict['input']['n_samples'] = inputDisp_array.size

    print("=== > Input data stat:\n", out_dict)
    inputDisp_arrayGaussianFilterPix, mean_gauss, std_gauss, rmse_gauss = filterRT.IterativeGaussianFit(
        sample=inputDisp_array,
        showDiagnose=showDiagnose,
        saveDiagnose=saveDiagnose,
        diagnoseSavePath=diagnoseSavingPath,
        imgIndex=imgIndex,
        title='Residuals before corrections',
        xlim=xlim, ylim=ylim)

    out_dict['input']['rmse_gauss'] = rmse_gauss
    out_dict['input']['std_gauss'] = std_gauss
    out_dict['input']['mean_gauss'] = mean_gauss
    print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)

    print("\n--- Fit a plane ")
    # return inputDisp_arrayPix
    # Fit plane
    correctionSurface, iter_mad,A = FitPlaneIter(inputDispArray=inputDisp_array, iterations=iterations)

    # correction with plane == deramping
    correctionSurface = correctionSurface.astype(np.int32)
    inputDisp_array = inputDisp_array.astype(np.int32)
    inputDisp_array -= correctionSurface

    cmap = 'seismic'
    # cmap = "rainbow"
    # cmap = 'gist_earth'
    if showDiagnose or saveDiagnose:
        print("Plotiing ....")
        # plot
        fig, axes = plt.subplots(ncols=3, figsize=(10, 3))
        vis_max = max(abs(val_min), abs(val_max))
        vis_min = -vis_max
        # vis_min =-0.5
        # vis_max = 0.5
        axes[0].imshow(inputDisp_array_before_correction, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[0].set_title('Input Image', fontsize=12)
        axes[0].tick_params(labelsize=12)

        # plot
        axes[1].imshow(correctionSurface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[1].set_title('Fitted plane', fontsize=12)
        axes[1].tick_params(labelsize=12)

        im = axes[2].imshow(inputDisp_array, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[2].set_title('After plane correction', fontsize=12)
        axes[2].tick_params(labelsize=12)
        fig.subplots_adjust(right=0.91)

        divider = make_axes_locatable(axes[2])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=12)
        axes[0].axis("off")
        axes[1].axis("off")
        axes[2].axis("off")
        if showDiagnose:
            # plt.show()
            plt.draw()
            plt.pause(1)
        if saveDiagnose:
            fig.savefig(diagnoseSavingPath + "img_" + str(imgIndex) + ".png")
        plt.close()

    # Get infos after plane fitting
    out_dict['plane_corrected'] = {}
    out_dict['plane_corrected']['iter_mad'] = iter_mad
    out_dict['plane_corrected']['mean'] = np.ma.mean(inputDisp_array)
    out_dict['plane_corrected']['std'] = np.ma.std(inputDisp_array)
    out_dict['plane_corrected']['rmse'] = np.sqrt(np.ma.mean(np.square(inputDisp_array)))
    out_dict['plane_corrected']['max'] = np.ma.max(inputDisp_array)
    out_dict['plane_corrected']['min'] = np.ma.min(inputDisp_array)

    inputDisp_arrayGaussianFilter, mean_gauss, std_gauss, rmse_gauss = filterRT.IterativeGaussianFit(
        sample=inputDisp_array,
        showDiagnose=showDiagnose,
        title='Residuals after Deramping',
        xlim=xlim, ylim=ylim)

    out_dict['plane_corrected']['rmse_gauss'] = rmse_gauss
    out_dict['plane_corrected']['std_gauss'] = std_gauss
    out_dict['plane_corrected']['mean_gauss'] = mean_gauss

    print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)

    return inputDisp_array, correctionSurface


if __name__ == '__main__':
    CorrSpotCorrelation()
    # # print("OffsetCorrection")
    # # path = "//home/cosi/2-Data/4-Shishper/3-PostProcessing_Joint_LC_SL/MochwareSide/0-EW_NS_L8/NS_Velocity_Filtered_Gaussian_Masked/"
    # # savingPath = "//home/cosi/2-Data/4-Shishper/3-PostProcessing_Joint_LC_SL/MochwareSide/0-EW_NS_L8/NS_Velocity_Filtered_Gaussian_Masked_offsetCorr/"
    # # MainCorrOffset(inputPath=path,
    # #                savingPath=savingPath, saveGrid=False)
    #
    # # Disp_2_VelocitiesMaps_WithGaussian(path1="//home/cosi/5-Temp_SA/test_coregis/Test_On_Shishper/1-Raw/",
    # #                       path2="//home/cosi/5-Temp_SA/test_coregis/Test_On_Shishper/2-Mod/",
    # #                       title1='EW before',title2="EW after",vMax=4,vMin=-4,transparency=False)
    # #
    # # animationFunction(path1="/home/cosi/2-Data/4-Shishper/2-PostProcessing_LC_SA/149/5-WithOffsetCorrection/4-NS-EW_gaussian/EW_Masked_Gaussian_Masked/",
    # #                       path2="/home/cosi/2-Data/4-Shishper/2-PostProcessing_LC_SA/149/5-WithOffsetCorrection/4-NS-EW_gaussian/NS_Masked_Gaussian_Masked/",
    # #                       title1='EW',title2="NS",vMax=4,vMin=-4,transparency=True)
    # imgPath = "//home/cosi/2-Data/5-Shaybah/S22Subsets/Correlation_64_32_6/"
    # savingPath = "/home/cosi/2-Data/5-Shaybah/S22Subsets/Correlation_64_32_6/NS_OffsetCorrection/"
    #
    # # output = np.ma.array(OffsetCorrectionPixelBased_SA(dispMap=imgPath,
    # #                                                    bandNumber=1,
    # #                                                    snrMap=imgPath,
    # #                                                    snrBandNb=3,
    # #                                                    snrThreshold=0.6,
    # #                                                    maskLargeValues=50,
    # #                                                    maskZeros=False,
    # #                                                    showDiagnose=True,
    # #                                                    saveDiagnose=False,
    # #                                                    diagnoseSavingPath=""))
    # #
    # # RT.WriteRaster(refRasterPath=imgPath, newRasterPath=os.path.join(savingPath, "NS.tif"),
    # #                Listarrays=[output],
    # #                numberOfBands=1)
    # imgs = FileRT.FilesInDirectory(path=imgPath)
    # imgsNewPath = []
    # for index, img_ in enumerate(imgs):
    #     if ".tif" in img_ and ".hdr" not in img_ and ".aux.xml" not in img_:
    #         imgsNewPath.append(imgPath + img_)
    #
    # for index, img_ in enumerate(imgsNewPath):
    #     print("\n ------------ image number :", index + 1, "---------------\n")
    #
    #     # output = OffsetCorrectionPixelBased(displacement_field=imgPath, master_slave=imgPath,
    #     #                                     ref_path=imgPath, snr_filed=None, saveDiagnose=saveGrid,
    #     #                                     diagnoseSavingPath=gridSavingPath, imgIndex=index + 1,
    #     #                                     maskLargeValues=False, maskZeros=True, showDiagnose=False)
    #     # rasterinfo = RT.GetRasterInfo(inputRaster=img_)
    #     # print(rasterinfo)
    #     corrInfo = Correlation_misc.CorrFileInfo(img_, oneDic=True)
    #
    #     title = corrInfo["MasterDate"] + "_" + corrInfo["MasterSatellite"] + "_" + corrInfo["SlaveDate"] + "_" + \
    #             corrInfo["MasterSatellite"] + "_NS.tif"
    #     output = np.ma.array(OffsetCorrectionPixelBased_SA(dispMap=img_,
    #                                                        bandNumber=2,
    #                                                        snrMap=img_,
    #                                                        snrBandNb=3,
    #                                                        snrThreshold=0.6,
    #                                                        maskLargeValues=50,
    #                                                        maskZeros=False,
    #                                                        showDiagnose=False,
    #                                                        saveDiagnose=False,
    #                                                        diagnoseSavingPath=""))
    #
    #     metaData = corrInfo
    #     metaData["CorrBand"] = "N/S"
    #     RT.WriteRaster(refRasterPath=img_, newRasterPath=os.path.join(savingPath, title),
    #                    Listarrays=[output],
    #                    numberOfBands=1, metaData=corrInfo)
