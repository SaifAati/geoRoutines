import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import norm
from scipy import ndimage
import geoRoutines.FilesCommandRoutine as FileRT
import geoRoutines.Routine as RT
import geoRoutines.Plotting.Plotting_Routine as Plot
import gdal
from pathlib import Path


def MaskDisplacementMap(maskFile, map2Mask, bandNumber=1, visualization=False):
    maskArray = RT.ImageAsArray(RT.GetRasterInfo(maskFile))
    print(maskArray)
    mask = maskArray > 0

    inputArray = RT.ImageAsArray(imageInfo=RT.GetRasterInfo(map2Mask), bandNumber=bandNumber)
    inputArray_masked = np.empty(inputArray.shape)
    inputArray_masked = np.ma.masked_array(np.ma.masked_invalid(inputArray), mask=mask)

    return inputArray_masked


def MultiScaleFilterWithSNR(inputData, outputData="", thresholdVelovity=7, thresholdMedian=1.5, step=3,
                            multiScaleFactor=2,
                            withSNR=True):
    import Multi_Scale_Filter as MF
    filter = MF.Filtering()

    workspacePath = Path(inputData).parents[1]

    if not outputData:
        workspacePath = Path(inputData).parents[1]
        folders = FileRT.FilesInDirectory(path=workspacePath)
        print(folders)
        if "MultiScaleFilter" in folders:
            tempPath = os.path.join(workspacePath, "MultiScaleFilter")
        else:
            tempPath = FileRT.CreateDirectory(directoryPath=workspacePath, folderName="MultiScaleFilter")

        outputData = FileRT.CreateDirectory(directoryPath=tempPath,
                                            folderName=os.path.basename(os.path.normpath(inputData)) + "_Filtered")

    filter.MultiScaleFilterWithSNR(inputData=inputData, outputData=outputData, thresholdVelovity=thresholdVelovity,
                                   thresholdMedian=thresholdMedian, step=step,
                                   multiScaleFactor=multiScaleFactor, withSNR=withSNR)

    return


def MaskDispWithSNR(snrArray, dispArray, snrThreshold=0.5, visualization=False):
    """

    :param snrArray:
    :param dispArray:
    :param snrThreshold:
    :param visualization:
    :return:
    """
    dispArray_origin = np.array(dispArray)
    stat_beforeMasking = {"Mean": np.nanmean(dispArray_origin), "Max": np.nanmax(dispArray_origin),
                          "Min": np.nanmin(dispArray_origin)}
    print("stat_beforeSNRMasking:", stat_beforeMasking)
    print('--- Masking displacement field with correlation coefficients (SNR)')
    ## Create mask
    snrField_Mask = snrArray < snrThreshold
    ## Apply mask
    dispArray = np.ma.masked_array(dispArray, mask=snrField_Mask)
    stat_afterMasking = {"Mean": np.ma.mean(dispArray), "Max": np.ma.max(dispArray), "Min": np.ma.min(dispArray)}
    print("stat_afterSNRMasking:", stat_afterMasking)

    ## Check if there are "enough" valid pixels
    if snrField_Mask.sum() / snrField_Mask.size > 0.5:
        warnings.warn("Less than 50% of the pixels are valid.", UserWarning)
        warnings.warn("The reason might be low correlation or small overlap among the input images.", UserWarning)
    del snrField_Mask
    if visualization:
        Plot.Plot_2_maps(imgArray1=dispArray_origin, imgArray2=dispArray,
                         backgroundImg=None, figureTitle="Masking Values using SNR")
    return dispArray


def MaskLargeValues(inputArray, maxDisplacement,fillValue=None, visualization=False):
    """
    :param inputArray:
    :param maxDisplacement:
    :param visualization:
    :return:
    """
    inputArray = np.ma.masked_invalid(inputArray)
    inputArray_origin = np.ma.masked_invalid(inputArray)
    # print(inputArray_origin)
    stat_beforeMasking = {"Mean": np.nanmean(inputArray_origin), "Max": np.nanmax(inputArray_origin),
                          "Min": np.nanmin(inputArray_origin)}
    print("stat_beforeLargeValuesMasking:", stat_beforeMasking)
    print('--- Masking large values, max displacement:', maxDisplacement, -maxDisplacement)
    mask_large_values = np.absolute(inputArray) > abs(maxDisplacement)
    if np.any(mask_large_values):
        print("Number of pixel masked:", mask_large_values.sum(),
              "Corresponding to a ratio:", mask_large_values.sum() * 100 / mask_large_values.size, "%")
        res = np.ma.masked_where(np.absolute(inputArray) > abs(maxDisplacement), inputArray)
        stat_afterMasking = {"Mean": np.nanmean(res), "Max": np.nanmax(res), "Min": np.nanmin(res)}
        print("stat_afterLargeValuesMasking:", stat_afterMasking)
        if visualization:
            Plot.Plot_2_maps(imgArray1=inputArray_origin, imgArray2=res,
                             backgroundImg=None,
                             vMin=-maxDisplacement, vMax=maxDisplacement)
        if fillValue != None:
            res = res.filled(fill_value=fillValue)

        return res
    else:
        print("No large values have been masked")
        return inputArray


def MaskZeros(inputArray):
    """
    Returns a masked array if max_displacement is exceeded anywhere on the map
    :param inputArray:
    :param maxDisplacement:
    :return:
    """
    print("--- Masking zeros ---")
    mask_zeros_values = np.ma.masked_where(inputArray == 0, inputArray)
    if np.any(mask_zeros_values):
        print("Number of pixel masked:", mask_zeros_values.sum(),
              "Corresponding to a ratio:", mask_zeros_values.sum() * 100 / mask_zeros_values.size, "%")
        return mask_zeros_values
    else:
        print("No zeros values have been masked")
        return inputArray


def IterativeGaussianFit(sample, max_iterations=100,
                         min_change_sigma=0.005,
                         sigma_level=2.575,
                         showDiagnose=False,
                         saveDiagnose=False,
                         diagnoseSavePath="",
                         imgIndex=None,
                         title='Sample distribution',
                        xLabel='Displacement [m]',
                         xlim=[], ylim=[0, 3],nbins =50,fontSize=16):
    """
    Iterative fits the sample population with a Gaussian and removes outliers outside the defined confidence interval
    :param sample: Array with the samples
    :param max_iterations: Maximum number of iterations
    :param min_change_sigma: Minimum change in sigma to continue fitting
    :param sigma_level: sigma_level*sigma > defines threshold for outliers, 2.575 corresponds to a 99% confidence level
    :param diagnose: if True a plot will be shown
    :param title: title of the upper figure
    :param xlim: lower limit, upper limit of the x - axis
    :param ylim: lower limit, upper limit of the y - axis
    :return: 1D array of the sample population with outliers and the RMSE
    """

    print("\n---  Start:Iterative gaussian Fit ")
    # # Remove mask and array to vector
    # if isinstance(sample, np.ma.MaskedArray):  # check if the sample was masked using the class numpy.ma.MaskedArray
    #     sample = sample.compressed()  ## return all the non-masked values as 1-D array
    #
    # else:
    #     if sample.ndim > 1:  # if the dimension of the array more than 1 transform it to 1-D array
    #         sample = sample.flatten()
    #
    # # Estimate initial sigma and RMSE
    # (mu, sigma) = norm.fit(sample)
    # print("mu, sigma", mu, sigma)
    # RMSE = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(sample))))
    dist = RT.Distribution(inputArray=sample)
    RMSE =dist.RMSE
    mu = dist.mu
    sigma = dist.sigma

    if showDiagnose or saveDiagnose:  # plot
        print("Plotting ....")
        # fig2 = plt.figure()  # figsize=(10, 10))
        fig, axes = plt.subplots(nrows=2, figsize=(10, 10))
        lower_bin = dist.mu - 3 * dist.sigma
        upper_bin = dist.mu + 3 * dist.sigma
        hist, bins = np.histogram(dist.sample, range=[lower_bin, upper_bin], density=False, bins=nbins)
        bars = axes[0].bar((bins[:-1] + bins[1:]) / 2, hist, align='center', width=(bins[1] - bins[0]), color="blue",
                       edgecolor='black')
        # sns.distplot(dist.sample, hist=False, kde=True, bins=nbins, color='darkblue',
        #              hist_kws={'edgecolor': 'black'}, kde_kws={'shade': True, 'linewidth': 4})
        title_ = title
        axes[0].set_title(title_)
        # axes[0].set_ylim(ylim[0], ylim[1])
        if not any(xlim):
            axes[0].set_xlim(float(dist.mean) - 4 * float(dist.std), float(dist.mean) + 4 * float(dist.std))
        else:
            axes[0].set_xlim(xlim[0], xlim[1])
        txt = "Mean=" + str(dist.mean) \
              + "\nStd=" + str(dist.std) \
              + "\nMedian=" + str(dist.median) \
              + "\nNMAD=" + str(dist.nmad)
        txt1 = axes[0].text(0.75, 0.9, txt,
                        horizontalalignment='left',
                        verticalalignment='top',
                        transform=axes[0].transAxes, fontsize=fontSize)
        axes[0].set_xlabel(xLabel, fontsize=fontSize)
        axes[0].set_ylabel(' Number of Samples ', fontsize=fontSize)



        txt2 = axes[1].text(0.9, 0.9, 'RMSE=' + str(dist.RMSE),
                            horizontalalignment='right',
                            verticalalignment='top',
                            transform=axes[1].transAxes)
        y = norm.pdf(bins, dist.mu, dist.sigma)

        axes[1].plot(bins, y, 'r--', linewidth=2)
        text2 = axes[1].set_title('Gaussian fit', fontsize=12)
        axes[1].set_ylim(0, y.max() + 0.1)
        if not any(xlim):
            axes[1].set_xlim(float(dist.mean) - 4 * float(dist.std), float(dist.mean) + 4 * float(dist.std))
        else:
            axes[1].set_xlim(xlim[0], xlim[1])


    for i in range(max_iterations):
        sigmaStart = dist.sigma  # update starting sigma
        # Remove outliers
        sample = np.extract(sample < mu + sigma_level * sigma, sample)  # remove upper tail
        sample = np.extract(sample > mu - sigma_level * sigma, sample)  # remove lower tail
        outliers = sample.sum()
        print("Number of outliers=", outliers, "with a ratio= ", outliers * 100 / sample.size, "%")
        # Re-estimate Gaussian and RMSE
        # (mu, sigma) = norm.fit(sample)
        # RMSE = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(sample))))
        dist = RT.Distribution(sample)
        mu = dist.mu
        RMSE = dist.RMSE
        sigma = dist.sigma
        if showDiagnose or showDiagnose:  # Update plot
            txt2.set_visible(False)



            hist, bins = np.histogram(dist.sample, range=[lower_bin, upper_bin], density=False, bins=nbins)
            # bars = axes[0].bar((bins[:-1] + bins[1:]) / 2, hist, align='center', width=(bins[1] - bins[0]),
            #                    color="blue",
            #                    edgecolor='black')
            y = norm.pdf(bins, mu, sigma)
            for rect, h in zip(bars, hist):
                rect.set_height(h)

            for line in axes[0].lines:
                line.set_xdata(bins)
                line.set_ydata(y)

            text2.set_text("Gaussian fit, iteration number: " + str(i))
            txt2 = axes[1].text(0.9, 0.9, 'RMSE=' + str(RMSE),
                                horizontalalignment='right',
                                verticalalignment='top',
                                transform=axes[1].transAxes)
            axes[1].set_ylim(0, y.max() + 0.1)
            if showDiagnose:
                # plt.pause(1)
                plt.show()
            if showDiagnose:
                fig.savefig(diagnoseSavePath + "Img_" + str(imgIndex) + "_GauussianFit_Iter_" + str(i) + ".png")

        if sigmaStart - sigma < min_change_sigma:
            if showDiagnose or showDiagnose:
                plt.close()
            print("--- END: Iterative gaussian Fit: number of iteration", i,
                  "Ratio of ouliers:", outliers * 100 / sample.size, "%, Rmse:", RMSE, "\n")
            return sample, mu, sigma, RMSE
    print("--- END: Iterative gaussian Fit \n")
    if showDiagnose or saveDiagnose:
        plt.close()
    return sample, mu, sigma, RMSE


def GaussianSmoothing_Batch(inputData, outputData="", sigma=1, mode="nearest"):
    imgs = FileRT.FilesInDirectory(inputData)
    workspacePath = Path(inputData).parents[1]

    if not outputData:
        workspacePath = Path(inputData).parents[1]
        folders = FileRT.FilesInDirectory(path=workspacePath)
        print(folders)
        if "GaussianFilter" in folders:
            tempPath = os.path.join(workspacePath, "GaussianFilter")
        else:
            tempPath = FileRT.CreateDirectory(directoryPath=workspacePath, folderName="GaussianFilter")

        outputData = FileRT.CreateDirectory(directoryPath=tempPath,
                                            folderName=os.path.basename(os.path.normpath(inputData)) + "_Gaussian")

    for im in range(len(imgs)):
        if ".tif" in imgs[im]:
            img_ = imgs[im]
            rasterInfo = RT.GetRasterInfo(inputRaster=inputData + img_)
            velocityMapArray = RT.ImageAsArray(imageInfo=rasterInfo, bandNumber=1)

            velocityMapFinal = gaussian_filter(velocityMapArray, sigma=sigma, mode=mode)
            ### save the smoothed velocity map as raster
            print("max=", velocityMapFinal[~np.isnan(velocityMapFinal)].max())
            print("min=", velocityMapFinal[~np.isnan(velocityMapFinal)].min())
            metaData = rasterInfo["MetaData"]
            metaData["GaussianFilterParams"] = " Sigma= " + str(sigma) + ", Mode= " + mode
            RT.WriteRaster(refRasterPath=inputData + img_, newRasterPath=outputData + img_,
                           Listarrays=[velocityMapFinal], numberOfBands=1, metaData=metaData)
    return


def Build_NaN_Mask(inputData, saving=True):
    rasterInfo = RT.GetRasterInfo(inputRaster=inputData)
    logicalMasksArr = []
    for band_ in range(rasterInfo["NbBands"]):
        bandAsArray = RT.ImageAsArray(imageInfo=rasterInfo, bandNumber=band_ + 1)
        temp = np.ma.masked_invalid(bandAsArray)
        print("Nb of NaNs=", temp.size - temp.count())
        mask_ = temp.mask
        logicalMask = ~mask_
        logicalMasksArr.append(logicalMask.astype(np.int))

    if saving:
        outputPath = os.path.join(os.path.dirname(inputData), "NaN_RasterMask")
        RT.WriteRaster(refRasterPath=inputData, newRasterPath=outputPath, Listarrays=logicalMasksArr,
                       numberOfBands=rasterInfo["NbBands"], metaData={"Info": " NaNs mask"})
        return outputPath
    else:
        return logicalMasksArr


def Convert_NaN_to_Value(inputData, value=0.0, saving=True):
    rasterInfo = RT.GetRasterInfo(inputRaster=inputData)
    convBandsArr = []
    for band_ in range(rasterInfo["NbBands"]):
        bandAsArray = RT.ImageAsArray(imageInfo=rasterInfo, bandNumber=band_ + 1)
        temp = np.ma.masked_invalid(bandAsArray)
        print("Nb of NaNs=", temp.size - temp.count())
        mask_ = temp.mask
        bandAsArray[mask_] = value
        convBandsArr.append(bandAsArray)

    if saving:
        outputPath = os.path.join(os.path.dirname(inputData), "Data_Converting_NaNs")
        RT.WriteRaster(refRasterPath=inputData, newRasterPath=outputPath, Listarrays=convBandsArr,
                       numberOfBands=rasterInfo["NbBands"], metaData={"Info": " Vonverting NaNs as " + str(value)})
        return outputPath
    else:
        return convBandsArr


def MedianFilter_raster(rasterPath, th=5, display=False, noData=None):
    rasterInfo = RT.GetRasterInfo(rasterPath)
    rasterArry = RT.ImageAsArray(rasterInfo)
    print(rasterArry.min())
    if noData != None:
        rasterArry = np.ma.masked_where(rasterArry == noData, rasterArry)

    im_med = ndimage.median_filter(rasterArry, th)
    if isinstance(rasterArry, np.ma.masked_array):
        mask = np.ma.getmask(rasterArry)
        im_med = np.ma.masked_array(im_med, mask=mask)
    RT.WriteRaster(refRasterPath=rasterPath,
                   newRasterPath=os.path.join(os.path.dirname(rasterPath),
                                              Path(rasterPath).stem + "_med" + str(th) + ".tif"), Listarrays=[im_med],
                   numberOfBands=1)
    if display:
        plt.figure()
        plt.subplot(121)
        plt.imshow(rasterArry, cmap=plt.cm.gray)  # , vmin=40, vmax=220)
        plt.axis('off')
        # plt.title('noisy', fontsize=20)
        plt.subplot(122)
        plt.imshow(im_med, cmap=plt.cm.gray)  # , vmin=40, vmax=220)
        plt.axis('off')
        # plt.title('Gaussian filter', fontsize=20)
        # plt.subplot(133)
        # plt.imshow(med_denoised, cmap=plt.cm.gray, vmin=40, vmax=220)
        # plt.axis('off')
        # plt.title('Median filter', fontsize=20)

        # plt.subplots_adjust(wspace=0.02, hspace=0.02, top=0.9, bottom=0, left=0,
        #                     right=1)
        plt.show()


if __name__ == '__main__':
    rasterPath = os.path.join("//home/cosicorr/0-WorkSpace/3D-Correlation_project/DEMS_correlation",
                              "1-Pre_dem_subset_HH.tif")
    # rasterPath = os.path.join("/home/cosicorr/0-WorkSpace/3D-Correlation_project/DEMS_correlation/Cosi-Corr",
    #                           "PRE_dem_subset_HH.tif")
    # MedianFilter_raster(rasterPath=rasterPath,th=9,display=False,noData=0)
    # path = "/home/cosi/2-Data/4-Shishper/3-PostProcessing_Joint_LC_SL/ShishperSide/3-PCA_envi/DataCube.tif"
    # Convert_NaN_to_Value(inputData=path)

    # a =np.array([[1,2,3],[2,np.nan,10]])
    # print(a)
    # temp = np.ma.masked_invalid(a)
    # mask = temp.mask
    #
    #
    # a[mask] = 0
    #
    # print(a)

    rasterPath = "/home/cosicorr/0-WorkSpace/Test_Data_ElMayor/Correlation_Set2_test/Test3/All/EW_Stack/EW_All_T3_3DD_det"
    rasterInfo = RT.GetRasterInfo(rasterPath)
    dataCube = RT.MultiBandsRaster2Array(rasterInfo)
    dataCube_filtered = []
    # temp = MaskLargeValues(dataCube[0,:,:],maxDisplacement=-1,fillValue=-32767)
    # print(temp)
    for i in range(dataCube.shape[0]):
        dataCube_filtered.append(MaskLargeValues(inputArray=dataCube[i, :, :], maxDisplacement=12,fillValue=-32767))
    #
    # print(dataCube_filtered)
    RT.WriteRaster(refRasterPath=rasterPath,
                   newRasterPath=os.path.join(os.path.dirname(rasterPath), Path(rasterPath).stem + "_filt.tif"),
                   Listarrays=dataCube_filtered, numberOfBands=len(dataCube_filtered),
                   descriptions=rasterInfo["BandInfo"],noData=-32767,dtype=gdal.GDT_Float32)
