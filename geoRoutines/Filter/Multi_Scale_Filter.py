import numpy as np
from osgeo import gdal, osr
import sys, os, warnings
from matplotlib import pyplot as plt

import geoRoutines.FilesCommandRoutine as FileRT
import geoRoutines.georoutines as geoRT

from pathlib import Path


class Filtering:
    def __init__(self):
        return

    def WeightedMeanFilter(self, velocityMapArray, snrMapArray, thresholdVelocity=7, thresholdMedian=1.1, step=3):
        weightArray = velocityMapArray * snrMapArray
        row, column = velocityMapArray.shape
        velocityMapFinal = np.zeros((row, column))

        windTmp = np.zeros((step, step))
        i = j = 0
        while i < row:
            while j < column:
                for a in range(step):
                    for b in range(step):
                        windTmp[a, b] = velocityMapArray[i + a, j + b]
                        windTmp[a, b] = weightArray[i + a, j + b]

                windTmp_ = windTmp

                for a in range(step):
                    for b in range(step):
                        if windTmp[a, b] > thresholdVelocity:
                            windTmp[a, b] = np.nan
                            windTmp[a, b] = 0

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    med = np.median(windTmp[~np.isnan(windTmp)])

                for a in range(step):
                    for b in range(step):
                        if windTmp[a, b] > med * thresholdMedian and med != np.nan:
                            windTmp[a, b] = np.nan
                            # windTmp[a, b] = 0

                # Computing the mean local average
                numberOfNaN = len(windTmp[np.isnan(windTmp)])
                if numberOfNaN == 0:
                    localAvg = windTmp.mean()

                if numberOfNaN != 0 and numberOfNaN < step ** 2:
                    localAvg = windTmp[~np.isnan(windTmp)].mean()

                if numberOfNaN == step ** 2:
                    localAvg = np.nan

                if localAvg != 0 and localAvg > 6:
                    print("windTmp=\n", windTmp)
                    print("mean=", localAvg)
                    print("median=", med)
                    print("\n")

                for a in range(step):
                    for b in range(step):
                        velocityMapFinal[i + a, j + b] = localAvg

                j = j + step
                # print("j=",j,"i=",i)
                if (j + step > column):
                    # print(True)
                    j = np.Infinity
            i = i + step
            j = 0
            # print("\ni=",i)
            if (i + step > row):
                # print(True,"Out")
                i = np.Infinity

        return (velocityMapFinal)

    # def WeightedMeanFilterV2(self, velocityMapArray, snrMapArray, thresholdVelocity=7, thresholdMedian=1.5, step=3):
    #     weightArray = velocityMapArray * snrMapArray
    #     row, column = velocityMapArray.shape
    #     velocityMapFinal = np.array(weightArray)
    #
    #     i = j = 0
    #     while i < row:
    #         while j < column:
    #             windTmp = np.zeros((step, step))
    #             windTmp_ = np.zeros((step, step))
    #             windTmp_VelocityThreshold = np.zeros((step, step))
    #             # Compute windTmp
    #             for a in range(step):
    #                 for b in range(step):
    #                     # windTmp[a, b] = velocityMapArray[i + a, j + b]
    #                     windTmp[a, b] = weightArray[i + a, j + b]
    #
    #             windTmp_ = np.array(windTmp)
    #             # Apply Velocity Threshold
    #             np.warnings.simplefilter("ignore", category=RuntimeWarning)
    #             windTmp[(~np.isnan(windTmp)) & (windTmp > thresholdVelocity)] = np.nan
    #             # for a in range(step):
    #             #     for b in range(step):
    #             #         if windTmp[a, b] > thresholdVelocity:
    #             #             windTmp[a, b] = np.nan
    #             #             # windTmp[a, b] = 0
    #             #             # windTmp[a,b]=
    #             windTmp_VelocityThreshold = np.array(windTmp)
    #
    #             # Compute median
    #             # with warnings.catch_warnings():
    #             np.warnings.simplefilter("ignore", category=RuntimeWarning)
    #             med = np.median(windTmp[~np.isnan(windTmp)])
    #
    #             # Apply statistical median threshold
    #             if ~np.isnan(med):
    #                 for a in range(step):
    #                     for b in range(step):
    #                         if ~np.isnan(windTmp[a, b]):
    #                             if windTmp[a, b] > med + thresholdMedian or windTmp[a, b] < med - thresholdMedian:
    #                                 windTmp[a, b] = np.nan
    #                             # windTmp[a, b] = 0
    #
    #             windTmp_AfterStatisticalFilter = np.array(windTmp)
    #             # Computing number of NaN values into Window
    #             numberOfNaN = len(windTmp[np.isnan(windTmp)])
    #
    #             # Build FinalVelocity by changing Filtered Values
    #             if numberOfNaN < step ** 2:
    #                 locAvg = np.mean(windTmp[~np.isnan(windTmp)])
    #                 for a in range(step):
    #                     for b in range(step):
    #                         if np.isnan(windTmp[a, b]):
    #                             # windTmp[a,b]=med
    #                             # velocityMapFinal[i + a, j + b] = med
    #                             windTmp[a, b] = locAvg
    #                             velocityMapFinal[i + a, j + b] = locAvg
    #
    #             if numberOfNaN == step ** 2:
    #                 locAvg = np.nan
    #                 for a in range(step):
    #                     for b in range(step):
    #                         velocityMapFinal[i + a, j + b] = np.nan
    #                         windTmp[a, b] = np.nan
    #
    #             ## prints
    #             # if numberOfNaN < step ** 2:
    #             #     print("windTmp first =\n",windTmp_,"\n")
    #             #     print("windTmp afterVelocityThreshold =\n", windTmp_VelocityThreshold)
    #             #     print("median = ",med)
    #             #     print("windTmp_AfterStatisticalFilter=\n",windTmp_AfterStatisticalFilter)
    #             #     print("number of Nan=",numberOfNaN)
    #             #     print("WindTmpNew=\n",windTmp,"\n\n")
    #             j = j + step
    #             # print("j=",j,"i=",i)
    #             if (j + step > column):
    #                 # print(True)
    #                 j = np.Infinity
    #         i = i + step
    #         j = 0
    #         # print("\ni=",i)
    #         if (i + step > row):
    #             # print(True,"Out")
    #             i = np.Infinity
    #
    #     return (velocityMapFinal)
    #
    # def WeightedMeanFilterV3(self, velocityMapArray, snrMapArray, thresholdVelocity, thresholdMedian, step):
    #     weightArray = velocityMapArray * snrMapArray
    #     row, column = velocityMapArray.shape
    #     velocityMapFinal = np.array(weightArray)
    #
    #     i = j = 0
    #     while i < row:
    #         while j < column:
    #             windTmp = np.zeros((step, step))
    #             windTmp_ = np.zeros((step, step))
    #             windTmp_VelocityThreshold = np.zeros((step, step))
    #             # Compute windTmp
    #             for a in range(step):
    #                 for b in range(step):
    #                     # windTmp[a, b] = velocityMapArray[i + a, j + b]
    #                     windTmp[a, b] = weightArray[i + a, j + b]
    #
    #             windTmp_ = np.array(windTmp)
    #             # Apply Velocity Threshold
    #             np.warnings.simplefilter("ignore", category=RuntimeWarning)
    #             windTmp[(~np.isnan(windTmp)) & (windTmp > thresholdVelocity)] = np.nan
    #
    #             windTmp_VelocityThreshold = np.array(windTmp)
    #
    #             # Compute median
    #             # with warnings.catch_warnings():
    #             np.warnings.simplefilter("ignore", category=RuntimeWarning)
    #             med = np.median(windTmp[~np.isnan(windTmp)])
    #
    #             # Apply statistical median threshold
    #             if ~np.isnan(med):
    #                 for a in range(step):
    #                     for b in range(step):
    #                         if ~np.isnan(windTmp[a, b]):
    #                             if windTmp[a, b] > med + thresholdMedian or windTmp[a, b] < med - thresholdMedian:
    #                                 windTmp[a, b] = np.nan
    #                             # windTmp[a, b] = 0
    #
    #             windTmp_AfterStatisticalFilter = np.array(windTmp)
    #             # Computing number of NaN values into Window
    #             numberOfNaN = len(windTmp[np.isnan(windTmp)])
    #
    #             # Build FinalVelocity by changing Filtered Values
    #             for a in range(step):
    #                 for b in range(step):
    #                     velocityMapFinal[i + a, j + b] = windTmp[a, b]
    #
    #             ## prints
    #             # if numberOfNaN < step ** 2:
    #             #     print("windTmp first =\n",windTmp_,"\n")
    #             #     print("windTmp afterVelocityThreshold =\n", windTmp_VelocityThreshold)
    #             #     print("median = ",med)
    #             #     print("windTmp_AfterStatisticalFilter=\n",windTmp_AfterStatisticalFilter)
    #             #     print("number of Nan=",numberOfNaN)
    #             #     print("WindTmpNew=\n",windTmp,"\n\n")
    #             j = j + step
    #             # print("j=",j,"i=",i)
    #             if (j + step > column):
    #                 # print(True)
    #                 j = np.Infinity
    #         i = i + step
    #         j = 0
    #         # print("\ni=",i)
    #         if (i + step > row):
    #             # print(True,"Out")
    #             i = np.Infinity
    #
    #     return (velocityMapFinal)

    def WeightedMeanFilterV4(self, velocityMapArray, snrMapArray, thresholdVelocity, thresholdMedian, step,
                             multiScaleFactor):
        """
        We add the multi-scale factor
        :param velocityMapArray:
        :param snrMapArray:
        :param thresholdVelocity:
        :param thresholdMedian:
        :param step:
        :return:
        """
        weightArray = velocityMapArray * snrMapArray
        row, column = velocityMapArray.shape
        velocityMapFinal = np.array(weightArray)

        i = j = 0
        while i < row:
            while j < column:
                windTmp = np.zeros((step, step))
                windTmp_ = np.zeros((step, step))
                windTmp_VelocityThreshold = np.zeros((step, step))
                # Compute windTmp
                for a in range(step):
                    for b in range(step):
                        # windTmp[a, b] = velocityMapArray[i + a, j + b]
                        windTmp[a, b] = weightArray[i + a, j + b]

                windTmp_ = np.array(windTmp)
                # Apply Velocity Threshold
                np.warnings.simplefilter("ignore", category=RuntimeWarning)
                windTmp[(~np.isnan(windTmp)) & (windTmp > thresholdVelocity)] = np.nan
                windTmp[(~np.isnan(windTmp)) & (windTmp < -thresholdVelocity)] = np.nan

                windTmp_VelocityThreshold = np.array(windTmp)

                # Compute median
                # with warnings.catch_warnings():
                np.warnings.simplefilter("ignore", category=RuntimeWarning)
                med = np.median(windTmp[~np.isnan(windTmp)])

                # Apply statistical median threshold
                if ~np.isnan(med):
                    for a in range(step):
                        for b in range(step):
                            if ~np.isnan(windTmp[a, b]):
                                if windTmp[a, b] > med + thresholdMedian or windTmp[a, b] < med - thresholdMedian:
                                    windTmp[a, b] = np.nan
                                # windTmp[a, b] = 0

                windTmp_AfterStatisticalFilter = np.array(windTmp)
                # Computing number of NaN values into Window
                numberOfNaN = len(windTmp[np.isnan(windTmp)])

                # Build FinalVelocity by changing Filtered Values
                for a in range(step):
                    for b in range(step):
                        velocityMapFinal[i + a, j + b] = windTmp[a, b]

                ## prints
                # if numberOfNaN < step ** 2:
                #     print("windTmp first =\n",windTmp_,"\n")
                #     print("windTmp afterVelocityThreshold =\n", windTmp_VelocityThreshold)
                #     print("median = ",med)
                #     print("windTmp_AfterStatisticalFilter=\n",windTmp_AfterStatisticalFilter)
                #     print("number of Nan=",numberOfNaN)
                #     print("WindTmpNew=\n",windTmp,"\n\n")
                j = j + step
                # print("j=",j,"i=",i)
                if (j + step > column):
                    # print(True)
                    j = np.Infinity
            i = i + step
            j = 0
            # print("\ni=",i)
            if (i + step > row):
                # print(True,"Out")
                i = np.Infinity

        newWindowSize = step + multiScaleFactor

        numberOfNaNPostFiltering = len(velocityMapFinal[np.isnan(velocityMapFinal)])
        print("numberOfNaN after Filtering=", numberOfNaNPostFiltering)
        count = 0
        velocityMapFinal_withboundry = np.empty((velocityMapFinal.shape[0] + 2 * int(newWindowSize / 2),
                                                 velocityMapFinal.shape[1] + 2 * int(newWindowSize / 2)))
        velocityMapFinal_withboundry[:] = np.nan
        print(velocityMapFinal.shape, velocityMapFinal_withboundry.shape)
        velocityMapFinal_IntilizeNan = np.array(velocityMapFinal)
        offSet = int(newWindowSize / 2)
        print("offset=", offSet)

        for l in range(offSet, velocityMapFinal.shape[0] + offSet):
            for c in range(offSet, velocityMapFinal.shape[1] + offSet):
                velocityMapFinal_withboundry[l, c] = velocityMapFinal[l - offSet, c - offSet]

        for l in range(offSet, velocityMapFinal.shape[0] + offSet):
            for c in range(offSet, velocityMapFinal.shape[1] + offSet):
                if np.isnan(velocityMapFinal_withboundry[l, c]):
                    count += 1

                    def crop_center(img, index, windowSize):
                        x, y = index
                        cropx, cropy = windowSize
                        startx = x - cropx // 2
                        starty = y - cropy // 2
                        return img[startx:x + cropx // 2 + 1, starty:y + cropy // 2 + 1]

                    windTmp__ = crop_center(velocityMapFinal_withboundry, (l, c), (newWindowSize, newWindowSize))
                    n = len(windTmp__[np.isnan(windTmp__)])
                    if n != newWindowSize ** 2:
                        moy = np.average(windTmp__[~np.isnan(windTmp__)])
                        # print(moy)
                        velocityMapFinal_IntilizeNan[l - offSet, c - offSet] = moy

        # from matplotlib import pyplot as plt
        # # plt.imshow(velocityMapFinal)
        # fig =plt.figure()
        # plt.imshow(velocityMapFinal_withboundry)
        # fig2 = plt.figure()
        # plt.imshow(velocityMapFinal_IntilizeNan)
        # plt.show()
        # print("count=",count)

        return (velocityMapFinal_IntilizeNan)

    def MultiScaleFilterWithSNR(self, inputData, outputData, thresholdVelovity, thresholdMedian, step,
                                multiScaleFactor=2, withSNR=False):
        """

        :param inputData:
        :param outputData:
        :param thresholdVelovity:
        :param thresholdMedian:
        :param step:
        :param multiScaleFactor:
        :param withSNR:
        :return:
        """
        imgs = FileRT.FilesInDirectory(inputData)
        bandNb = geoRT.RasterInfo(inputRaster=os.path.join(inputData, imgs[0])).nbBand
        print(bandNb)
        if bandNb == 2:
            # for img_ in imgs :
            for im in range(len(imgs)):
                img_ = imgs[im]
                print("\n", img_)
                rasterInfo = geoRT.RasterInfo(inputRaster=os.path.join(inputData, img_))

                velocityMapArray = rasterInfo.ImageAsArray(bandNumber=1)
                snrMapArray = rasterInfo.ImageAsArray(bandNumber=2)

                # velocityMapFinal = WeightedMeanFilterV3(velocityMapArray=velocityMapArray, snrMapArray=snrMapArray,
                #                                       thresholdVelocity=thresholdVelovity, thresholdMedian=thresholdMedian, step=step)

                # velocityMapFinal = WeightedMeanFilterV3(velocityMapArray=velocityMapArray, snrMapArray=snrMapArray,
                #                                       thresholdVelocity=thresholdVelovity, thresholdMedian=thresholdMedian,
                #                                         step=step)

                velocityMapFinal = self.WeightedMeanFilterV4(velocityMapArray=velocityMapArray, snrMapArray=snrMapArray,
                                                             thresholdVelocity=thresholdVelovity,
                                                             thresholdMedian=thresholdMedian,
                                                             step=step, multiScaleFactor=multiScaleFactor)
                ### Save the smoothed velocity map as raster
                print("max=", velocityMapFinal[~np.isnan(velocityMapFinal)].max())
                print("min=", velocityMapFinal[~np.isnan(velocityMapFinal)].min())

                metaData = rasterInfo["MetaData"]
                metaData["MultiScaleFilterParams"] = "VThr_" + str(thresholdVelovity) + "_MedianThr_" + str(
                    thresholdMedian) + "_MultiscaleFactor_" + str(multiScaleFactor) + "_Step_" + str(step)

                if withSNR == False:

                    geoRT.WriteRaster(refRasterPath=os.path.join(inputData, img_),
                                      newRasterPath=os.path.join(outputData, img_),
                                      Listarrays=[velocityMapFinal], numberOfBands=1, metaData=metaData)
                else:
                    geoRT.WriteRaster(refRasterPath=os.path.join(inputData, img_),
                                      newRasterPath=os.path.join(outputData, img_),
                                      Listarrays=[velocityMapFinal, snrMapArray], numberOfBands=2,
                                      descriptions=["Multiscale Filtered Map", "SNR"], metaData=metaData)

                # break  # to be removed
            return
        else:
            print("Raster don't have SNR band")
            return

    def MultiScaleFilterWithSNR_(self, inputData, outputFolder, snrData, thresholdVelovity, thresholdMedian, step,
                                 multiScaleFactor=2, withSNR=False):
        """

        :param inputData:
        :param outputData:
        :param thresholdVelovity:
        :param thresholdMedian:
        :param step:
        :param multiScaleFactor:
        :param withSNR:
        :return:
        """

        # for img_ in imgs :
        for img_, snr_ in zip(inputData, snrData):
            # img_ = imgs[im]
            print("\n", img_)
            rasterInfo = geoRT.RasterInfo(inputRaster=img_)

            velocityMapArray = rasterInfo.ImageAsArray(bandNumber=1)
            snrMapArray = geoRT.RasterInfo(snr_).ImageAsArray(bandNumber=1)

            velocityMapFinal = self.WeightedMeanFilterV4(velocityMapArray=velocityMapArray, snrMapArray=snrMapArray,
                                                         thresholdVelocity=thresholdVelovity,
                                                         thresholdMedian=thresholdMedian,
                                                         step=step, multiScaleFactor=multiScaleFactor)
            ### Save the smoothed velocity map as raster
            print("max=", velocityMapFinal[~np.isnan(velocityMapFinal)].max())
            print("min=", velocityMapFinal[~np.isnan(velocityMapFinal)].min())

            metaData = rasterInfo["MetaData"]
            metaData["MultiScaleFilterParams"] = "VThr_" + str(thresholdVelovity) + "_MedianThr_" + str(
                thresholdMedian) + "_MultiscaleFactor_" + str(multiScaleFactor) + "_Step_" + str(step)

            if withSNR == False:

                geoRT.WriteRaster(refRasterPath=img_,
                                  newRasterPath=os.path.join(outputFolder, Path(img_).stem + "_filtered.tif"),
                                  Listarrays=[velocityMapFinal], numberOfBands=1, metaData=metaData)
            else:
                geoRT.WriteRaster(refRasterPath=img_,
                                  newRasterPath=os.path.join(outputFolder, Path(img_).stem + "_filtered.tif"),
                                  Listarrays=[velocityMapFinal, snrMapArray], numberOfBands=2,
                                  descriptions=["Multiscale Filtered Map", "SNR"], metaData=metaData)

            # break  # to be removed
        return
