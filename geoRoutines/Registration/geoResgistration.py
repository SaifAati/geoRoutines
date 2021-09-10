import sys, warnings

import geospatialroutine.georoutines as geoRT
import geospatialroutine.Routine as RT
import geospatialroutine.Plotting.Plotting_Routine as pltRT
import matplotlib.pyplot as plt
import os, sys
# import geospatialroutine.DEM_routine.DEM_profile as demRT
import geospatialroutine.Registration.Offset_correction as offCorrRT
import geospatialroutine.FilesCommandRoutine as fileRT
import numpy as np
from pathlib import Path
from p_tqdm import p_map
import geospatialroutine.FootPrint as fpRT
import geospatialroutine.Routine_Lyr as lyrRT
import gdal
import geospatialroutine.Filter.Filtering_Routine as filterRT


class cDEMRegistration:
    ## Private members
    iBaseArray = None
    iTargetArray = None
    iBaseDEMInfo = None
    iTargetDEMInfo = None
    intersectionCoords = 0
    Lon = []
    Lat = []
    oGSD = None
    oFolder = None
    iBaseCrop = None
    iTargetCrop = None
    workFolder = None
    doDemArray = None

    # constructor
    def __init__(self, iBaseDEMPath, iTargetDEMPath, oFolder=None, oGSD=None, extentWindow=[]):

        """
        Note: add the option to compute the correction using data in the window extent
        Args:
            iBaseDEMPath:
            iTargetDEMPath:
            oFolder:
            oGSD:
            extentWindow:
        """
        self.iBaseDEMPath = iBaseDEMPath
        self.iBaseDEMPath = self.iBaseDEMPath
        self.iTargetDEMPath = iTargetDEMPath
        self.iTargetDEMPath_UTM = self.iTargetDEMPath
        self.extentWindow = extentWindow
        self.oGSD = oGSD
        self.oFolder = oFolder
        if self.oFolder == None:
            self.workFolder = fileRT.CreateDirectory(directoryPath=os.path.dirname(self.iTargetDEMPath),
                                                     folderName=Path(self.iBaseDEMPath).stem + "_" +
                                                                Path(self.iTargetDEMPath).stem + "_registration",
                                                     cal="y")
        else:
            self.workFolder = fileRT.CreateDirectory(directoryPath=self.oFolder,
                                                     folderName=Path(self.iBaseDEMPath).stem + "_" +
                                                                Path(self.iTargetDEMPath).stem + "_registration",
                                                     cal="y")
            ## ## need to check how to convert meter to arcsecend
            ## for instance the average GSD is be computed

    def Decimalmod(self, value, param, precision=None):
        if precision == None:
            precision = 1e-5

        result = value % param

        if (np.abs(result) < precision) or (param - np.abs(result) < precision):
            result = 0
        return result

    def SetoGSD(self):
        """

        Returns:

        """
        if self.oGSD == None:
            xRes = np.mean([self.iTargetDEMInfo_UTM.pixelWidth, self.iTargetDEMInfo_UTM.pixelWidth])
            yRes = np.mean(
                [np.abs(self.iTargetDEMInfo_UTM.pixelHeight), np.abs(self.iTargetDEMInfo_UTM.pixelHeight)])
            self.oGSD = (xRes, yRes)
        return

    def GetOverlapAreaOfRasters(self, display=False):
        """
        Return the intersection area between a list of rasters
        :param rasterPathList: list of the raster paths
        :return: 0 if no intersection
                 geojson format of the intersection are, and write on the same directory of the images a temp.geojson
                 that contain the intersection are footprint
        """
        path = os.path.dirname(self.iBaseDEMPath)
        fpPathList = []
        fpTempFolder = fileRT.CreateDirectory(os.path.dirname(path), "Temp_SA_FP", cal="y")
        for img_ in [self.iBaseDEMPath, self.iTargetDEMPath]:
            extent, fpPath = fpRT.RasterFootprint(rasterPath=img_, writeFp=True,
                                                  savingPath=os.path.join(fpTempFolder, Path(img_).stem))
            fpPathList.append(fpPath)

        overlay = lyrRT.Intersection(fp1=fpPathList[0], fp2=fpPathList[1], dispaly=display)

        if overlay.intersection != 0:
            # # print(overlay.res_inter)
            tempIntersect = fpRT.WriteJson(features=overlay.res_inter,
                                           outputFile=os.path.join(fpTempFolder, "Temp"))
            intersection = overlay.intersection
            index = 2
            while intersection != 0 and index < len(fpPathList):

                overlay = lyrRT.Intersection(fp1=tempIntersect, fp2=fpPathList[index], dispaly=display)
                intersection = overlay.intersection
                index += 1
                if intersection != 0:
                    tempIntersect = fpRT.WriteJson(features=overlay.res_inter,
                                                   outputFile=os.path.join(fpTempFolder, "Temp"))

            if intersection == 0 and index != len(fpPathList) - 1:
                print("\nNo overlapping between images !!!")
                return 0
            else:
                from shapely.geometry import mapping
                interDic = mapping(overlay.res_inter["geometry"])
                self.intersectionCoords = interDic["features"][0]["geometry"]["coordinates"][0]
                print("\nFinal intersection:", self.intersectionCoords)

        else:
            print("\nNo overlapping between images !!!")

    def ComputeoMapGrid(self):
        self.oRes = self.oGSD
        if self.Decimalmod(self.upleftEW, self.oRes[0], self.oRes[0] / 1000) != 0:
            if (self.upleftEW - (self.upleftEW % self.oRes[0])) > self.upleftEW:
                self.oUpleftEW = self.upleftEW - (self.upleftEW % self.oRes[0])
            else:

                self.oUpleftEW = (self.upleftEW - (self.upleftEW % self.oRes[0]) + self.oRes[0])

        if self.Decimalmod(self.lowerEW, self.oRes[0], self.oRes[0] / 1000) != 0:
            if (self.lowerEW - (self.lowerEW % self.oRes[0])) < self.lowerEW:
                self.olowerEW = self.lowerEW - (self.lowerEW % self.oRes[0])
            else:
                self.olowerEW = (self.lowerEW - (self.lowerEW % self.oRes[1]) + self.oRes[1])

        if self.Decimalmod(self.upleftNS, self.oRes[1], self.oRes[1] / 1000) != 0:
            if (self.upleftNS - (self.upleftNS % self.oRes[1])) < self.upleftNS:
                self.oUpleftNS = self.upleftNS - (self.upleftNS % self.oRes[1])
            else:
                self.oUpleftNS = (self.upleftNS - (self.upleftNS % self.oRes[1]) + self.oRes[1])

        if self.Decimalmod(self.lowerNS, self.oRes[1], self.oRes[1] / 1000) != 0:
            if (self.lowerNS - (self.lowerNS % self.oRes[1])) > self.lowerNS:
                self.olowerNS = self.lowerNS - (self.lowerNS % self.oRes[1])
            else:
                self.olowerNS = (self.lowerNS - (self.lowerNS % self.oRes[1]) + self.oRes[1])

            return

    def CropSameGSDandExtent(self, vrt=True):
        """

        Args:
            vrt:

        Returns:

        """

        oFileList = []

        self.upleftEW = min(self.xMap)
        self.upleftNS = max(self.yMap)
        self.lowerNS = max(self.xMap)
        self.lowerEW = min(self.yMap)
        self.ComputeoMapGrid()
        for img_ in [self.iBaseDEMInfo_UTM, self.iTargetDEMInfo_UTM]:
            if self.oFolder == None:
                if vrt == True:
                    outputFile = os.path.join(self.workFolder,
                                              Path(img_.rasterPath).stem + "_crop_GSD_" + str(self.oGSD[0]) + ".vrt")
                    oFileList.append(outputFile)
                else:
                    outputFile = os.path.join(self.workFolder,
                                              Path(img_.rasterPath).stem + "_crop_GSD_" + str(self.oGSD[0]) + ".tif")
                    oFileList.append(outputFile)

            if vrt == True:
                gdal.Translate(destName=oFileList[-1],
                               srcDS=gdal.Open(img_.rasterPath),
                               # projWin=[self.Lon[0], self.Lat[1], self.Lon[2], self.Lat[0]],
                               projWin=[self.oUpleftEW, self.oUpleftNS, self.olowerNS, self.olowerEW],
                               xRes=self.oGSD[0], yRes=-1 * self.oGSD[1])

            else:
                format = "GTiff"
                gdal.Translate(destName=oFileList[-1],
                               srcDS=gdal.Open(img_.rasterPath),
                               # projWin=[self.Lon[0], self.Lat[1], self.Lon[2], self.Lat[0]],
                               projWin=[self.oUpleftEW, self.oUpleftNS, self.olowerNS, self.olowerEW],
                               xRes=self.oGSD[0], yRes=-1 * self.oGSD[1])

        self.iBaseCrop = oFileList[0]
        self.iTargetCrop = oFileList[1]
        return

    def ComputeDoD(self):
        baseDemCropInfo = geoRT.RasterInfo(self.iBaseCrop)
        targetDemCropInfo = geoRT.RasterInfo(self.iTargetCrop)
        baseDemCropArray = baseDemCropInfo.ImageAsArray()
        targetDemCropArray = targetDemCropInfo.ImageAsArray()

        baseDemCropArray_ma = np.ma.masked_where(baseDemCropArray == -32767, baseDemCropArray)
        targetDemCropArray_ma = np.ma.masked_where(targetDemCropArray == -32767, targetDemCropArray)
        if targetDemCropArray.shape != baseDemCropArray:
            warnings.warn("Inputs dont have the same shape! ==> Reshaping !!")
            print(baseDemCropArray_ma.shape, targetDemCropArray_ma.shape)
            if baseDemCropArray_ma.shape[0] > targetDemCropArray_ma.shape[0]:
                baseDemCropArray_ma = baseDemCropArray_ma[0:targetDemCropArray_ma.shape[0], :]
            if baseDemCropArray_ma.shape[0] < targetDemCropArray_ma.shape[0]:
                targetDemCropArray_ma = targetDemCropArray_ma[0:baseDemCropArray_ma.shape[0], :]
            if baseDemCropArray_ma.shape[1] > targetDemCropArray_ma.shape[1]:
                baseDemCropArray_ma = baseDemCropArray_ma[:, 0:targetDemCropArray_ma.shape[1]]
            if baseDemCropArray_ma.shape[1] < targetDemCropArray_ma.shape[1]:
                targetDemCropArray_ma = targetDemCropArray_ma[:, 0:baseDemCropArray_ma.shape[1]]

        self.doDemArray = baseDemCropArray_ma - targetDemCropArray_ma
        return

    def ComputeCorrection(self, noData=-32767, maskLargeValues=False, max_displacement=None, maskZeros=False,
                          showDiagnose=False, saveDiagnose=False):
        """"""

        # demRT.Plot_DEM(dtmPath=inputPath, vmin=-100, vmax=100, hillshde=False)
        # pltRT.PlotDistribution(inputArray=self.doDemArray, xlim=[])
        ## Mask large displacements
        self.doDemArray = np.ma.masked_where(self.doDemArray == noData, self.doDemArray)
        if maskLargeValues == True:
            self.doDemArray = filterRT.MaskLargeValues(inputArray=self.oDemArray, maxDisplacement=max_displacement)

        if maskZeros == True:
            self.doDemArray = filterRT.MaskZeros(inputArray=self.doDemArray)

        # print(inputDisp_array_before_correction)
        ### Mask nan values
        # if isinstance(inputDisp_arrayPix, np.ma.MaskedArray):
        self.doDemArray = np.ma.masked_invalid(self.doDemArray)
        # print(inputDisp_array_before_correction)
        ## Get statistical infos on the input
        val_mean = np.ma.mean(self.doDemArray)
        val_std = np.ma.std(self.doDemArray)

        val_min = val_mean - val_std
        val_max = val_mean + val_std
        print("===========================================")
        print("val_mean=%f,val_std=%f" % (val_mean, val_std))
        print("val_min=%f,val_max=%f" % (val_min, val_max))
        print("===========================================")

        # out_dict['input']['rmse_gauss'] = rmse_gauss
        # out_dict['input']['std_gauss'] = std_gauss
        # out_dict['input']['mean_gauss'] = mean_gauss
        # print("--- Statistic about displacement map after Gaussian filter=\n", out_dict)

        print("\n--- Fit a plan ")
        # return inputDisp_arrayPix
        # Fit plane
        _, _, corrParam = offCorrRT.FitPlaneIter(inputDispArray=self.doDemArray, iterations=10)
        print(corrParam)
        iTargetDEMArray = self.iTargetDEMInfo.ImageAsArray()
        iTargetDEMArray_ma = np.ma.masked_where(iTargetDEMArray == -32767, iTargetDEMArray)
        shape = iTargetDEMArray_ma.shape
        x_size = shape[1]
        y_size = shape[0]
        nx = np.linspace(-1, 1, x_size)
        ny = np.linspace(-1, 1, y_size)
        x, y = np.meshgrid(nx, ny)
        # Construct design matrix
        x_fl = x.flatten()
        y_fl = y.flatten()
        z_ones = np.ones([x.size, 1])
        X = np.hstack((np.reshape(x_fl, ([len(x_fl), 1])), np.reshape(y_fl, ([len(y_fl), 1])), z_ones))
        oValues = np.dot(X, corrParam)
        oCrrectionGrid = np.reshape(oValues, shape)

        # correction with plane == deramping
        oCrrectionGrid = oCrrectionGrid.astype(np.int32)
        iTargetDEMArray_ma = iTargetDEMArray_ma.astype(np.int32)
        oTargetDEMArray_ma = iTargetDEMArray_ma + oCrrectionGrid
        # pltRT.PlotDistribution(inputArray=inputDisp_array, xlim=[])
        if self.oFolder == None:
            geoRT.WriteRaster(oRasterPath=os.path.join(self.workFolder,
                                                       Path(self.iTargetDEMPath).stem + "_alig_r" + Path(
                                                           self.iBaseDEMPath).stem) + ".tif",
                              geoTransform=self.iTargetDEMInfo.geoTrans,
                              arrayList=[oTargetDEMArray_ma], epsg=self.iTargetDEMInfo.EPSG_Code, noData=-32767)
        else:
            geoRT.WriteRaster(oRasterPath=os.path.join(self.oFolder,
                                                       Path(self.iTargetDEMPath).stem + "_alig_r" + Path(
                                                           self.iBaseDEMPath).stem) + ".tif",
                              geoTransform=self.iTargetDEMInfo.geoTrans,
                              arrayList=[oTargetDEMArray_ma], epsg=self.iTargetDEMInfo.EPSG_Code, noData=-32767)

        if showDiagnose or saveDiagnose:
            print("Plotiing ....")
            cmap = 'seismic'
            # cmap = "rainbow"
            # cmap = 'gist_earth'
            # plot
            fig, axes = plt.subplots(ncols=3, figsize=(10, 3))
            vis_max = max(abs(val_min), abs(val_max))
            vis_min = -vis_max
            # vis_min =-0.5
            # vis_max = 0.5

            axes[0].imshow(oCrrectionGrid, cmap=cmap, vmin=vis_min, vmax=vis_max)
            axes[0].set_title('Fitted plane', fontsize=12)
            axes[0].tick_params(labelsize=12)

            vis_min = np.ma.mean(iTargetDEMArray_ma) - 2 * np.ma.std(iTargetDEMArray_ma)
            vis_max = np.ma.mean(iTargetDEMArray_ma) + 2 * np.ma.std(iTargetDEMArray_ma)
            axes[1].imshow(iTargetDEMArray_ma, cmap=cmap, vmin=vis_min, vmax=vis_max)
            axes[1].set_title('Input Image', fontsize=12)
            axes[1].tick_params(labelsize=12)

            vis_min = np.ma.mean(oTargetDEMArray_ma) - 2 * np.ma.std(oTargetDEMArray_ma)
            vis_max = np.ma.mean(oTargetDEMArray_ma) + 2 * np.ma.std(oTargetDEMArray_ma)
            im = axes[2].imshow(oTargetDEMArray_ma, cmap=cmap, vmin=vis_min, vmax=vis_max)
            axes[2].set_title('After plane correction', fontsize=12)
            axes[2].tick_params(labelsize=12)
            fig.subplots_adjust(right=0.91)

            # from mpl_toolkits.axes_grid1 import make_axes_locatable
            # divider = make_axes_locatable(axes[2])
            # cax = divider.append_axes("right", size="5%", pad=0.05)
            # cb = fig.colorbar(im, cax=cax)
            # cb.ax.tick_params(labelsize=12)
            # axes[0].axis("off")
            # axes[1].axis("off")
            # axes[2].axis("off")
            if showDiagnose:
                plt.show()
                # plt.draw()
                # plt.pause(1)
            # if saveDiagnose:
            #     fig.savefig(diagnoseSavingPath + "img_" + str(imgIndex) + ".png")
            plt.close()
        return

    def Registration(self):
        self.iBaseDEMInfo = geoRT.RasterInfo(self.iBaseDEMPath)
        self.iTargetDEMInfo = geoRT.RasterInfo(self.iTargetDEMPath)
        self.GetOverlapAreaOfRasters(display=False)
        if self.intersectionCoords != 0:
            for coord_ in self.intersectionCoords:
                self.Lon.append(coord_[0])
                self.Lat.append(coord_[1])
            print(self.Lat)
            print(self.Lon)
            utmEpsg = geoRT.ComputeEpsg(lon=self.Lon[0], lat=self.Lat[0])
            coord_UTM = geoRT.ConvCoordMap1ToMap2_Batch(X=self.Lat, Y=self.Lon, sourceEPSG=4326, targetEPSG=utmEpsg)
            self.xMap = coord_UTM[0]
            self.yMap = coord_UTM[1]
            if self.iBaseDEMInfo.EPSG_Code != utmEpsg:
                self.iBaseDEMPath_UTM = os.path.join(os.path.dirname(self.iBaseDEMInfo.rasterPath),
                                                     Path(self.iBaseDEMInfo.rasterPath).stem) + "_UTM.tif"
                gdal.Warp(self.iBaseDEMPath_UTM, self.iBaseDEMPath, dstSRS='EPSG:' + str(utmEpsg))

            if self.iTargetDEMInfo.EPSG_Code != utmEpsg:
                self.iTargetDEMPath_UTM = os.path.join(os.path.dirname(self.iTargetDEMInfo.rasterPath),
                                                       Path(self.iTargetDEMInfo.rasterPath).stem) + "_UTM.tif"
                gdal.Warp(self.iTargetDEMPath_UTM, self.iTargetDEMPath, dstSRS='EPSG:' + str(utmEpsg))
            self.iBaseDEMInfo_UTM = geoRT.RasterInfo(self.iBaseDEMPath_UTM)
            self.iTargetDEMInfo_UTM = geoRT.RasterInfo(self.iTargetDEMPath_UTM)
            self.SetoGSD()
            print(self.oGSD)
            self.CropSameGSDandExtent()
            self.ComputeDoD()
            self.ComputeCorrection()

        else:
            return

        return


def DerampingDem(inputPath, oFolder, iterations=20, showDiagnose=False,
                 saveDiagnose=False, diagnoseSavingPath="", maskLargeValues=False,
                 max_displacement=500):
    """

    Args:
        inputPath: 
        derampingPath: 
        iterations:
    Notes:
        - Add the correction to the raster metadata 
    """
    # inputPath = os.path.join(path, "3-PS2017_PS2019.tif")
    rasterInfo = geoRT.RasterInfo(inputPath)
    rasterArray = rasterInfo.ImageAsArray()

    rasterArray = np.ma.masked_where(rasterArray == -32767, rasterArray)
    # demRT.Plot_DEM(dtmPath=inputPath, vmin=-100, vmax=100, hillshde=False)
    # pltRT.PlotDistribution(inputArray=rasterArray, xlim=[])

    output,correction = offCorrRT.Deramping(
        displacement_field=inputPath,
        ref_path=inputPath,
        bandNumber=1,
        snr_filed=None,
        snr_threshold=0.333,
        iterations=iterations,
        zoom=1,
        showDiagnose=showDiagnose,
        saveDiagnose=saveDiagnose,
        diagnoseSavingPath=diagnoseSavingPath,
        imgIndex=None,
        maskLargeValues=maskLargeValues,
        max_displacement=max_displacement,
        maskZeros=False,

        xlim=[], ylim=[0, 3])

    metaData = rasterInfo.metaData
    oRasterPath = os.path.join(oFolder, Path(inputPath).stem + "_deramped.tif")
    RT.WriteRaster(refRasterPath=inputPath,
                   newRasterPath=oRasterPath,
                   Listarrays=[output],
                   numberOfBands=1, metaData=metaData)
    return oRasterPath, output, correction


def DoD(demPath1, demPath2, oFolder, oRes=10):
    dem1Info = geoRT.RasterInfo(demPath1)
    dem2Info = geoRT.RasterInfo(demPath2)

    print(dem1Info.rasterHeight, dem1Info.rasterHeight, dem1Info.pixelWidth, dem1Info.EPSG_Code)
    print(dem2Info.rasterHeight, dem2Info.rasterHeight, dem2Info.pixelWidth, dem2Info.EPSG_Code)
    demObj = cDEMRegistration(iBaseDEMPath=demPath1, iTargetDEMPath=demPath2, oGSD=(oRes, oRes))
    demObj.iBaseDEMInfo = dem1Info
    demObj.iTargetDEMInfo = dem2Info
    demObj.iBaseDEMPath_UTM = demPath1
    demObj.iTargetDEMPath_UTM = demPath2
    demObj.GetOverlapAreaOfRasters(display=False)
    if demObj.intersectionCoords != 0:
        for coord_ in demObj.intersectionCoords:
            demObj.Lon.append(coord_[0])
            demObj.Lat.append(coord_[1])
        #     print(self.Lat)
        #     print(self.Lon)
        utmEpsg = geoRT.ComputeEpsg(lon=demObj.Lon[0], lat=demObj.Lat[0])
        print(utmEpsg)
        coord_UTM = geoRT.ConvCoordMap1ToMap2_Batch(X=demObj.Lat, Y=demObj.Lon, sourceEPSG=4326, targetEPSG=utmEpsg)
        demObj.xMap = coord_UTM[0]
        demObj.yMap = coord_UTM[1]
        if demObj.iBaseDEMInfo.EPSG_Code != utmEpsg:
            demObj.iBaseDEMPath_UTM = os.path.join(os.path.dirname(demObj.iBaseDEMInfo.rasterPath),
                                                   Path(demObj.iBaseDEMInfo.rasterPath).stem) + "_UTM.tif"
            gdal.Warp(demObj.iBaseDEMPath_UTM, demObj.iBaseDEMPath, dstSRS='EPSG:' + str(utmEpsg))

        if demObj.iTargetDEMInfo.EPSG_Code != utmEpsg:
            demObj.iTargetDEMPath_UTM = os.path.join(os.path.dirname(demObj.iTargetDEMInfo.rasterPath),
                                                     Path(demObj.iTargetDEMInfo.rasterPath).stem) + "_UTM.tif"
            gdal.Warp(demObj.iTargetDEMPath_UTM, demObj.iTargetDEMPath, dstSRS='EPSG:' + str(utmEpsg))
        demObj.iBaseDEMInfo_UTM = geoRT.RasterInfo(demObj.iBaseDEMPath_UTM)
        demObj.iTargetDEMInfo_UTM = geoRT.RasterInfo(demObj.iTargetDEMPath_UTM)
        demObj.SetoGSD()
        print(demObj.oGSD)
        demObj.CropSameGSDandExtent()
        demObj.ComputeDoD()

        # plt.imshow(demObj.doDemArray,cmap="RdYlBu",vmax=50,vmin=-50)
        # print(demObj.doDemArray)
        oDoDArray = demObj.doDemArray.filled(fill_value=-32767)
        # print(oDoDArray)
        geoRT.WriteRaster(
            oRasterPath=os.path.join(oFolder, "DoD_" + Path(demPath1).stem + "_" + Path(demPath2).stem) + ".tif",
            geoTransform=[demObj.upleftEW, oRes, 0, demObj.upleftNS, 0, -oRes],
            arrayList=[oDoDArray], epsg=utmEpsg, noData=-32767)


if __name__ == '__main__':
    # pathDems = "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/PS2_DSM_Metashape"
    # # imgs =  ["DEM_2016_WGs84_onlyshisperSide.tif","DEM_2016.tif"]#  ["DEM_2017.tif", "DEM_2019.tif", "DEM_2020.tif",
    # # iDem = [os.path.join(pathDems, item) for item in imgs]
    #
    # iDem = fileRT.GetFilesBasedOnExtension(path=pathDems, disp=True)
    #
    # errorList = []
    # rDemPath = "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/0-AsterDEM.tif"
    # for img_ in iDem:
    #     try:
    #         demRegistration = cDEMRegistration(iBaseDEMPath=rDemPath,
    #                                        iTargetDEMPath=img_, oGSD=(10, 10))
    #         demRegistration.Registration()
    #     except:
    #
    #         errorList.append(Path(img_).stem)
    #
    # print("ErrorList:", errorList)

    # demPath1= "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/PS2_DSM_Metashape/0-AsterDEM_DEM_2019_3.75m_registration/DEM_2019_3.75m_alig_r0-AsterDEM.tif"
    #
    # demPath2List = [
    #
    #     "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/PS2_DSM_Metashape/0-AsterDEM_DEM_2020_8m_PS2_wgs84_registration/DEM_2020_8m_PS2_wgs84_alig_r0-AsterDEM.tif"
    # ]
    # for dem2 in demPath2List:
    #     DoD(demPath1=demPath1,
    #     demPath2=dem2,
    #     oFolder="/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/DoD_py")

    # # imgs = fileRT.GetFilesBasedOnExtension(path)
    # # DerampingDem(imgs[0],derampingPath)
    # # outputList = len(imgs)*[derampingPath]
    # # print(outputList)
    # # p_map(DerampingDem, imgs,outputList,num_cpus=56)
    #
    # # imgs = fileRT.GetFilesBasedOnExtension(derampingPath)
    # # for inputPath in imgs:
    # #     demRT.Plot_DEM(dtmPath=inputPath, vmin=-50, vmax=50, hillshde=False,saveFigure=True)
    # # break
    # path = "/home/cosicorr/0-WorkSpace/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D"
    # imgs = fileRT.GetFilesBasedOnExtension(path)
    # for inputPath in imgs:
    #     # rasterInfo = geoRT.RasterInfo(inputPath)
    #     # rasterArray = rasterInfo.ImageAsArray()
    #     # rasterArray = np.ma.masked_where(rasterArray == -32767, rasterArray)
    #     # pltRT.PlotDistribution(inputArray=rasterArray, xlim=[])
    #     demRT.Plot_DEM(dtmPath=inputPath, vmin=1000, vmax=7000, hillshde=True, saveFigure=False)
    #     # break

    demPath = "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/DoD_py/DoD_DEM_2017_3.8m_alig_r0-AsterDEM_DEM_2019_3.75m_alig_r0-AsterDEM.tif"
    demRT.Plot_DEM(dtmPath=demPath, vmin=-50, vmax=50, cmap="RdYlBu", hillshde=False)
