import os
import numpy as np
import matplotlib.pyplot as plt
import rasterio
import gdal
from rasterio.plot import show
from pathlib import Path
import geospatialroutine.Plotting.Plotting_Routine as PltRT
from matplotlib.ticker import FormatStrFormatter
import geospatialroutine.Routine_Lyr as LyrRT
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import geospatialroutine.FilesCommandRoutine as FileRT
import geospatialroutine.Filter.Filtering_Routine as filterRT

import geospatialroutine.Routine_Lyr as RTLyr
# import geospatialroutine.Routine as RT
import geospatialroutine.georoutines as geoRT


class ExtractProfile:
    def __init__(self, rasterPath, profilePath, **kwargs):
        """

        :param rasterPath:
        :param profilePath:
        :param kwargs: bandNumber offset=False, bandNumber=1, windwSize=0, maxVal=10
        """
        kwargs_ = {"bandNumber": 1, "windowSize": 3, "maxVal": 0, "offset": False, "center": 1000, "center_plus": 80,
                   "center_minus": 80, "fontSize": 14}
        if kwargs:
            kwargs_.update(kwargs)
        print("kwargs:", kwargs_)
        self.kwargs = kwargs_
        self.rasterPath = rasterPath
        self.profilePath = profilePath
        self.windowSize = self.kwargs["windowSize"]
        self.rasterInfo = geoRT.RasterInfo(self.rasterPath)
        self.bandNumber = self.kwargs["bandNumber"]
        self.epsg = self.rasterInfo.EPSG_Code
        self.profileValues = []
        self.coordinates = []
        self.maxVal = self.kwargs["maxVal"]
        self.Extract()
        if self.kwargs["offset"] == True:
            self.offset = None
            self.sigma = None
            self.centering = self.kwargs["center"]
            self.centering_plus = self.kwargs["center_plus"]
            self.centering_minus = self.kwargs["center_minus"]
            self.display = True
            self.OffsetCompute()

    def Extract(self):
        """

        """
        # RTLyr.GetVectorInfo(profile)
        shpRasterPath = RTLyr.RasterizeVector(shpInput=self.profilePath, refRaster=self.rasterPath,
                                              output=os.path.join(os.path.dirname(self.profilePath),
                                                                  Path(self.profilePath).stem + "_temp.tif"))
        rasterArray = self.rasterInfo.ImageAsArray(self.bandNumber)
        rasterArray = np.ma.masked_invalid(rasterArray)
        if self.maxVal != 0:
            rasterArray = np.ma.masked_where(rasterArray < -self.maxVal, rasterArray)
            rasterArray = np.ma.masked_where(rasterArray > self.maxVal, rasterArray)
            # rasterArray = filterRT.MaskLargeValues(inputArray=rasterArray,maxDisplacement=self.maxVal)
        shpInfo = geoRT.RasterInfo(shpRasterPath)
        shpArray = shpInfo.ImageAsArray() / 255
        shpArray = shpArray.astype(int)

        temp = np.ma.masked_invalid(rasterArray) * shpArray

        indexes = np.argwhere(temp != 0)
        print(indexes)
        w = self.windowSize
        for index in indexes:
            sub = rasterArray[index[0] - w:index[0] + w + 1, index[1] - w:index[1] + w + 1]
            meanVal = np.nanmean(sub)
            self.profileValues.append(meanVal)
            (xMap, yMap) = self.rasterInfo.Pixel2Map(x=index[0], y=index[1])
            self.coordinates.append((xMap, yMap))
            # print(self.coordinates)
        os.remove(shpRasterPath)

        return

    def cumDistances(self):
        """
        Note: we need to handle the fact that the epsg code is 4326 for distance computing 
        Returns:

        """
        distances = [0]
        for index, coord_ in enumerate(self.coordinates):

            if index == len(self.coordinates) - 1:
                break
            # print(coord_)
            distance = np.sqrt(((coord_[0] - self.coordinates[index + 1][0]) ** 2) + (
                    (coord_[1] - self.coordinates[index + 1][1]) ** 2))
            distances.append(distance)
        return np.cumsum(distances)

    def PlotProfile(self, ylim=[]):
        fig1, (ax1) = plt.subplots(1, 1)
        ax1.tick_params(direction='in', top=True, right=True, which="both", axis="both",
                        labelsize=self.kwargs["fontSize"])
        ax1.plot(self.cumDistances(), self.profileValues, linewidth="3")
        ax1.grid()
        ax1.minorticks_on()
        ax1.set_xlabel("Distance along profile [m]", fontsize=self.kwargs["fontSize"])
        ax1.set_ylabel("Displacement [m]", fontsize=self.kwargs["fontSize"])
        if any(ylim):
            ax1.set_ylim(ylim[0], ylim[-1])
        ax1.set_xlim(self.cumDistances()[0], self.cumDistances()[-1])
        plt.show()
        return

    def OffsetCompute(self):
        after_indexList = []
        before_indexList = []
        for index, val in enumerate(self.cumDistances()):
            if val >= self.centering + self.centering_plus:
                after_indexList.append(index)
            if val < self.centering - self.centering_minus:
                before_indexList.append(index)

        x_after = np.array([self.cumDistances()[i] for i in after_indexList])
        y_after = np.array([self.profileValues[i] for i in after_indexList])
        res_after = np.polyfit(x_after, y_after, 1, full=True)
        coef_after = res_after[0]
        sigma_after = np.sqrt(res_after[1] / (y_after.size - coef_after.size))
        # print(coef_after, sigma_after)
        x_before = np.array([self.cumDistances()[i] for i in before_indexList])
        y_before = np.array([self.profileValues[i] for i in before_indexList])
        res_before = np.polyfit(x_before, y_before, 1, full=True)
        coef_before = res_before[0]
        sigma_before = np.sqrt(res_before[1] / (y_before.size - coef_before.size))
        # print(coef_before, sigma_before)

        self.offset = np.abs(
            (coef_after[0] * x_after[0] + coef_after[1]) - (coef_before[0] * x_before[-1] + coef_before[1]))
        self.sigma = (sigma_before + sigma_after) / 2
        print(self.offset, self.sigma)
        ##### ploting ##################
        if self.display:
            fig1 = plt.figure(figsize=(10, 5))
            ax1 = fig1.add_subplot(1, 1, 1)
            ax1.plot(self.cumDistances(), self.profileValues, linewidth="3")
            ax1.plot(x_after, coef_after[0] * x_after + coef_after[1], color="k")
            ax1.plot(x_before, coef_before[0] * x_before + coef_before[1], color="k")
            ax1.tick_params(direction='in', top=True, right=True, which="both", axis="both",
                            labelsize=self.kwargs["fontSize"])

            ax1.axvline(self.centering + self.centering_plus, linewidth=1, linestyle="--", color='k')
            ax1.axvline(self.centering - self.centering_minus, linewidth=1, linestyle="--", color='k')
            ax1.grid()
            ax1.minorticks_on()
            ax1.set_xlabel("Distance along profile [m]", fontsize=self.kwargs["fontSize"])
            ax1.set_ylabel("Displacement [m]", fontsize=self.kwargs["fontSize"])
            ax1.set_title(Path(self.rasterPath).stem + "_band:" + str(self.bandNumber))

            _, xMax = ax1.get_xlim()
            yMin, yMax = ax1.get_ylim()
            text = '$\Delta=$' + "%.2f" % self.offset + "m\n"
            text += '$\sigma=$' + "%.2f" % self.sigma + "m"
            ax1.text(xMax - self.kwargs["center"], (yMin + yMax) / 2, text, bbox=dict(facecolor='white', alpha=0.5),
                     fontsize=self.kwargs["fontSize"])
            ax1.set_xlim(self.cumDistances()[0], self.cumDistances()[-1])
            # plt.savefig(os.path.join(
            #     "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/SVG_Figures_ICA_PCA/Dz/Reconstruct_4_7_ICA/Temp",
            #     Path(self.rasterPath).stem + "_band:" + str(self.bandNumber) + ".svg"))
            plt.show()
        return


def StackingProfileFromCosiCorr():
    path = "//media/cosicorr/storage/Saif/1-Ridgecrest/3D_approch_Ridgecrest/DEMS_correlation/FigureForPaper/Profiles/Profile1_"
    # pathRaster = "/media/cosicorr/storage/Saif/1-Ridgecrest/3D_approch_Ridgecrest/DEMS_correlation/FigureForPaper/"
    # profile1 = os.path.join(path, "StackProfile/Profile1.geojson")
    # profile2 = os.path.join(path, "StackProfile/Profile2.geojson")
    profile1Txt = os.path.join(path, "Profile.txt")
    # profile2Txt = os.path.join(path, "StackProfile//Profile2.txt")

    data1 = np.loadtxt(profile1Txt, comments=";")
    # print(data1)

    data1_subset = data1[16:, :]
    ew_offset = "%.2f" % (data1[5, 0])
    ew_sigma = "%.2f" % (data1[5, 1])
    ns_offset = "%.2f" % (data1[10, 0])
    ns_sigma = "%.2f" % (data1[10, 1])
    print(data1_subset)
    print(ew_offset, ew_sigma, ns_offset, ns_sigma)
    centering = 350

    fig1 = plt.figure(figsize=(10, 5))
    ax1 = fig1.add_subplot(111)
    ax1.plot(data1_subset[:, 0] - centering, data1_subset[:, 1], label="E/W", color="k", linewidth=3)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    # plt.xlabel("",fontsize=14)
    ax1.set_ylabel("Displacement [m]", fontsize=16)
    ax1.set_xlim(-350, 350)
    # plt.ylim(-1.5, 1.5)

    offset_txt = r"$\Delta$ [m]=" + str(ew_offset) + "\n" + r"$\sigma$ [m]=" + str(ew_sigma)
    plt.text(ax1.get_xlim()[1] - 500, 0.5, offset_txt, color="k", bbox=dict(facecolor='white', alpha=0.8),
             fontsize=14)
    plt.legend(fontsize=16)

    plt.grid()

    fig2 = plt.figure(figsize=(10, 5))
    ax2 = fig2.add_subplot(111)
    ax2.plot(data1_subset[:, 0] - centering, data1_subset[:, 2], linewidth=3, color="k", label="N/S")
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    # plt.xlabel("",fontsize=14)
    plt.ylabel("Displacement [m]", fontsize=16)
    plt.legend(fontsize=14)

    ax2.set_xlim(-350, 350)
    # plt.ylim(-1.5, 1.5)
    offset_txt = r"$\Delta$ [m]=" + str(ns_offset) + "\n" + r"$\sigma$ [m]=" + str(ns_sigma)
    plt.text(ax2.get_xlim()[1] - 500, 0.5, offset_txt, color="k", bbox=dict(facecolor='white', alpha=0.8),
             fontsize=14)
    plt.grid()

    kwargs = {'format': 'GTiff', "outputType": gdal.GDT_Float32, "workingType": gdal.GDT_Float32}
    # gdal.Warp(os.path.join(path, "temp.tif"), pathRaster, dstSRS='EPSG:4326', **kwargs)

    # fig3, (ax) = plt.subplots(1,1 , figsize=(18, 18), dpi=100)
    #
    # dest = rasterio.open(os.path.join(path, "temp.tif"))
    # # dest = rasterio.open(pathRaster)
    #
    # im = show(dest.read(2), ax=ax, transform=dest.transform, cmap='Spectral', title="N/S", vmin=-3, vmax=3)
    # fpDF1 = LyrRT.ReadVector(profile1)
    # fpDF2 = LyrRT.ReadVector(profile2)
    #
    # fpDF1.plot(ax=ax, facecolor='none',edgecolor= "k",  linewidth=3.0)
    # fpDF2.plot(ax=ax, facecolor='none',edgecolor= "b",  linewidth=3.0)
    #
    # lon_formatter = LongitudeFormatter(zero_direction_label=True,number_format='.2f')
    # lat_formatter = LatitudeFormatter(number_format='.2f')
    # ax.xaxis.set_major_formatter(lon_formatter)
    # ax.yaxis.set_major_formatter(lat_formatter)
    # plt.xticks(rotation=45,fontsize=14)
    # # plt.xlabel("Lon",fontsize=14)
    # # plt.ylabel("Lat",fontsize=14)
    # plt.yticks(fontsize=14)

    plt.show()


def Profile_forUncertainty():
    path = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/6p4_3D/Correlation"
    profile1 = os.path.join(path, "Set1_DZ_Profile1 .txt")
    profile2 = os.path.join(path, "Set2_DZ_Profile1 .txt")
    profile3 = os.path.join(path, "Set3_DZ_Profile1 .txt")
    profile4 = os.path.join(path, "Set4_DZ_Profile1 .txt")
    data1 = np.loadtxt(profile1)
    data2 = np.loadtxt(profile2)
    data3 = np.loadtxt(profile3)
    data4 = np.loadtxt(profile4)

    fig1 = plt.figure(figsize=(10, 5))
    plt.plot(data1[:, 0], data1[:, 1], label="Set1_Dz", color="k", linewidth=2)
    plt.plot(data2[:, 0], data2[:, 1], label="Set2_Dz", color="g", linewidth=2)
    plt.plot(data3[:, 0], data3[:, 1], label="Set3_Dz", color="r", linewidth=2)
    plt.plot(data4[:, 0], data4[:, 1], label="Set4_Dz", color="b", linewidth=2)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    # plt.xlabel("",fontsize=14)
    plt.ylabel("Displacement(m)", fontsize=16)
    plt.legend(fontsize=16)
    plt.grid()
    plt.show()


def Plot_Profile_Qgis(xLabel="Distance along profile (m)", yLabel="Displacement (m)", fontsize=18, linewidth=3):
    # """
    ## Profile 1
    path = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/CompareDiffrentResiltsSets/Mw7p1_Lidar/ProfileCC'/Dz"
    # file1 = os.path.join(path, "ProfileCC'_Couple33_PreDEM.txt")
    # file2 = os.path.join(path, "ProfileCC'_dDSM.txt")
    # file3= os.path.join(path,  "Dz_ProfileCC'_Couple33_SRTM.txt")
    # profileList = [file1, file2]#,file3]
    profileList = FileRT.GetFilesBasedOnExtension(path, "*.txt")
    offset = 900
    # colors = ["chocolate","olive","teal"]
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)
    colors = ["gray", "blue", "aqua", "purple", "red", "orange", "peru", "green"]
    for index, profile_ in enumerate(profileList):
        data = np.loadtxt(profile_)
        # ax.plot(data[:, 0], data[:, 1], ".",label=Path(profile_).stem, linewidth=linewidth)
        label_ = Path(profile_).stem
        ax.plot(data[:, 0] - offset, data[:, 1], label=label_.split("_")[2], linewidth=linewidth, color=colors[index])

    ax.grid("on")
    # ax.set_xlim(-800, 800)
    # ax.set_ylim(-1, 2)
    ax.set_xlabel(xLabel, fontsize=fontsize)
    ax.set_ylabel(yLabel, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.show()
    return


############################################### Examples ###############################################################
def PlotProfilesDZ():
    profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileAA'.kml"
    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/Profile_Shps/ProfileBB'.kml"
    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileBB'.kml"
    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileCC'.kml"

    rasterPath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/Dz/ICA/Gaussian/ICA_reconst_4_7/ICA/DZ_reconstruct_4_7_ICA_Gaussina_Reconstruct_PC1_2"
    # rasterPath = "//home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/Dz/allDz.tif"
    rasterInfo = RT.GetRasterInfo(rasterPath, True)
    nbBands = rasterInfo["NbBands"]
    profilObj = []
    labels = []
    ############### based oo a single raster with multi band ####################################
    for band_ in range(nbBands):
        profile = ExtractProfile(rasterPath=rasterPath,
                                 profilePath=profilePath, bandNumber=band_ + 1, offset=True, windwSize=3)
        profilObj.append(profile)
        # profile.PlotProfile()
        labels.append(rasterInfo["BandInfo"][band_])

    referenceICA = "//home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/Dz/ICA/Gaussian/ICA_reconst_4_7/ICA/DZ_reconstruct_4_7_ICA_Gaussina_Reconstruct_PC1_2_Stack_median.tif"
    profilObj.append(
        ExtractProfile(rasterPath=referenceICA, profilePath=profilePath, bandNumber=1, offset=True, windwSize=3))
    labels.append("Reference_ICA")

    diffRastePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Compare_DSM_diff_and_couple33_PreDEM/DSM_diff_subset_.tif"
    diffRasterInfo = RT.GetRasterInfo(diffRastePath, True)
    profilObj.append(
        ExtractProfile(rasterPath=diffRastePath, profilePath=profilePath, bandNumber=1, offset=True, windwSize=3))
    labels.append("Diff_DSMs")

    fig1, (ax1) = plt.subplots(1, 1)
    centering = 600
    for index, profileObj_ in enumerate(profilObj):
        ax1.tick_params(direction='in', top=True, right=True, which="both", axis="both", labelsize=14)
        ax1.plot(profileObj_.cumDistances() - centering, profileObj_.profileValues, linewidth="3",
                 label=labels[index])
    ax1.grid()
    ax1.minorticks_on()
    ax1.set_xlabel("Distance along profile [m]", fontsize=14)
    ax1.set_ylabel("Displacement [m]", fontsize=14)
    # ax1.set_ylim(-1.8, 1)
    # ax1.set_xlim(self.cumDistances()[0], self.cumDistances()[-1])
    plt.legend()
    plt.show()
    return


def PlotProfilesEW():
    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileAA'.kml"
    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/Profile_Shps/ProfileBB'.kml"
    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileBB'.kml"

    profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileCC'.kml"
    rasterPath = "//home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/EW/allEW.tif"
    # rasterPath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/EW/ICA/allEW_ICA_Reconstruct_PC1&PC3"
    rasterInfo = RT.GetRasterInfo(rasterPath, True)
    nbBands = rasterInfo["NbBands"]
    profilObj = []
    labels = []
    ############### based oo a single raster with multi band ####################################
    for band_ in range(nbBands):
        profile = ExtractProfile(rasterPath=rasterPath,
                                 profilePath=profilePath, bandNumber=band_ + 1, offset=True)
        profilObj.append(profile)
        # profile.PlotProfile()
        labels.append(rasterInfo["BandInfo"][band_])

    # diffRastePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/Compare_DSM_diff_and_couple33_PreDEM/DSM_diff_subset_.tif"
    # diffRasterInfo = RT.GetRasterInfo(diffRastePath,True)
    # profilObj.append(ExtractProfile(rasterPath=diffRastePath, profilePath=profilePath, bandNumber= 1))
    # labels.append("Diff_DSMs")

    referenceICA = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/EW/ICA/allEW_ICA_Reconstruct_PC1&PC3_Stack_median.tif"
    profilObj.append(ExtractProfile(rasterPath=referenceICA, profilePath=profilePath, bandNumber=1, offset=True))
    labels.append("Reference_EW_ICA")
    centering = 780
    fig1, (ax1) = plt.subplots(1, 1)
    for index, profileObj_ in enumerate(profilObj):
        ax1.tick_params(direction='in', top=True, right=True, which="both", axis="both", labelsize=14)
        ax1.plot(profileObj_.cumDistances() - centering, profileObj_.profileValues, linewidth="3",
                 label=labels[index])
    ax1.grid()
    ax1.minorticks_on()
    ax1.set_xlabel("Distance along profile [m]", fontsize=14)
    ax1.set_ylabel("Displacement [m]", fontsize=14)

    # ax1.set_ylim(-2.5, 1)
    # ax1.set_xlim(self.cumDistances()[0], self.cumDistances()[-1])
    plt.legend()
    plt.show()
    # """
    return


def PlotProfilesNS():
    profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileAA'.kml"
    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/Profile_Shps/ProfileBB'.kml"
    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileBB'.kml"

    # profilePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Profile_Shps/_ProfileCC'.kml"

    rasterPath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/NS/ICA/allNS_ICA_Reconstruct_PC1_PC2"
    # rasterPath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/NS/allNS.tif"
    rasterInfo = RT.GetRasterInfo(rasterPath, True)
    nbBands = rasterInfo["NbBands"]
    profilObj = []
    labels = []
    ############### based oo a single raster with multi band ####################################
    for band_ in range(nbBands):
        profile = ExtractProfile(rasterPath=rasterPath,
                                 profilePath=profilePath, bandNumber=band_ + 1, offset=True)
        profilObj.append(profile)
        # profile.PlotProfile()
        labels.append(rasterInfo["BandInfo"][band_])

    # diffRastePath = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/Compare_DSM_diff_and_couple33_PreDEM/DSM_diff_subset_.tif"
    # diffRasterInfo = RT.GetRasterInfo(diffRastePath,True)
    # profilObj.append(ExtractProfile(rasterPath=diffRastePath, profilePath=profilePath, bandNumber= 1))
    # labels.append("Diff_DSMs")
    #
    referenceICA = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/NS/ICA/allNS_ICA_Reconstruct_PC1_PC2_Stack_median.tif"
    profilObj.append(ExtractProfile(rasterPath=referenceICA, profilePath=profilePath, bandNumber=1, offset=True))
    labels.append("Reference_NS_ICA")

    centering = 600

    fig1, (ax1) = plt.subplots(1, 1)
    for index, profileObj_ in enumerate(profilObj):
        ax1.tick_params(direction='in', top=True, right=True, which="both", axis="both", labelsize=14)
        ax1.plot(profileObj_.cumDistances() - centering, profileObj_.profileValues, linewidth="3",
                 label=labels[index])
    ax1.grid()
    ax1.minorticks_on()
    ax1.set_xlabel("Distance along profile [m]", fontsize=14)
    ax1.set_ylabel("Displacement [m]", fontsize=14)
    # ax1.set_ylim(-2.5, 1)
    # ax1.set_xlim(self.cumDistances()[0], self.cumDistances()[-1])
    plt.legend()
    plt.show()
    return


if __name__ == '__main__':
    # StackingProfileFromCosiCorr()
    # Plot_Profile_Qgis()
    # PlotProfilesDZ()

    path = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Sets/Couple_12/1-CCD_correction_classicMethod/NoCCDCorr/Correlation_sameExtent"
    profilePath = os.path.join(path, "CCD_Misalignement_Profile.kml")
    rasterPath = os.path.join(path, "16SEP08185521_NoCCDCorr_VS_18JUN16213928_NoCCD_corr_detrended_subset.tif")
    profile = ExtractProfile(rasterPath=rasterPath,
                             profilePath=profilePath, bandNumber=2, offset=False)
    profile.PlotProfile()
