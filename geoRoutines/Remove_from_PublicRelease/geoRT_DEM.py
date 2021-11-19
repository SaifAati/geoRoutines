import os, pandas, rasterio
import sys

import earthpy.spatial as es
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from rasterio.plot import show, show_hist
from scipy import stats
from pathlib import Path

from geospatialroutine.georoutines import GetWindow
import geospatialroutine.Plotting.Plotting_Routine as pltRT
import geospatialroutine.FilesCommandRoutine as FileRT
import geospatialroutine.georoutines as geoRT

import gdal
import matplotlib.transforms as mtransforms


def PlotProfileFromCSV(csvPath, xLabel="Distance along profile (km)", yLabel="Elevation (m)", nbProfile=2,
                       xlim=[-30, 30], nbins=50, linewidth=2, fontsize=14):
    """

    :param csvPath:
    :param nbProfile:   2 if the csv has two 2 profiles to be displayed
                        1 if the CSV has 1 profile to be displayed
    :return:
    """

    data = pandas.read_csv(csvPath)
    print("-----------", os.path.basename(csvPath))
    print(data)
    print(data.keys())

    fig = plt.figure(figsize=(6.5, 5))
    ax = fig.add_subplot(1, 1, 1)
    disKm = data['Distance (Total)'] / 1000

    if nbProfile == 2:
        # disKm = data[list(data.keys())[7]] / 1000
        prof1 = ax.plot(disKm, data[list(data.keys())[4]], label=data.keys()[4], linewidth=linewidth)
        prof2 = ax.plot(disKm, data[list(data.keys())[5]], label=data.keys()[5], linewidth=1.5)

        # prof3 = ax.plot(disKm, data[list(data.keys())[5]], label="SRTM DEM",linewidth=3,linestyle="--",color="blue")
        # prof3 = ax.plot(disKm, data[list(data.keys())[4]],label="DEM 2019",linewidth=3,linestyle="-",color="red")

    if nbProfile == 3:
        # disKm = data[list(data.keys())[11]] / 1000
        prof1 = ax.plot(disKm, data[list(data.keys())[4]], label=data.keys()[4], linewidth=linewidth)
        prof2 = ax.plot(disKm, data[list(data.keys())[5]], label=data.keys()[5], linewidth=linewidth)
        prof3 = ax.plot(disKm, data[list(data.keys())[6]], label=data.keys()[6], linewidth=linewidth)

    if nbProfile == 4:
        # disKm = data[list(data.keys())[11]] / 1000
        prof1 = ax.plot(disKm, data[list(data.keys())[4]], label=data.keys()[4], linewidth=linewidth)
        prof2 = ax.plot(disKm, data[list(data.keys())[5]], label=data.keys()[5], linewidth=linewidth)
        prof3 = ax.plot(disKm, data[list(data.keys())[6]], label=data.keys()[6], linewidth=linewidth)
        prof4 = ax.plot(disKm, data[list(data.keys())[7]], label=data.keys()[7], linewidth=linewidth)

        # prof3 = ax.plot(disKm, data[list(data.keys())[5]], label="SRTM DEM",linewidth=3,linestyle="--",color="blue")
        # prof3 = ax.plot(disKm, data[list(data.keys())[4]],label="DEM 2019",linewidth=3,linestyle="-",color="red")
    if nbProfile == 1:
        # print(data['Distance (Total)'])
        # disKm = data['Distance (Total)'] /1000 #data[list(data.keys())[7]] / 1000
        yData = data[list(data.keys())[3]]
        prof4 = ax.plot(disKm, yData, linewidth=linewidth, color="blue")
        print(np.mean(yData), np.std(yData))
        PlotDistribution(inputArray=yData, title="Elevation Difference", xLabel=yLabel, xlim=xlim,
                         nbins=nbins)

    ax.legend(fontsize=fontsize)
    ax.set_xlabel(xLabel, fontsize=fontsize)
    ax.set_ylabel(yLabel, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # ax.set_ylabel("Elevation difference (m)")
    # ax.set_title("Elevation  profile \nalong center flow line \nbetween GE WV DEM and SRTM DEM")
    # plt.show()

    # coordDistance2Plot.append(flowLineCoordPix[index_])
    #
    # # coordDistance2Plot = np.asarray(coordDistance2Plot)
    ax.grid("on")
    plt.show()

    return


def Plot_DEM(dtmPath, vmin=None, vmax=None, title=None, hillshde=True, saveFigure=False, cmap=None,
             applyGeoTransfrom=False, windowGeo=[], axes=True):
    """
    
    Args:
        dtmPath: 
        vmin: 
        vmax: 
        title: 
        hillshde: 
        saveFigure: 
        cmap: 
        applyGeoTransfrom: 
        windowGeo: 
        axes: 

    Returns:

    """
    # cmap = "hsv"
    # cmap = "twilight"
    # cmap = "twilight_shifted"
    # cmap ="seismic"
    # cmap = "bwr"
    # cmap = "RdBu"
    # cmap = "RdYlBu"

    # cmap = "gist_rainbow_r"

    if cmap == None:
        cmap = "RdYlBu"
        cmap = "terrain"
    # cmap = "gist_rainbow_r"
    fig = plt.figure(figsize=(10, 5))
    ax = fig.gca()
    if not title:
        title = Path(dtmPath).stem
    # with rasterio.open(dtmPath) as src:
    #     array = src.read(1)
    dtmInfo = geoRT.RasterInfo(dtmPath)
    array = geoRT.RasterInfo(dtmPath).ImageAsArray()
    xLim, yLim = dtmInfo.Pixel2Map_Batch(X=[0, dtmInfo.rasterWidth], Y=[dtmInfo.rasterHeight, 0])
    elevation = np.copy(array)
    if windowGeo:
        xLim = [windowGeo[0], windowGeo[1]]
        yLim = [windowGeo[2], windowGeo[3]]

    elevation_ma = np.ma.masked_where(elevation <= -32767, elevation, copy=True)
    print("xLim:{}".format(xLim))
    print("yLim:{}".format(yLim))
    # print(np.max(elevation), np.min(elevation))
    # elevation = np.ma.masked_array(elevation, mask=mask2)
    # elevation = np.ma.masked_array(elevation, mask=mask1)
    # elevation[elevation==-32767] = np.nan

    trans = mtransforms.Affine2D(matrix=geoRT.GeoTransfomAsArray(geoTrans=dtmInfo.geoTrans))
    transData = trans + ax.transData
    if vmin and vmax:
        if applyGeoTransfrom:
            im = ax.imshow(elevation_ma, interpolation=None, cmap=cmap, vmin=vmin, vmax=vmax, clip_on=True,
                           transform=transData)
        else:
            if windowGeo:
                from geospatialroutine.georoutines import GetWindow
                win = GetWindow(iRasterPath=dtmPath, windowGeo=windowGeo)
                elevation_ma = np.copy(elevation_ma[win[1]:win[1] + win[-1], win[0]:win[0] + win[2]])
                elevation_ma = np.ma.masked_where(elevation_ma <= -32767, elevation_ma, copy=True)
                im = ax.imshow(elevation_ma, interpolation=None, cmap=cmap, vmin=vmin, vmax=vmax)

            else:
                im = ax.imshow(elevation_ma, interpolation=None, cmap=cmap, vmin=vmin, vmax=vmax)

    else:
        if applyGeoTransfrom:
            im = ax.imshow(elevation_ma, interpolation=None, cmap=cmap, vmin=vmin, vmax=vmax, clip_on=True,
                           transform=transData)
        else:
            if windowGeo:
                from geospatialroutine.georoutines import GetWindow
                win = GetWindow(iRasterPath=dtmPath, windowGeo=windowGeo)
                elevation_ma = np.copy(elevation_ma[win[1]:win[1] + win[-1], win[0]:win[0] + win[2]])
                im = ax.imshow(elevation_ma, interpolation=None, cmap=cmap)

            else:
                im = ax.imshow(elevation_ma, interpolation=None, cmap=cmap)
    if hillshde:
        hillshade = es.hillshade(elevation_ma, azimuth=90, altitude=10)
        if applyGeoTransfrom:
            ax.imshow(hillshade, cmap="Greys", clip_on=True, transform=transData, alpha=0.5)

        else:

            ax.imshow(hillshade, cmap="Greys", alpha=0.5)

    if applyGeoTransfrom:
        ax.set_xlim(xLim)
        ax.set_ylim(yLim)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # pltRT.Colorbar_2(im)
    pltRT.ColorBar_(ax=ax, mapobj=im, cmap=cmap, vmin=-100, vmax=100, size=14, label="Elevation[m]")
    ax.set_title(title)
    # fig2 = plt.figure()
    # show_hist(elevation, bins=50, lw=0.0, stacked=False, alpha=0.3, histtype='stepfilled', title="Histogram")
    if axes:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        if applyGeoTransfrom and dtmInfo.EPSG_Code == 4326:
            ax.set_xlabel("Lon [$^\circ$]", fontsize=14)
            ax.set_ylabel("Lat [$^\circ$]", fontsize=14)
        xVals = ax.get_xticks()
        yVals = ax.get_yticks()
        ax.set_xticks(xVals[1::3])
        ax.set_yticks(yVals[1::3])
    else:
        ax.axis('off')
    if saveFigure:
        fig.savefig(os.path.join(os.path.dirname(dtmPath), Path(dtmPath).stem + ".svg"), dpi=400)

    plt.show()
    return


def Plot_DEM_withMask(dtmPath, maskPath, basemap, vmin=None, vmax=None, title=None):
    cmap = "hsv"
    # cmap = "twilight"
    # cmap = "twilight_shifted"
    # cmap ="seismic"
    # cmap = "bwr"
    # cmap = "RdBu"
    cmap = "RdYlBu"
    # cmap = "gist_rainbow_r"
    # cmap = "terrain"
    fig = plt.figure()
    ax = fig.gca()
    if not title:
        title = os.path.basename(dtmPath)
    with rio.open(maskPath) as masksrc:
        mask_ = masksrc.read(1)
        print(mask_)
        # plt.imshow(mask_)
        # plt.show()
        mask = mask_ == 0
    with rio.open(dtmPath) as src:
        elevation = src.read(1)
        # Set masked values to np.nan
        # elevation[elevation < 0.0] = np.nan
        mask2 = elevation == 255
        # mask2 = elevation == 0
        # mask =elevation ==0
        # mask = elevation <= 0
        # mask2 = elevation> 65500
        elevation = np.ma.masked_array(elevation, mask=mask)
        elevation = np.ma.masked_array(elevation, mask=mask2)
        # elevation = np.ma.masked_array(elevation, mask=mask1)
        # elevation[elevation==-32767] = np.nan
    with rio.open(basemap) as src_basemap:
        basemap_ = src_basemap.read(1)

    print(elevation[0, 0])
    hillshade = es.hillshade(elevation, azimuth=90, altitude=10)
    hillshade_basemap = es.hillshade(basemap_, azimuth=360, altitude=16)
    # print(hillshade)
    # ep.plot_bands(elevation, ax=ax, cmap="terrain", title="HH DEM", )
    # ax.imshow(hillshade, cmap="Greys", alpha=0.5)

    # ep.plot_bands(elevation, ax=ax, cmap="terrain", title="HH DEM", )
    if vmin and vmax:
        im = ax.imshow(elevation, cmap=cmap, vmin=vmin, vmax=vmax)
        # ep.plot_bands(elevation, ax=ax, cmap="terrain", title="HH DEM", )

        show(elevation, ax=ax, cmap=cmap, title=title, transform=src.transform, vmin=vmin, vmax=vmax)

    else:
        im = ax.imshow(elevation, cmap=cmap)
        # ep.plot_bands(elevation, ax=ax, cmap="terrain", title="HH DEM", )
        show(elevation, ax=ax, cmap=cmap, title=title, transform=src.transform)
    # show(hillshade, ax=ax, cmap="Greys", alpha=0.5, transform=src.transform)
    show(hillshade_basemap, ax=ax, cmap="Greys", alpha=0.6, title=title, transform=src_basemap.transform, vmin=vmin,
         vmax=vmax)
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    colorbar(im)
    plt.show()
    # fig2 = plt.figure()
    # show_hist(elevation, bins=50, lw=0.0, stacked=False, alpha=0.3, histtype='stepfilled', title="Histogram")

    return


def Plot_ListDEM(dtmList, vmin=2000, vmax=6000, title=None):
    cmap = "hsv"
    cmap = "twilight"
    cmap = "twilight_shifted"
    # cmap ="seismic"
    # cmap = "bwr"
    # cmap = "RdBu"
    # cmap = "RdYlBu"
    cmap = "gist_rainbow_r"
    elevationList = []
    hillshadeList = []
    srcList = []
    for dtm in dtmList:
        with rio.open(dtm) as src:
            elevation = src.read(1)
            # Set masked values to np.nan
            # elevation[elevation < 0.0] = np.nan
            # mask1 = elevation == -32767
            mask = elevation <= 0.1
            mask2 = elevation >= 65535
            elevation = np.ma.masked_array(elevation, mask=mask)
            elevation = np.ma.masked_array(elevation, mask=mask2)
            # elevation[elevation==-32767] = np.nan
            hillshade = es.hillshade(elevation, azimuth=90, altitude=10)
            elevationList.append(elevation)
            hillshadeList.append(hillshade)
            srcList.append(src)
            print(elevation[0, 0], hillshade[0, 0])

    # #print(hillshade)
    # #ep.plot_bands(elevation, ax=ax, cmap="terrain", title="HH DEM", )
    # #ax.imshow(hillshade, cmap="Greys", alpha=0.5)
    # #ep.plot_bands(elevation, ax=ax, cmap="terrain", title="HH DEM", )
    for elevation, hillshade, src in zip(elevationList, hillshadeList, srcList):
        fig = plt.figure()
        ax = fig.gca()
        if vmin and vmax:
            im = ax.imshow(elevation, cmap=cmap, vmin=vmin, vmax=vmax)
            # ep.plot_bands(elevation, ax=ax, cmap="terrain", title="HH DEM", )
            show(elevation, ax=ax, cmap=cmap, title=title, transform=src.transform, vmin=vmin, vmax=vmax)
        else:
            im = ax.imshow(elevation, cmap=cmap)
            # ep.plot_bands(elevation, ax=ax, cmap="terrain", title="HH DEM", )
            show(elevation, ax=ax, cmap=cmap, title=title, transform=src.transform)
        print(src.transform)
        show(hillshade, ax=ax, cmap="Greys", alpha=0.5, transform=src.transform)
        colorbar(im)

        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)

    plt.show()
    # fig2 = plt.figure()
    # show_hist(elevation, bins=50, lw=0.0, stacked=False, alpha=0.3, histtype='stepfilled', title="Histogram")

    return


def Plot_DoD(dodPath, vmin=-50, vmax=50, cmap="RdYlBu", title=None, basemap=None, basemapHillshade=False,
             extent=[], maskShapeFile=None, noData=-32767, saveFig=False, plotDistribution=True, maxDisplacement=50):
    """
    Extent: same extent as qgis
    """
    rasterInfo = geoRT.RasterInfo(dodPath)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.gca()
    if not title:
        title = Path(dodPath).stem
    array = rasterInfo.ImageAsArray()


    if maskShapeFile:
        import geospatialroutine.Routine_Lyr as lyrRT
        oRastershp = os.path.join(os.path.dirname(maskShapeFile), Path(maskShapeFile).stem + "_rasterized.vrt")

        lyrRT.RasterizeVector(shpInput=maskShapeFile,
                              refRaster=dodPath,
                              output=oRastershp)
        maskArray = geoRT.RasterInfo(oRastershp).ImageAsArray()
        mask = maskArray < 1

    if extent:
        # from rasterio.windows import Window
        win = GetWindow(iRasterPath=dodPath, windowGeo=extent)

        # elevation = src.read(1, window=Window(win[0],win[1],win[2],win[3]))
        array = np.copy(array[win[1]:win[1] + win[-1], win[0]:win[0] + win[2]])
        if maskShapeFile:
            mask = np.copy(mask[win[1]:win[1] + win[-1], win[0]:win[0] + win[2]])

    if plotDistribution:
        import geospatialroutine.Filter.Filtering_Routine as filterRT
        if maskShapeFile:
            distArray = np.ma.masked_array(array, mask=~mask)

        if maxDisplacement:
            distArray = filterRT.MaskLargeValues(inputArray=distArray, maxDisplacement=maxDisplacement)
        
        pltRT.PlotDistribution(inputArray=distArray, xlim=[], title=title,
                               xLabel="Elevation Difference [m]", nbins=50)

    if basemap:
        basemapInfo = geoRT.RasterInfo(basemap)
        basemapArray = basemapInfo.ImageAsArray()
        if basemapInfo.EPSG_Code != rasterInfo.EPSG_Code or basemapInfo.pixelWidth != rasterInfo.pixelWidth:
            print(" Basemap and input raster do not have the same grid")
            print(" Resampling the base map to the same grid ")
            newBasemap = os.path.join(os.path.dirname(basemap), Path(basemap).stem + "_resampled.vrt")
            gdal.Warp(newBasemap, basemap, dstSRS='EPSG:' + str(rasterInfo.EPSG_Code), xRes=rasterInfo.pixelWidth,
                      yRes=rasterInfo.pixelHeight)
            basemapInfo = geoRT.RasterInfo(newBasemap)
            basemapArray = basemapInfo.ImageAsArray()
            basemap = newBasemap
        if extent:
            win = GetWindow(iRasterPath=basemap, windowGeo=extent)
            print(win)
            basemapArray = np.copy(basemapArray[win[1]:win[1] + win[-1], win[0]:win[0] + win[2]])
            print(basemapArray.shape)
        if basemapHillshade:
            hillshade = es.hillshade(basemapArray, azimuth=133, altitude=45)
            ax.imshow(hillshade, cmap="Greys", alpha=0.5)
        else:
            ax.imshow(basemapArray, cmap="Greys")
    if maskShapeFile:
        array_ma = np.ma.masked_array(array, mask=mask)
        array_ma = np.ma.masked_where(array_ma == noData, array_ma, copy=True)
    else:
        array_ma = np.ma.masked_where(array == noData, array, copy=True)





    im = ax.imshow(array_ma, interpolation=None, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title)

    xStep = int(array_ma.shape[1] / 4)
    yStep = int(array_ma.shape[0] / 4)
    xAxisVals = np.arange(xStep, array_ma.shape[1] - xStep, step=xStep)
    yAxisVals = np.arange(yStep, array_ma.shape[0] - yStep, step=yStep)

    if rasterInfo.EPSG_Code != 4326 and len(xAxisVals) == len(yAxisVals):
        xyMapCoords = rasterInfo.Pixel2Map_Batch(X=xAxisVals, Y=yAxisVals)
        latLongCoords = geoRT.ConvCoordMap1ToMap2_Batch(X=xyMapCoords[0], Y=xyMapCoords[1], targetEPSG=4326,
                                                        sourceEPSG=rasterInfo.EPSG_Code)
        ax.set_xticks(xAxisVals)
        lats = [round(val, 2) for val in latLongCoords[0]]
        longs = [round(val, 2) for val in latLongCoords[1]]
        ax.set_xticklabels(longs)
        ax.set_yticks(yAxisVals)
        ax.set_yticklabels(lats)

    ax.set_xlabel("Lon [$^\circ$]", fontsize=14)
    ax.set_ylabel("Lat [$^\circ$]", fontsize=14)
    # plt.show()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pltRT.ColorBar_(ax=ax, mapobj=im, cmap=cmap, vmin=vmin, vmax=vmax, size=14, label=" Elev. Diff. [m]")
    if saveFig:
        plt.savefig(os.path.join(os.path.dirname(dodPath), Path(dodPath).stem + ".svg"), dpi=400)
    else:
        plt.show()
    return


if __name__ == '__main__':
    # # # main()
    # # # flowLineCoord = ReadShapeFile("/home/cosi/2-Data/4-Shisper/0-Shp_ROI_KML/FlowLine/ShishperFlowLine/ShisperCenterFlowLine.shp")
    # # path = "E:\OneDrive - California Institute of Technology\\04-Projects\\01-Shisper_Project\\0-Results\\2-DEM Results"
    # # csv1 = os.path.join(path, "2-Elevation_profile_DEM2019_DEM2014_Mochware.csv")
    # # csv2 = os.path.join(path, "1-Elevation_profile_DEM2019_DEM2014_Shisper.csv")
    # # csv3 = os.path.join(path, "3-ElevationDIfference_2019-2014_ElevationVSdistance_shisper.csv")
    # # csv4 = os.path.join(path, "4-ElevationDIfference_2019-2014_ElevationVSdistance_Mochwar.csv")
    # # csvPath = "//media/storage/Saif/VB/Sharefiles/Profile_stableArea_GEWV_SRTM.csv"
    # tempPath = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-3D project-Planet\DEMS_results_figures4Paper\\Figures"
    # # data = PlotProfileFromCSV(csvPath=csv2,nbProfile=2,labels=["SRTM","DEM2019"])
    # # # PlotDem(data)
    # # # PlotDEMv2()
    # # # PlotDEMv3()
    # # # PlotProfilesAndHHDEM(csvPath1=csv1,csvPath2=csv2)
    # # # path = "H:\\2-Data\\4-Shisper\\4-Shisper_3D\\2-DEM_worlSpace_Planet\PS_WorkSpace\\1-2017\DEMs"
    # # path = "H:\\2-Data\\4-Shisper\\4-Shisper_3D\\2-DEM_worlSpace_Planet\PS_WorkSpace\\3-2019"
    # csvFiles = FileRT.GetFilesBasedOnExtension(path=tempPath, filter="*.csv")
    # # PlotProfileFromCSV(csvPath=csvFiles[6], nbProfile=2, xLabel="Distance (km)",xlim=[-25,25],nbins=75)
    # #
    #
    # # csvFiles = FileRT.GetFilesBasedOnExtension(path=tempPath, filter="*.csv")
    # # PlotProfileFromCSV(csvPath=csvFiles[1], nbProfile=1, xLabel="Distance (km)")
    #
    #
    # ## Plot DEM
    # path  = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-3D project-Planet\DEMS_results_figures4Paper\\Figures"
    # dtmList = FileRT.GetFilesBasedOnExtension(path,filter="6-*.tif")
    # # for dem in  dtmList:
    # path = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-3D project-Planet\DEMS_results_figures4Paper\Figures\DEM_figures\Elevation_mask_shp\\masked.tif"
    # maskFile = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-3D project-Planet\DEMS_results_figures4Paper\Figures\DEM_figures\Elevation_mask_shp\\Mask.tif"
    # basemap = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-3D project-Planet\DEMS_results_figures4Paper\Figures\DEM_figures\Elevation_mask_shp\\basemapDEM.tif"
    #
    # # Plot_DEM_withMask(dtmPath=path,maskPath=maskFile,vmin=-50,vmax=50,basemap=basemap)
    # # Plot_ListDEM(dtmList)
    #
    # skysatDTM = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-3D project-Planet\SkySat\Results_for_Paper\\DEM_3m_MUL_Left_Nadir_27moffsetCorrection_.tif"
    # skysatDiff = "G:\SkySatData\Morenci_Mine_AZ\Morenci_Mine_AZ\L1B\MUL_basic_analytic_converted\DF_2019-3DEP_27offsetCorrected.tif"
    # sybset = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-3D project-Planet\SkySat\Results_for_Paper\\imageSpace_subset.tif"
    # Plot_DEM(dtmPath=sybset)#,vmax=200,vmin=-200)
    # skysatProfiles = "G:\SkySatData\Morenci_Mine_AZ\Morenci_Mine_AZ\L1B\MUL_basic_analytic_converted"
    # csvFiles = FileRT.GetFilesBasedOnExtension(path=skysatProfiles, filter="*.csv")
    # # PlotProfileFromCSV(csvPath=csvFiles[3], nbProfile=3, xLabel="Distance (km)", xlim=[-20, 20], nbins=50)

    path = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-3D project-Planet\SkySat\Results_for_Paper\\image-object.tif"
    # demInfo = GetRasterInfo(inputRaster=path)
    # demArray = ImageAsArray(demInfo)
    # PlotDistribution(demArray)
    # path = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-Planet_projects\\1-3D project-Planet\Results_figures4Paper\SkySat_Results_for_Paper"
    # csv = os.path.join(path, "Profile_ImageSapce_objectSpace_2.csv")
    # PlotProfileFromCSV(csvPath=csv)
    # path = "E:\OneDrive - California Institute of Technology\\04-Projects\\6-Planet_projects\\1-3D project-Planet\Results_figures4Paper\SkySat_Results_for_Paper\SkySat_DEM_morinci_results\WV_DEMs"
    # diff1_noInterpolation = os.path.join(path,"Diff_WVFilled_Skysat.tif")
    #
    # # diff1_noInterpolation = os.path.join(path, "temp.tif")
    # diff1_noInterpolation = os.path.join(path, "Diff_WVnonFilled_Skysat_.tif")
    # dtm = os.path.join(path,"19JAN08_DEM_NOFilleGaps_fromdensePointCloud.tif")
    # dtm = os.path.join(path,"19JAN08_DEM_NOFilleGaps_fromdensePointCloud_subset.tif")
    #
    # # path ="E:\OneDrive - California Institute of Technology\\04-Projects\\6-Planet_projects\\1-3D project-Planet\Results_figures4Paper\SkySat_Results_for_Paper\SkySat_DEM_morinci_results\WV_DEMs\\18_oct"
    # # diff1_noInterpolation = os.path.join(path, "Diff.tif")
    # rasterInfo = GetRasterInfo(diff1_noInterpolation,printInfo=True)
    #
    # array= ImageAsArray(rasterInfo)
    # array = np.ma.masked_where(array == -32767, array)
    # array = np.ma.masked_where(array > 30, array)
    # array = np.ma.masked_where(array <-30, array)
    #
    # # import rasterio
    # # from rasterio.plot import show, show_hist
    # #
    # # src = rasterio.open(path)
    # # show(src)
    # # show_hist((src, 1), bins=500, lw=0.0, stacked=False, alpha=0.3, histtype='stepfilled', title = "Histogram",vmin=-50,vmax=50)
    # # PlotDistribution(inputArray=array,xlim=[-20,20],nbins=35,xLabel="Elevation difference (m)",title="")
    # Plot_DEM(dtmPath=diff1_noInterpolation,vmin=-50,vmax=50)
    ####################################################################################################################
    ##                                   Example Plot DoD plotting                                                   ###
    ####################################################################################################################
    demPath = "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/DoD_py/DoD_DEM_2017_3.8m_alig_r0-AsterDEM_DEM_2019_3.75m_alig_r0-AsterDEM.tif"
    basemap = "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/0-AsterDEM_UTM.tif"
    shpMask = "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/ShisperMask.geojson"

    demPathList = FileRT.GetFilesBasedOnExtension(
        "/home/cosicorr/0-WorkSpace/PlanetProject/AutoCal_Planet_workSpace/Planet_PS/Results_PlanetLabs_PS/Shisper_3D/SelectedData4Paper/DoD_py")
    tilteList = ["Aster_Ps2016", "Aster_Ps2017", "Aster_Ps2019", "Aster_Ps2020", "Ps2016_Ps2017", "Ps2016_Ps2019",
                 "Ps2016_Ps2020", "Ps2017_Ps2019", "Ps2017_Ps2020", "Ps2019_Ps2020"]
    for i, demPath in enumerate(demPathList):
        Plot_DoD(dodPath=demPath, vmin=-50, vmax=50, cmap="RdYlBu", title=tilteList[i], basemap=basemap,
                 basemapHillshade=True,
                 extent=[451413.1483, 472407.3137, 4019782.7064, 4036926.6407], maskShapeFile=shpMask, saveFig=True)
