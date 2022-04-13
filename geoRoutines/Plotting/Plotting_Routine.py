# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2020

import os, gdal
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sns
import earthpy.spatial as es
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm
from pathlib import Path

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

from geoRoutines.georoutines import GetWindow
import geoRoutines.FilesCommandRoutine as FileRT
import geoRoutines.georoutines as geoRT


# from geoRoutines.Remove_from_PublicRelease.Routine import Distribution
# import geoRoutines.Remove_from_PublicRelease.Routine as RT

def BackGroundImg(backGroundImg):
    # Create back ground Image
    subsetImgInfo = RT.GetRasterInfo(inputRaster=backGroundImg)
    if subsetImgInfo["NbBands"] == 1:
        subsetImgAsRaster = RT.ImageAsArray(imageInfo=subsetImgInfo)
        subsetImgAsRasterNorm = RT.NormalizeImage(img=subsetImgAsRaster)
        return subsetImgAsRasterNorm
    if subsetImgInfo["NbBands"] == 3:
        for i in range(subsetImgInfo["NbBands"]):
            red = RT.ImageAsArray(imageInfo=subsetImgInfo, bandNumber=1)
            cropBand1Norm = RT.NormalizeImage(red)

            cropBand2 = RT.ImageAsArray(imageInfo=subsetImgInfo, bandNumber=2)
            cropBand2Norm = RT.NormalizeImage(cropBand2)

            cropBand3 = RT.ImageAsArray(imageInfo=subsetImgInfo, bandNumber=3)
            cropBand3Norm = RT.NormalizeImage(cropBand3)

            imgFinal = np.dstack((cropBand1Norm, cropBand2Norm, cropBand3Norm))
            return imgFinal


def ColorBar(ax, mapobj, width="3%", height="50%", label="Disp.[m]", ticks=None, size=16):
    """
    'upper right': 1,
    'upper left': 2,
    'lower left': 3,
    'lower right': 4,
    'right': 5,
    'center left': 6,
    'center right': 7,
    'lower center': 8,
    'upper center': 9,
    'center': 10
    """
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax,
                       width=width,  # width = 5% of parent_bbox width
                       height=height,  # height : 50%
                       loc=6,
                       bbox_to_anchor=(1, 0, 1, 1),  # [x,y,height,width]
                       bbox_transform=ax.transAxes,
                       borderpad=0.3)
    # cbar = plt.colorbar(imageArray, pad=0.1, cax=axins)  # , ax=axes.ravel().tolist()
    if ticks:
        cbar = plt.colorbar(mapobj, cax=axins, ticks=ticks)
    else:
        cbar = plt.colorbar(mapobj, cax=axins)

    cbar.set_label(label=label, size=size)
    cbar.ax.tick_params(labelsize=size)

    return


def ColorBar_(ax, mapobj, cmap, vmin, vmax, label="Disp.[m]", width="3%", height="50%",
              orientation='vertical',  # horizontal',
              bounds=None,
              extend="neither", size=8):
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    # fig, ax = plt.subplots(figsize=(6, 1))
    # fig.subplots_adjust(bottom=0.5)

    # cmap = mpl.cm.cool
    axins = inset_axes(ax,
                       width=width,  # width = 5% of parent_bbox width
                       height=height,  # height : 50%
                       loc=6,
                       bbox_to_anchor=(1, 0, 1, 1),  # [x,y,height,width]
                       bbox_transform=ax.transAxes,
                       borderpad=0.3)
    if bounds == None:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)

    cbar = plt.colorbar(mapobj,
                        cax=axins, orientation=orientation, label=label, spacing='uniform')
    cbar.set_label(label=label, size=size)
    cbar.ax.tick_params(labelsize=size)

    return


def Colorbar_2(mapobj, size="5%", pad=0.09):
    """Adjust colorbar height to match the matplotlib axis height.

    NOTE: This function requires matplotlib v 3.0.1 or greater or v 2.9 or
    lower to run properly.

    Parameters
    ----------
    mapobj : matplotlib axis object
        The image that the colorbar will be representing as a matplotlib axis
        object.
    size : char (default = "3%")
        The percent width of the colorbar relative to the plot.
    pad : int (default = 0.09)
        The space between the plot and the color bar.

    Returns
    -------
    matplotlib.pyplot.colorbar

        Matplotlib color bar object with the correct width that matches the
        y-axis height.

    """

    try:
        ax = mapobj.axes
    except AttributeError:
        raise AttributeError(
            "The colorbar function requires a matplotlib axis object. "
            "You have provided a {}.".format(type(mapobj))
        )
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=size, pad=pad, )
    return fig.colorbar(mapobj, cax=cax)


########################################################################################################################
def PlotMatlabFig(imgPath):
    from scipy.io import loadmat
    x = loadmat(imgPath)
    print(x)

    return


def PlotImg(imgPath):
    import matplotlib.image as mpimg
    img = mpimg.imread(imgPath)

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    widths = [1]
    heights = [1, 2]
    spec = fig.add_gridspec(nrows=1, ncols=1)
    ax1 = fig.add_subplot(spec[0, 0])
    ax1.imshow(img)
    plt.show()

    return


def Plot_Raster(imgArray1, figureTitle="", initAnimation=False):
    """

    :param imgArray1:
    :param figureTitle:
    :param initAnimation: if true, this function will be used to initiate the animation
                          if False, this function will be used to visualize a single image
    :return:
    """
    fig = plt.figure()
    if figureTitle:
        title = fig.suptitle(figureTitle, fontsize=10)

    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])  # row,col

    # ax1.axis("off")
    # ax2.axis("off")
    im1 = ax1.imshow(imgArray1, "gray", rasterized=True)

    if not initAnimation:
        plt.show()
    if initAnimation:
        ax1.axis("off")
        return im1, fig


def DispDiplacmentMap(displacmentRasterPath, bandNumber=1, backgroundImg=None, transparency=True, vMin=0, vMax=3,
                      label="Displacement (m/day)"):
    """

    :param displacmentRasterPath:
    :param backgroundImg:
    :param transparency:
    :param vMin:
    :param vMax:
    :return:
    """
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])  # row,col
    if backgroundImg:
        # Create a grid
        backgroundImgArray = BackGroundImg(backGroundImg=backgroundImg)
        ax1.imshow(backgroundImgArray)

    my_cmap = plt.get_cmap('rainbow')
    # my_cmap = plt.get_cmap('gist_earth')
    # my_cmap= plt.get_cmap("RdBu")
    # my_cmap = plt.get_cmap('seismic')
    # my_cmap = plt.get_cmap('gray')

    imageAsArray = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=displacmentRasterPath), bandNumber=bandNumber)
    if transparency:
        print("Plot with applying mask")
        import Filtering_Routine
        imageAsArray = Filtering_Routine.MaskZeros(inputArray=imageAsArray)
    im = ax1.imshow(imageAsArray, cmap=my_cmap, vmin=vMin, vmax=vMax, origin="upper")
    # cbar = plt.colorbar(im, label=label)
    ColorBar(ax=ax1, imageArray=im)
    savingFolder = "//home/cosi/2-Data/4-Shishper/FigForPaper/PCA/"
    # plt.savefig(os.path.join(savingFolder, "velocity.pdf"), format="pdf")
    plt.show()


def DispDiplacmentMap_fromArray(displacmentRasterArray, bandNumber=1, backgroundImg=None, vMin=0, vMax=3,
                                label="Displacement (m/day)"):
    """

    :param displacmentRasterPath:
    :param backgroundImg:
    :param transparency:
    :param vMin:
    :param vMax:
    :return:
    """
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])  # row,col
    if backgroundImg:
        # Create a grid
        backgroundImgArray = BackGroundImg(backGroundImg=backgroundImg)
        ax1.imshow(backgroundImgArray)

    my_cmap = plt.get_cmap('rainbow')
    # my_cmap = plt.get_cmap('gist_earth')
    # my_cmap= plt.get_cmap("RdBu")
    # my_cmap = plt.get_cmap('seismic')
    # my_cmap = plt.get_cmap('gray')

    imageAsArray = np.array(displacmentRasterArray)
    im = ax1.imshow(imageAsArray, cmap=my_cmap, vmin=vMin, vmax=vMax, origin="upper")
    # cbar = plt.colorbar(im, label=label)
    ColorBar(ax=ax1, imageArray=im)

    plt.show()


def DispProfileOnRaster(rasterInfo, rasterArray, profileFilePath, vMin=-1.5, vMax=1.5, transparency=True):
    """

    :param path:
    :param vMin:
    :param vMax:
    :return:
    """

    profileArray = np.loadtxt(profileFilePath)
    print(profileArray)
    xPixList = []
    yPixList = []
    y = []
    listIndex = [199, 428, 655, 880]
    xPixLine = []
    yPixLine = []

    for index, xMap in enumerate(profileArray[:, 1]):
        xPix, yPix = RT.Map2Pixel(imageInfo=rasterInfo, x=xMap, y=profileArray[index, 2])
        xPixList.append(xPix)
        yPixList.append(yPix)
        y.append(yPix + profileArray[index, -1])
        if index in listIndex:
            xPixTmp, yPixTmp = RT.Map2Pixel(imageInfo=rasterInfo, x=xMap, y=profileArray[index, 2])
            xPixLine.append(xPixTmp)
            yPixLine.append(yPixTmp)

    yData = profileArray[:, -1]

    yCCD = []
    xCCD = []
    for index_, index in enumerate(listIndex):
        print(index_)
        if index_ == 0:
            xCCD.append(xPixList[0:index])
            yCCD.append(yData[0:index])
        else:
            xCCD.append(xPixList[listIndex[index_ - 1]:index])
            yCCD.append(yData[listIndex[index_ - 1]:index])

    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)

    colorList = ["b", "g", "c", "m", "k", "r", "y"]
    ax1.plot(xPixList, yData)
    ax1.set_xlim(0, rasterInfo["Dimension"][0])
    for number, xList in enumerate(xCCD):
        ax1.plot(xList, yCCD[number], colorList[number])

    for xPixLine_ in xPixLine:
        ax1.axvline(x=xPixLine_, c="k", ls="--")
    ax1.grid()

    my_cmap = plt.get_cmap('RdYlBu')
    im = ax2.imshow(rasterArray, cmap=my_cmap, origin="upper", vmin=vMin, vmax=vMax)
    ax2.plot(xPixList, y, 'r-')
    ax2.plot(xPixLine, yPixLine, "b|")
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical', label=' EW (m)')

    plt.show()

    return


def MapDistribution(inputMap, xlim=[-250, 250], ylim=[0, 3]):
    rasterInfo = RT.GetRasterInfo(inputMap)
    sample = RT.ImageAsArray(rasterInfo)
    sample = np.ma.masked_invalid(sample)
    # Remove mask and array to vector
    if isinstance(sample, np.ma.MaskedArray):  # check if the sample was masked using the class numpy.ma.MaskedArray
        sample = sample.compressed()  ## return all the non-masked values as 1-D array
    else:
        if sample.ndim > 1:  # if the dimension of the array more than 1 transform it to 1-D array
            sample = sample.flatten()

    # Estimate initial sigma and RMSE
    (mu, sigma) = norm.fit(sample)
    print("mu, sigma", mu, sigma)
    RMSE = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(sample))))
    print("RMSE=", RMSE)
    max = '%.3f' % (np.nanmax(sample))
    min = '%.3f' % (np.nanmin(sample))
    std = '%.3f' % (np.nanstd(sample))
    mean = '%.3f' % (np.nanmean(sample))
    print(max, min, std, mean)

    print("Plotting ....")
    lower_bin = mu - 3 * sigma
    upper_bin = mu + 3 * sigma
    gs = gridspec.GridSpec(1, 1)
    fig = plt.figure()  # figsize=(10, 10))
    ax1 = plt.subplot(gs[0, 0])  # row 0, col 0
    hist, bins = np.histogram(sample, range=[lower_bin, upper_bin], density=False, bins=1000)
    bars = ax1.bar((bins[:-1] + bins[1:]) / 2, hist, align='center', width=(bins[1] - bins[0]))
    title_ = "Displacement distribution before Offset Correction\n(full image)"
    ax1.set_title(title_)
    # axes[0].set_ylim(ylim[0], ylim[1])
    ax1.set_xlim(xlim[0], xlim[1])
    txt = "mean=" + str(mean) + "\nmax=" + str(max) + "\nmin=" + str(min) + "\nRMSE=" + str(RMSE) + "\nstd=" + str(std)
    txt1 = ax1.text(0.8, 0.9, txt,
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1.transAxes)
    plt.xlabel('Displacement [m/day]')
    plt.ylabel(' Number of Samples ')

    # savingFolder = "/home/cosi/2-Data/4-Shishper/FigForPaper/Offest_Correction/"
    # plt.savefig(os.path.join(savingFolder, "distribution_before.eps"), format="eps")
    plt.show()


def Plot_2_maps(imgArray1, imgArray2, backgroundImg=None, histogram=True, animationRecording=False, figTitle="fig",
                title1="",
                title2="", figureTitle="", animation=False, vMin=-1, vMax=1, transparency=False):
    imgArray1 = np.ma.masked_invalid(imgArray1)
    imgArray2 = np.ma.masked_invalid(imgArray2)
    fig = plt.figure()
    if figureTitle:
        title = fig.suptitle(figureTitle, fontsize=10)
    if not histogram:
        gs = gridspec.GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])  # row,col
        ax2 = fig.add_subplot(gs[0, 1])  # row,col

    else:
        gs = gridspec.GridSpec(2, 2)
        ax1 = fig.add_subplot(gs[0, 0])  # row,col
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])

    ax1.set_title(title1, fontsize=10)
    ax2.set_title(title2, fontsize=10)
    if backgroundImg:
        subsetImgAsRasterNorm = BackGroundImg(backGroundImg=backgroundImg)
        ax1.imshow(subsetImgAsRasterNorm, cmap="gray")
        ax2.imshow(subsetImgAsRasterNorm, cmap="gray")
    # ax1.axis("off")
    # ax2.axis("off")

    # my_cmap = plt.get_cmap('gist_ncar')
    # my_cmap = plt.get_cmap('jet')
    # my_cmap = plt.get_cmap('rainbow')
    # my_cmap = plt.get_cmap('Reds')
    my_cmap = plt.get_cmap('seismic')
    my_cmap = plt.get_cmap('RdBu')

    if transparency == True:
        imageAsArray1 = np.ma.masked_inside(imgArray1, 0.00001, -0.00001)
        imageAsArray2 = np.ma.masked_inside(imgArray2, 0.00001, -0.00001)

    vmin = vMin
    vmax = vMax

    im1 = ax1.imshow(imgArray1, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)
    # im1 = ax1.pcolormesh(imgArray1, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)
    im2 = ax2.imshow(imgArray2, cmap=my_cmap, vmin=vmin, vmax=vmax, rasterized=True)
    if histogram:
        ax3.hist(imgArray1.ravel(), bins=1000, range=(np.nanmin(imgArray1), np.nanmax(imgArray1)), fc='k', ec='k')
        his1, bin_edges1 = np.histogram(imgArray1, bins=1000, range=(np.nanmin(imgArray1), np.nanmax(imgArray1)))
        ax3.plot(bin_edges1[0:-1], his1)  # <- or here
        ax4.hist(imgArray1.ravel(), bins=1000, range=(np.min(imgArray2), np.max(imgArray2)), fc='k', ec='k')
        his2, bin_edges2 = np.histogram(imgArray2, bins=1000, range=(np.nanmin(imgArray2), np.nanmax(imgArray2)))
        ax4.plot(bin_edges2[0:-1], his2)  # <- or here
    ##Legend
    ColorBar(ax=ax2, mapobj=im2)
    if not animation:
        plt.show()
    if animation:
        ax1.axis("off")
        ax2.axis("off")
        return im1, im2, fig
    return


def PlotPlaneFitting(arrayBeforeCorrection, arrayAfterCorrection, planeArray, save=False, show=True, vMin=None,
                     vMax=None, diagnoseSavingPath=None, imgIndex=None, histogram=True):
    fig = plt.figure(figsize=(10, 3))

    if not histogram:
        gs = gridspec.GridSpec(1, 3)
        ax1 = fig.add_subplot(gs[0, 0])  # row,col
        ax2 = fig.add_subplot(gs[0, 1])  # row,col
        ax3 = fig.add_subplot(gs[0, 2])

    else:
        gs = gridspec.GridSpec(2, 3)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
        ax4 = fig.add_subplot(gs[1, 0])
        ax5 = fig.add_subplot(gs[1, 1])
        ax6 = fig.add_subplot(gs[1, 2])
    cmap = 'seismic'
    cmap = plt.get_cmap('RdBu')
    # cmap = "rainbow"
    # cmap = 'gist_earth'

    print("Plotiing ....")

    if not vMax:
        vis_max = np.nanmax(arrayAfterCorrection)
    else:
        vis_max = vMax

    if not vMin:
        vis_min = np.nanmin(arrayAfterCorrection)
    else:
        vis_min = vMin
    # vis_min =-0.5
    # vis_max = 0.5
    ax1.imshow(arrayBeforeCorrection, cmap=cmap, vmin=vis_min, vmax=vis_max)
    ax1.set_title('Input Image', fontsize=12)
    ax1.tick_params(labelsize=12)

    # plot
    ax2.imshow(planeArray, cmap=cmap, vmin=vis_min, vmax=vis_max)
    ax2.set_title('Fitted plane', fontsize=12)
    ax2.tick_params(labelsize=12)

    im = ax3.imshow(arrayAfterCorrection, cmap=cmap, vmin=vis_min, vmax=vis_max)
    ax3.set_title('After plane correction', fontsize=12)
    ax3.tick_params(labelsize=12)
    fig.subplots_adjust(right=0.91)

    # divider = make_axes_locatable(axes[2])
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # cb = fig.colorbar(im, cax=cax)
    # cb.ax.tick_params(labelsize=12)
    ColorBar(ax=ax3, mapobj=im)
    ax1.axis("off")
    ax2.axis("off")
    ax3.axis("off")
    if histogram:
        ax4.hist(arrayBeforeCorrection.ravel(), bins=1000,
                 range=(np.nanmin(arrayBeforeCorrection), np.nanmax(arrayBeforeCorrection)), fc='k', ec='k')
        his1, bin_edges1 = np.histogram(arrayBeforeCorrection, bins=1000,
                                        range=(np.nanmin(arrayBeforeCorrection), np.nanmax(arrayBeforeCorrection)))
        ax4.plot(bin_edges1[0:-1], his1)

        ax5.hist(planeArray.ravel(), bins=1000, range=(np.min(planeArray), np.max(planeArray)), fc='k', ec='k')
        his2, bin_edges2 = np.histogram(planeArray, bins=1000, range=(np.nanmin(planeArray), np.nanmax(planeArray)))
        ax5.plot(bin_edges2[0:-1], his2)

        ax6.hist(arrayAfterCorrection.ravel(), bins=1000,
                 range=(np.min(arrayAfterCorrection), np.max(arrayAfterCorrection)), fc='k', ec='k')
        his3, bin_edges3 = np.histogram(arrayAfterCorrection, bins=1000,
                                        range=(np.nanmin(arrayAfterCorrection), np.nanmax(arrayAfterCorrection)))
        ax6.plot(bin_edges3[0:-1], his3)  # <- or here

        # ax4.axis("off")
        # ax5.axis("off")
        # ax3.axis("off")
    if show:
        plt.show()
    if save:
        plt.draw()
        plt.pause(1)
        fig.savefig(diagnoseSavingPath + "img_" + str(imgIndex) + ".png")
    plt.close()

    return


def Animation_2_Maps(im1, im2, array1List, array2List, fig, figTitleList, saveAnimation=False, animationSavingPath="",
                     animationTitle="", animationType="mp4"):
    def updatefig(i):
        print("Frame =", i + 1)

        im1.set_array(array1List[i])
        im2.set_array(array2List[i])
        # if not oneDate:
        #     txtTitle = dateList[i] + " (" + satType[i] + ")" + " diff=" + str(difDays[i]) + "  ID:" + str(i + 1)
        # else:
        #
        #     if satType:
        #         txtTitle = dateList[i] + " (" + satType[i] + ")" + "  ID:" + str(i + 1)
        #     else:
        #         txtTitle = dateList[i]  # +  "  ID:" + str(i + 1)
        #
        # txt = ax3.text(170, 15, txtTitle, color="red", fontsize=8, bbox=dict(facecolor='white', alpha=1))
        fig.suptitle(figTitleList[i], fontsize=10)
        return im1, im2
        ## define the postion and parameters of the color bar

    # plt.colorbar(stackImg, label=' Velocity (m/day)')
    # axins = inset_axes(ax3,
    #                    width="3%",  # width = 5% of parent_bbox width
    #                    height="30%",  # height : 50%
    #                    loc='lower left',
    #                    bbox_to_anchor=(1.1, 0.4, 1, 1),
    #                    bbox_transform=ax3.transAxes,
    #                    borderpad=0.1)
    # cbar = plt.colorbar(stackImg, pad=0.1, cax=axins, ticks=[vMin, 1, 2, vMax])  # , ax=axes.ravel().tolist()
    # cbar.set_label(label=' Velocity (m/day)', size=10, labelpad=-50)
    # cbar.set_ticklabels(['0', '1', '2', '>3'])
    # cbar.ax.tick_params(labelsize=10)

    ani = animation.FuncAnimation(fig, updatefig, frames=len(array1List), interval=1000,
                                  repeat_delay=3000)

    if saveAnimation == True:
        if animationType == "mp4":
            ani.save(os.path.join(animationSavingPath, animationTitle + "." + animationType), fps=0.5,
                     dpi=800)  # , writer="imagemagick")#,for GIF
        if animationType == "gif":
            ani.save(os.path.join(animationSavingPath, animationTitle + "." + animationType), dpi=100,
                     writer="imagemagick")  # ,for GIF
    plt.show()

    return


def Animation_Rasters_Maps(im1, array1List, fig, figTitleList, saveAnimation=False, animationSavingPath="",
                           animationTitle=""):
    def updatefig(i):
        if i >= len(array1List):
            i = len(array1List) - 1
        print("Frame =", i)
        im1.set_array(array1List[i])
        fig.suptitle(figTitleList[i], fontsize=10)
        return im1,

    ani = animation.FuncAnimation(fig, updatefig, frames=len(array1List) + 3, interval=1500)
    # repeat_delay=3000)

    if saveAnimation == True:
        # writer = animation.writers['avconv'](fps=0.5)
        ani.save(os.path.join(animationSavingPath, animationTitle), dpi=100, writer="imagemagick")  # ,for GIF

    plt.show()

    return


############################################### Distribution ###########################################################


def PlotDistribution(inputArray, xlim=[], ylim=[0, 3], xLabel='Displacement [m]',
                     title="Displacement distribution before Offset Correction\n(full image)", nbins=50, fontSize=16,
                     svgFig=""):
    """
    https://towardsdatascience.com/histograms-and-density-plots-in-python-f6bda88f5ac0
    """

    dist = geoRT.cgeoStat(inputArray=inputArray)
    print("Plotting ....")
    lower_bin = dist.mu - 4 * dist.sigma
    upper_bin = dist.mu + 4 * dist.sigma
    gs = gridspec.GridSpec(1, 1)
    fig2 = plt.figure()  # figsize=(10, 10))
    ax1 = plt.subplot(gs[0, 0])  # row 0, col 0
    hist, bins = np.histogram(dist.sample, range=[lower_bin, upper_bin], density=False, bins=nbins)
    bars = ax1.bar((bins[:-1] + bins[1:]) / 2, hist, align='center', width=(bins[1] - bins[0]), color="blue",
                   edgecolor='black')
    # sns.distplot(dist.sample, hist=False, kde=True, bins=nbins, color='darkblue',
    #              hist_kws={'edgecolor': 'black'}, kde_kws={'shade': True, 'linewidth': 4})
    title_ = title
    ax1.set_title(title_)
    # axes[0].set_ylim(ylim[0], ylim[1])
    if not any(xlim):
        ax1.set_xlim(float(dist.mean) - 4 * float(dist.std), float(dist.mean) + 4 * float(dist.std))
    else:
        ax1.set_xlim(xlim[0], xlim[1])
    txt = "Mean=" + str(dist.mean) \
          + "\nStd=" + str(dist.std) \
          + "\nMedian=" + str(dist.median) \
          + "\nNMAD=" + str(dist.nmad) \
          + "\nRMSE=" + str(dist.RMSE)
    txt1 = ax1.text(0.75, 0.9, txt,
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1.transAxes, fontsize=fontSize)
    ax1.set_xlabel(xLabel, fontsize=fontSize)
    ax1.set_ylabel(' Number of Samples ', fontsize=fontSize)
    # ax1.xtickes(fontsize=25)
    if svgFig != "":
        # savingFolder = "/home/cosi/2-Data/4-Shishper/FigForPaper/Offest_Correction/"
        plt.savefig(svgFig, dpi=400)

    plt.show()

    return


def PlotDistribution_Batch(inputArrayList, labels, colors, xlim=[-5, 5], ylim=[0, 3], xLabel='Displacement [m]',
                           title="Displacement distribution before Offset Correction\n(full image)", nbins=None,
                           fontSize=16, saveFig=""):
    distList = []
    for array_ in inputArrayList:
        dist = Distribution(inputArray=array_)
        distList.append(dist)
    fig = plt.figure()  # figsize=(10, 10))
    gs = gridspec.GridSpec(1, 1)
    ax1 = plt.subplot(gs[0, 0])  # row 0, col 0
    print("Plotting ....")

    offset = 3.5
    lowerBinList = []
    upperBinList = []

    for index, dist_ in enumerate(distList):
        lower_bin = dist_.mu - offset * dist_.sigma
        upper_bin = dist_.mu + offset * dist_.sigma
        lowerBinList.append(lower_bin)
        upperBinList.append(upper_bin)
        nbins_ = nbins
        if nbins == None:
            ## We use the Struge's rule to compute the optima number of bins
            nbins_ = int(1 + 3.322 * np.log(len(dist_.sample)))
            print("nbins=", nbins_)

        # sns.distplot(dist_.sample, hist=False, kde=True, bins=nbins_,
        #              hist_kws={"range": [lower_bin, upper_bin]},
        #              kde_kws={'shade': False, 'linewidth': 3, "color": colors[index]}, label=labels[index])
        sns.kdeplot(dist_.sample, label=labels[index], fill=True, alpha=0.4, linewidth=3, color=colors[index])
        # , hist=False, kde=True, bins=nbins_,
        # hist_kws={"range": [lower_bin, upper_bin]},
        # kde_kws={'shade': False, 'linewidth': 3, "color": colors[index]}, label=labels[index])
    title_ = title
    ax1.set_title(title_)
    # ax1.set_xlim(min(lowerBinList), max(upperBinList))
    ax1.set_xlim(xlim[0], xlim[-1])
    ax1.set_xlabel(xLabel, fontsize=fontSize)
    ax1.set_ylabel('Density ', fontsize=fontSize)
    # ax1.set_xticks(fontsize=18)
    ax1.tick_params(axis='x', labelsize=fontSize)
    ax1.tick_params(axis='y', labelsize=fontSize)
    # ax1.set_ytickes(fontsize=18)
    ax1.grid(True)
    plt.legend()
    if saveFig != "":
        plt.savefig(saveFig, dpi=400)
    else:
        plt.show()

    return


def Plot_scatter():
    import matplotlib.lines as pltL2D
    import pandas
    lim = 1.5
    # data = np.load("/media/cosicorr/storage/Saif2/3D_approch_Ridgecrest/0-3D_Correlation_Results/OrthoEval.txt",allow_pickle=True)
    data = pandas.read_csv("/media/cosicorr/storage/Saif2/3D_approch_Ridgecrest/0-3D_Correlation_Results/OrthoEval.txt",
                           sep="\\t")
    print(data["xError"])
    font = 12
    # Fixing random state for reproducibility

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.04

    rect_scatter = [left, bottom, width, height]

    # start with a rectangular Figure
    fig = plt.figure(figsize=(10, 10))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)

    # ax_histx.set_xlabel("XMap Error ")
    # ax_scatter.set_xlabel("Easting error (m) ", fontsize=font)
    # ax_scatter.set_ylabel("Northing Error (m)", fontsize=font)

    # the scatter plot:
    lines = []
    labels = []
    ax_scatter.scatter(data["xError"], data["yError"], s=35, color="k", linewidths=0.01)
    lines.append(pltL2D.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="k",
                               markersize=5))
    labels.append("nb. Ortho (" + str(len(data["xError"])) + ")")

    # now determine nice limits by hand:
    binwidth = 0.25

    ax_scatter.set_xlim((-lim, lim))
    ax_scatter.set_ylim((-lim, lim))

    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    ax_scatter.spines['left'].set_position('center')
    ax_scatter.spines['bottom'].set_position("center")
    # ax1.spines['bottom'].set_position(('axes', 20))

    # Eliminate upper and right axes
    ax_scatter.spines['right'].set_color('none')
    ax_scatter.spines['top'].set_color('none')
    ax_scatter.grid(True)

    # add text box for the statistics
    stats = ("meanEastingError=" + str(np.mean(np.asarray(list(data["xError"])))) + "\n" + "meanEastingError=" + str(
        np.mean(np.asarray(list(data["yError"])))) + "\n")
    ax_scatter.text(0.99, 1, stats, fontsize=font, transform=ax_scatter.transAxes, horizontalalignment='right')

    ax_scatter.legend(lines, labels, fontsize=font)
    plt.show()


def VisualizeDispDistribution(dispMapInfo, debug=False, vmin=None, vmax=None, factor=1):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import scipy.stats
    import seaborn as sns
    ewArray = dispMapInfo.ImageAsArray(bandNumber=1)
    nsArray = dispMapInfo.ImageAsArray(bandNumber=2)
    ewStat = geoRT.cgeoStat(inputArray=ewArray, displayValue=debug)
    nsStat = geoRT.cgeoStat(inputArray=nsArray, displayValue=debug)

    if vmin == None and vmax == None:
        nsvMin = float(nsStat.mean) - factor * float(nsStat.std)
        ewvMin = float(ewStat.mean) - factor * float(ewStat.std)
        nsvMax = float(nsStat.mean) + factor * float(nsStat.std)
        ewvMax = float(ewStat.mean) + factor * float(ewStat.std)

        vmin = min(ewvMin, nsvMin)
        vmax = max(ewvMax, nsvMax)

    fontSize = 18
    dpi = 600
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(1, 1)
    ax1 = plt.subplot(gs[0, 0])  # row 0, col 0
    xLabel = "Displacements [m]"
    sns.set()
    # # print("Plotting ....")

    ax1.set_title('Displacement density distribution')
    color = "#000080"
    sns.kdeplot(ewStat.sample, fill=True, alpha=0.4, linewidth=5, bw_adjust=1, label="EW disp.", color=color)
    sns.histplot(data=ewStat.sample, stat="density", bins=200, color=color, kde=False, ax=ax1, alpha=0.5, linewidth=1,
                 edgecolor=color)
    color = "#6B8E23"
    sns.histplot(data=nsStat.sample, stat="density", bins=200, kde=False, ax=ax1, alpha=0.5, linewidth=1, color=color,
                 edgecolor=color)
    sns.kdeplot(nsStat.sample, fill=True, alpha=0.4, linewidth=5, bw_adjust=1, label="NS disp.", color=color)

    ax1.set_xlabel(xLabel, fontsize=fontSize + 2)
    ax1.set_ylabel('Density ', fontsize=fontSize + 2)
    #
    ax1.tick_params(axis='x', labelsize=fontSize + 1)
    ax1.tick_params(axis='y', labelsize=fontSize + 1)
    ax1.set_xlim(vmin, vmax)
    ax1.grid(True)
    plt.legend(fontsize=fontSize)
    # plt.savefig(os.path.join(iFolder, "Displacement_density_distribution.svg"), dpi=600)
    plt.show()


def VisualizeCorrelation(iCorrPath,
                         ewArray,
                         nsArray,
                         snrArray=[],
                         vmin=None,
                         vmax=None,
                         title=False,
                         cmap="RdYlBu",
                         save=True,
                         show=False,
                         factor=1, dpi=100):
    # import geospatialroutine.georoutines as geoRT
    # corrInfo = geoRT.RasterInfo(iCorrPath)
    ewStat = geoRT.cgeoStat(inputArray=ewArray, displayValue=False)
    nsStat = geoRT.cgeoStat(inputArray=nsArray, displayValue=False)
    # nsMean = float(nsStat.mean)
    nsvMin, ewvMin = vmin, vmin
    nsvMax, ewvMax = vmax, vmax
    if vmin == None:
        nsvMin = float(nsStat.mean) - factor * float(nsStat.std)
        ewvMin = float(ewStat.mean) - factor * float(ewStat.std)

    if vmax == None:
        nsvMax = float(nsStat.mean) + factor * float(nsStat.std)
        ewvMax = float(ewStat.mean) + factor * float(ewStat.std)

    if len(snrArray) != 0:
        fig, axs = plt.subplots(1, 3, figsize=(16, 9))
        axs[0].imshow(snrArray, cmap="gray", vmin=0.3, vmax=0.95)
        axs[1].imshow(ewArray, cmap=cmap, vmin=vmin, vmax=vmax)
        im1 = axs[2].imshow(nsArray, cmap=cmap, vmin=vmin, vmax=vmax)
        for ax, title_ in zip(axs, ["SNR", "East/West", "North/South"]):
            ax.axis('off')
            ax.set_title(title_)
        # ColorBar_(ax=axs[-1], mapobj=im1, cmap=cmap, vmin=vmin, vmax=vmax, orientation="vertical")
    else:
        fig, axs = plt.subplots(1, 2, figsize=(16, 9))

        imEW = axs[0].imshow(ewArray, cmap=cmap, vmin=ewvMin, vmax=ewvMax)
        imNS = axs[1].imshow(nsArray, cmap=cmap, vmin=nsvMin, vmax=nsvMax)
        for ax, title_ in zip(axs, ["East/West", "North/South"]):
            ax.axis('off')
            ax.set_title(title_)
        ColorBar_(ax=axs[0], mapobj=imEW, cmap=cmap, vmin=ewvMin, vmax=ewvMax, orientation="vertical")
        ColorBar_(ax=axs[1], mapobj=imNS, cmap=cmap, vmin=nsvMin, vmax=nsvMax, orientation="vertical")

    if not title:
        fig.suptitle(Path(iCorrPath).stem)
    else:
        fig.suptitle(title)
    if save:
        plt.savefig(os.path.join(os.path.dirname(iCorrPath), Path(iCorrPath).stem + ".png"), dpi=dpi)
    if show:
        plt.show()
    fig.clear()
    plt.close(fig)
    return


def VisualizeDisplacment(dispPath, vmin=-50, vmax=50, cmap="RdYlBu", title=None, basemap=None, basemapHillshade=False,
                         extent=[], maskShapeFile=None, noData=-32767, saveFig=False, bandNumber=1, coordStep=4,
                         fontSize=12):
    """
    Extent: same extent as qgis
    """
    rasterInfo = geoRT.RasterInfo(dispPath)
    fig = plt.figure()  # figsize=(20, 20))
    ax = fig.gca()
    if not title:
        title = Path(dispPath).stem
    array = rasterInfo.ImageAsArray(bandNumber=bandNumber)

    if maskShapeFile:
        import geospatialroutine.Routine_Lyr as lyrRT
        oRastershp = os.path.join(os.path.dirname(maskShapeFile), Path(maskShapeFile).stem + "_rasterized.vrt")

        lyrRT.RasterizeVector(shpInput=maskShapeFile,
                              refRaster=dispPath,
                              output=oRastershp)
        maskArray = geoRT.RasterInfo(oRastershp).ImageAsArray()
        mask = maskArray < 1

    if extent:
        # from rasterio.windows import Window
        win = GetWindow(iRasterPath=dispPath, windowGeo=extent)

        # elevation = src.read(1, window=Window(win[0],win[1],win[2],win[3]))
        array = np.copy(array[win[1]:win[1] + win[-1], win[0]:win[0] + win[2]])
        if maskShapeFile:
            mask = np.copy(mask[win[1]:win[1] + win[-1], win[0]:win[0] + win[2]])
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
    plt.xticks(fontsize=fontSize)
    plt.yticks(fontsize=fontSize)
    plt.tick_params(axis='y', direction="in")
    plt.tick_params(axis='x', direction="in")
    plt.xticks(rotation=45)
    ColorBar_(ax=ax, mapobj=im, cmap=cmap, vmin=vmin, vmax=vmax)

    xStep = int(array_ma.shape[1] / coordStep)
    yStep = int(array_ma.shape[0] / coordStep)
    xAxisVals = np.arange(xStep, array_ma.shape[1] - xStep, step=xStep)
    yAxisVals = np.arange(yStep, array_ma.shape[0] - yStep, step=yStep)
    print(len(xAxisVals), len(yAxisVals))
    if rasterInfo.EPSG_Code != 4326:
        if len(xAxisVals) > len(yAxisVals):
            xAxisVals = xAxisVals[0:len(yAxisVals)]
        if len(yAxisVals) > len(xAxisVals):
            yAxisVals = yAxisVals[0:len(xAxisVals)]
        xyMapCoords = rasterInfo.Pixel2Map_Batch(X=xAxisVals, Y=yAxisVals)
        latLongCoords = geoRT.ConvCoordMap1ToMap2_Batch(X=xyMapCoords[0], Y=xyMapCoords[1], targetEPSG=4326,
                                                        sourceEPSG=rasterInfo.EPSG_Code)
        ax.set_xticks(xAxisVals)
        lats = [round(val, 2) for val in latLongCoords[0]]
        longs = [round(val, 2) for val in latLongCoords[1]]
        ax.set_xticklabels(longs)
        ax.set_yticks(yAxisVals)
        ax.set_yticklabels(lats)
    plt.tight_layout()
    # plt.show()
    if saveFig:
        plt.savefig(os.path.join(os.path.dirname(dispPath), Path(dispPath).stem + ".png"), dpi=600)
    else:
        plt.show()
    return


# ===========================================PLOT DATA ================================================================#
def Plot_residuals(sample, yLabel='Error [pix]', title=None, save=True, oFolder=None):
    stat = geoRT.cgeoStat(inputArray=sample, displayValue=False)

    print(title + ":mean:{:.3f} ,median:{:.3f}, std:{:.3f}, RMSE:{:.3f}".format(float(stat.mean),
                                                                                float(stat.median),
                                                                                float(stat.std),
                                                                                float(stat.RMSE)))
    fig, ax = plt.subplots()
    ax.scatter(np.arange(0, sample.shape[0], 1), sample, c="k")
    ax.plot(np.arange(0, sample.shape[0], 1), sample)
    ax.axhline(y=float(stat.mean), color='r', linewidth=4, linestyle='--')
    ax.grid()

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.tick_params(which='both', width=2, direction="in")
    ax.set_xlabel('#GCPs')
    ax.set_ylabel(yLabel)
    if title is not None:
        ax.set_title(title)
    if save == True:
        plt.savefig(os.path.join(oFolder, title + ".svg"), dpi=300)
    return


def Plot_residuals_before_after(sampleIni, sampleFinal, yLabel="Error [pix]", save=True, oFolder=None, title=None,
                                show=False):
    stat = geoRT.cgeoStat(inputArray=sampleIni, displayValue=False)
    print(
        title + "Ini :mean:{:.3f}, std:{:.3f}, RMSE:{:.3f}".format(float(stat.mean), float(stat.std), float(stat.RMSE)))
    fig, ax = plt.subplots()
    ax.scatter(np.arange(0, sampleIni.shape[0], 1), sampleIni, c="k")
    ax.plot(np.arange(0, sampleIni.shape[0], 1), sampleIni, label="Initial")
    ax.axhline(y=float(stat.mean), color='r', linewidth=4, linestyle='--')

    stat = geoRT.cgeoStat(inputArray=sampleFinal, displayValue=False)
    print(
        title + "Final :mean:{:.3f}, std:{:.3f}, RMSE:{:.3f}".format(float(stat.mean), float(stat.std),
                                                                     float(stat.RMSE)))

    ax.scatter(np.arange(0, sampleFinal.shape[0], 1), sampleFinal, c="g")
    ax.plot(np.arange(0, sampleFinal.shape[0], 1), sampleFinal, label="Optimized")
    ax.axhline(y=float(stat.mean), color='g', linewidth=4, linestyle='--')

    ax.grid()

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.tick_params(which='both', width=2, direction="in")
    ax.set_xlabel('#GCPs')
    ax.set_ylabel(yLabel)
    if title is not None:
        ax.set_title(title)
    ax.legend()
    if save == True:
        plt.savefig(os.path.join(oFolder, title + ".svg"), dpi=300)
    if show == True:
        plt.show()
    plt.close(fig)
    return


########################################### Examples ###################################################################
def Example():
    # # corrPath = "/home/cosi/2-Data/1-Ridgecrest/0-Sentinel2_WorkSpace/2_Correlate_SL_images_With_no_Ground_displacement/S2A_2019-08-02_VS_2019_08_12/Cosi-Corr/PostProcessing/Corr_NLMF_1p5-41-5"
    # # rasterInfo = RT.GetRasterInfo(inputRaster=corrPath)
    # # ewArray = RT.ImageAsArray(imageInfo=rasterInfo, bandNumber=1)
    # #
    # # profileFilePath = "/home/cosi/2-Data/1-Ridgecrest/0-Sentinel2_WorkSpace/2_Correlate_SL_images_With_no_Ground_displacement/S2A_2019-08-02_VS_2019_08_12/Cosi-Corr/Profiles/CCD_Profile/Profile_CCD.txt"
    # # DispProfileOnRaster(rasterInfo=rasterInfo, rasterArray=ewArray, profileFilePath=profileFilePath)
    # tmp_path1 = "/home/cosi/2-Data/4-Shishper/DoubleCheck/"
    #
    # img1 = "Master_2018-07-10_S2A_B04_vs_Slave_2018-07-15_S2B_B04_EW_Velocity_Masked.tif"
    # img2 = "Master_2018-07-10_S2A_B04_vs_Slave_2018-07-15_S2B_B04_NS_Velocity_Masked.tif"
    # img = tmp_path1 + img2
    # tci = "/home/cosi/2-Data/4-Shishper/DownSamplingImg/TCI_subset_2_L8.tif"
    # # DispDiplacmentMap(displacmentRasterPath=img, backgroundImg=tci, transparency=True, vMin=-3, vMax=3)
    # # MapDistribution(inputMap=img)
    # pngPath = "/home/cosi/2-Data/4-Shishper/4-Results/2-Png/1-VelocityMaps_png/Png/"
    # # PlotImg(imgPath=pngPath+"2013-05-18 (LC) Time Span:16 ID:1.png")
    #
    # inputFolder = "/home/cosi/2-Data/4-Shishper/4-Results/2-Png/Shishper_Grid/Shishper_original_velocities/Data/"
    # tci = "/home/cosi/2-Data/4-Shishper/DownSamplingImg/TCI_subset_2_L8.tif"
    # savingPath = "//home/cosi/2-Data/4-Shishper/4-Results/2-Png/1-VelocityMaps_png/"
    #
    # # SavePlots(inputFolder=inputFolder, backgroundImg=tci, savingPath=savingPath, saveFig=True,
    # #                   transparency=True)
    # # SaveImgsOnGrids(inputFolder=inputFolder, savingFolder=savingPath, gridCouple=(3, 4), backgroundImg=tci,
    # #                 transparency=True, titleType=2)
    #
    # # Plot_scatter()
    #
    # ################################################# Start ############################################################
    # # path = "/home/cosicorr/0-WorkSpace/Temp/Diff_between_3DEP_and_PreDEM.tif"
    # # path = "//home/cosicorr/0-WorkSpace/3D-Correlation_project/Compare_DSM_diff_and_couple33_PreDEM"
    # #
    # # imgPath1 = os.path.join(path, "DSM_diff_subset.tif")
    # # tempArray1 = RT.ImageAsArray(RT.GetRasterInfo(imgPath1), bandNumber=1)
    # # array1 = np.ma.masked_where(tempArray1 == -32767.0, tempArray1)
    # # imgPath2 = os.path.join(path, "Couple33_PreDEM_3DDispalcment_subset.tif")
    # # tempArray2 = RT.ImageAsArray(RT.GetRasterInfo(imgPath2), bandNumber=3)
    # # array2 = np.ma.masked_where(tempArray2 == -32767.0, tempArray2)
    # # imgPath3 = os.path.join(path, "Diff_Couple33_Dz_dDSM.tif")
    # # tempArray3 = RT.ImageAsArray(RT.GetRasterInfo(imgPath3), bandNumber=1)
    # # array3 = np.ma.masked_where(tempArray3 == -32767.0, tempArray3)
    # #
    # #
    # # PlotDistribution_Batch(inputArrayList=[array1, array2, array3],
    # #                        labels=[Path(imgPath1).stem, Path(imgPath2).stem, Path(imgPath3).stem], nbins=None)
    # # # PlotDistribution(inputArray=array2,title="diff_Dz_Couple33 and dDSMs",xlim=[-2,2])
    # ################################################# END ############################################################
    #
    # ################################################# Start ############################################################
    # # path = "//home/cosicorr/0-WorkSpace/3D-Correlation_project/CompareDiffrentResiltsSets/Mw7p1_Lidar/Subsets"
    # # imgList = FileRT.GetFilesBasedOnExtension(path)
    # arrayList = []
    # labels = []
    # bandNumber = 3
    # index = 0
    # report = {"Couple": [], "mean": [], "std": [], "RMSE": []}
    # # path = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/PCA_ICA/ENVI/Dz/ICA/Gaussian/ICA_reconst_4_7/ICA/DZ_residuals.tif"
    # path = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/7p1_3DDA/Sets/Couple_12/1-CCD_correction_classicMethod/NoCCDCorr/Correlation_sameExtent/16SEP08185521_NoCCDCorr_VS_18JUN16213928_NoCCD_corr_detrended_subset.tif"
    # rasterInfo = RT.GetRasterInfo(path)
    # # for img_ in imgList:
    # for band in range(rasterInfo["NbBands"]):
    #     if band == 2:
    #         break
    #     # tempArray1 = RT.ImageAsArray(RT.GetRasterInfo(img_), bandNumber=bandNumber)
    #     tempArray1 = RT.ImageAsArray(rasterInfo, band + 1)
    #     # array = np.ma.masked_where(tempArray1 < -100, tempArray1)
    #     array = np.ma.masked_outside(tempArray1, -4, 4)
    #     arrayList.append(array)
    #     # labels.append(Path(img_).stem.split("-")[1])
    #     labels.append(rasterInfo["BandInfo"][band])
    #     dist_ = Distribution(inputArray=array)
    #     # report["Couple"].append(Path(img_).stem.split("-")[1])
    #     report["Couple"].append(rasterInfo["BandInfo"][band])
    #     report["mean"].append(dist_.mean)
    #     report["std"].append(dist_.std)
    #     report["RMSE"].append(dist_.RMSE)
    #
    # import pandas
    #
    # data = pandas.DataFrame.from_dict(report)
    # print(data)
    # # data.to_csv(os.path.join(path,"band_"+str(bandNumber)+"_Report.csv"), index=False)
    # data.to_csv(os.path.join(os.path.dirname(path), "Residuals_Report.csv"), index=False)
    # PlotDistribution_Batch(inputArrayList=arrayList, labels=labels, nbins=None,
    #                        colors=["gray", "blue", "aqua", "purple", "red", "orange", "peru", "green"],
    #                        title="DZ Residual Density ")
    # # PlotDistribution(inputArray=array2,title="diff_Dz_Couple33 and dDSMs",xlim=[-2,2])
    # ################################################# END ############################################################

    return


if __name__ == '__main__':
    Example()
