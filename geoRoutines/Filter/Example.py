# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2020
import matplotlib.animation as animation

from matplotlib import pyplot as plt

from pathlib import Path
from geoRoutines.Plotting.Plotting_Routine import Colorbar_2, ColorBar
import geoRoutines.Filter.Multi_Scale_Filter as MultScaleFilter
import geoRoutines.FilesCommandRoutine as FileRT
import geoRoutines.georoutines as geoRT

def plot(ewImgList, nsImgList):
    fig, (axEw, axNs) = plt.subplots(1, 2, figsize=(15, 10))

    # my_cmap = plt.get_cmap('rainbow')
    # my_cmap = plt.get_cmap('seismic')
    my_cmap = plt.get_cmap('RdYlBu')
    # my_cmap = plt.get_cmap('gray')

    ims = []

    for imgEw_, imgNS_ in zip(ewImgList, nsImgList):
        imageAsArray_EW = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=imgEw_))
        imageAsArray_NS = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=imgNS_))

        imEw = axEw.imshow(imageAsArray_EW, cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
        imNS = axNs.imshow(imageAsArray_NS, cmap=my_cmap, vmin=-5, vmax=5, origin="upper")

        ew_txtTitle = "E/W_" + Path(imgEw_).stem.split("_")[1]
        ns_txtTitle = "N/S_" + Path(imgNS_).stem.split("_")[1]
        txt_ew = axEw.text(800, 50, ew_txtTitle, color="k", bbox=dict(facecolor='white', alpha=0.8))
        txt_ns = axNs.text(800, 50, ns_txtTitle, color="k", bbox=dict(facecolor='white', alpha=0.8))
        # txt = axEw.set_title( txtTitle)
        ims.append([imEw, imNS, txt_ew, txt_ns])

    ColorBar(ax=axNs, mapobj=imNS, label="Displacement [m]")

    ani = animation.ArtistAnimation(fig, ims, interval=300, blit=True, repeat_delay=1000)
    # if saveAnimation == True:
    #     ani.save(animationSavingFile + ".mp4", fps=0.3, dpi=400)  # writer="imagemagick",for GIF
    plt.show()

    return


def plot_withDist(ewImgList, nsImgList, saving=None):
    import matplotlib.animation as animation

    fig = plt.figure(figsize=(20, 15))
    axEw = fig.add_subplot(221)
    axNs = fig.add_subplot(222)
    axEWDist = fig.add_subplot(223)
    axNSDist = fig.add_subplot(224)

    # my_cmap = plt.get_cmap('rainbow')
    # my_cmap = plt.get_cmap('seismic')
    my_cmap = plt.get_cmap('RdYlBu')
    # my_cmap = plt.get_cmap('gray')
    #
    ims = []
    count = 1
    ewData = []
    nsData = []
    ewDataFlat = []
    nsDataFlat = []
    ewTitle = []
    nsTitle = []
    ewObjList = []
    nsObjList = []
    for imgEw_, imgNS_ in zip(ewImgList, nsImgList):
        print("------------", count, "/" + str(len(ewImgList)) + "-----------")
        imageAsArray_EW = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=imgEw_))
        imageAsArray_NS = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=imgNS_))
        dist_EW_obj = RT.Distribution(inputArray=imageAsArray_EW)
        dist_NS_obj = RT.Distribution(inputArray=imageAsArray_NS)
        ewObjList.append(dist_EW_obj)
        nsObjList.append(dist_NS_obj)
        ewData.append(imageAsArray_EW)
        nsData.append(imageAsArray_NS)
        ewDataFlat.append(imageAsArray_EW.flatten())
        nsDataFlat.append(imageAsArray_NS.flatten())

        ew_txtTitle = "E/W_" + Path(imgEw_).stem.split("_")[1]
        ns_txtTitle = "N/S_" + Path(imgNS_).stem.split("_")[1]
        ewTitle.append(ew_txtTitle)
        nsTitle.append(ns_txtTitle)

        count += 1

    imEw = axEw.imshow(ewData[0], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
    imNS = axNs.imshow(nsData[0], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
    axEw.set_title(ewTitle[0])
    axNs.set_title(nsTitle[0])
    axNs.axis('off')
    axEw.axis('off')
    ColorBar(ax=axNs, mapobj=imNS, label="Displacement [m]")
    kwargs = dict(histtype='stepfilled', alpha=0.3, bins=200)

    ewDist = axEWDist.hist(ewDataFlat[0], **kwargs)
    nsDist = axNSDist.hist(nsDataFlat[0], **kwargs)

    axEWDist.set_xlim(-75, 75)
    axEWDist.set_ylim(0, 300000)
    axNSDist.set_xlim(-75, 75)
    axNSDist.set_ylim(0, 300000)
    ew_txt = "mean=" + str(ewObjList[0].mean) + "\n" + "sigma=" + str(ewObjList[0].sigma_) + "\n" + "min=" + str(
        ewObjList[0].min)+"\n" + "max=" + str(ewObjList[0].max)
    axEWDist.text(20, axEWDist.get_ylim()[1] - 70000, ew_txt, color="k", bbox=dict(facecolor='white', alpha=0.8),
                  fontsize=14)
    ns_txt = "mean=" + str(nsObjList[0].mean) + "\n" + "sigma=" + str(nsObjList[0].sigma_)+ "\n" + "min=" + str(
        nsObjList[0].min)+"\n" + "max=" + str(nsObjList[0].max)
    axNSDist.text(20, axEWDist.get_ylim()[1] - 70000, ns_txt, color="k", bbox=dict(facecolor='white', alpha=0.8),
                  fontsize=14)

    def updatefig(i):
        print("Frame =", i + 1)
        axNSDist.cla()
        axEWDist.cla()
        axEw.cla()
        axNs.cla()

        imEw = axEw.imshow(ewData[i], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
        imNS = axNs.imshow(nsData[i], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
        axEw.set_title(ewTitle[i])
        axNs.set_title(nsTitle[i])
        axNs.axis('off')
        axEw.axis('off')

        ewDist = axEWDist.hist(ewDataFlat[i], **kwargs)
        nsDist = axNSDist.hist(nsDataFlat[i], **kwargs)
        axEWDist.set_xlim(-75, 75)
        axEWDist.set_ylim(0, 300000)
        axNSDist.set_xlim(-75, 75)
        axNSDist.set_ylim(0, 300000)
        ew_txt = "mean=" + str(ewObjList[i].mean) + "\n" + "sigma=" + str(ewObjList[i].sigma_)+ "\n" + "min=" + str(
        ewObjList[i].min)+"\n" + "max=" + str(ewObjList[i].max)
        axEWDist.text(20, axEWDist.get_ylim()[1] - 70000, ew_txt, color="k", bbox=dict(facecolor='white', alpha=0.8),
                      fontsize=14)
        ns_txt = "mean=" + str(nsObjList[i].mean) + "\n" + "sigma=" + str(nsObjList[i].sigma_)+"\n" + "min=" + str(
        nsObjList[i].min)+"\n" + "max=" + str(nsObjList[i].max)
        axNSDist.text(20, axEWDist.get_ylim()[1] - 70000, ns_txt, color="k", bbox=dict(facecolor='white', alpha=0.8),
                      fontsize=14)
        return imEw, imNS, ewDist, nsDist

    ani = animation.FuncAnimation(fig, updatefig, frames=len(ewData), interval=300, repeat_delay=3000)
    plt.show()

    #
    # ani = animation.ArtistAnimation(fig, ims, interval=300, blit=True, repeat_delay=1000)
    if saving:
        # Set up formatting for the movie files
        print(saving)
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=1, metadata=dict(artist='saif@caltech.edu'), bitrate=7800)
        # line_ani.save('lines.mp4', writer=writer) 
        # ani.save(saving ,writer="imagemagick")#,for GIF
        ani.save(saving, writer)  # dpi=100, writer="imagemagick")
        # plt.show()


def plot_beforeandafter(ewImgList, nsImgList, ewListAfter, nsListafter, saving=None):
    import matplotlib.animation as animation

    fig = plt.figure(figsize=(20, 15))
    axEw = fig.add_subplot(221)
    axNs = fig.add_subplot(222)
    axEw_after = fig.add_subplot(223)
    axNs_after = fig.add_subplot(224)

    # my_cmap = plt.get_cmap('rainbow')
    # my_cmap = plt.get_cmap('seismic')
    my_cmap = plt.get_cmap('RdYlBu')
    # my_cmap = plt.get_cmap('gray')
    #
    ims = []
    count = 1
    ewData = []
    nsData = []
    ewDataFlat = []
    nsDataFlat = []
    ewTitle = []
    nsTitle = []
    ewObjList = []
    nsObjList = []
    ewData_after = []
    nsData_after = []
    ewDataFlat_after = []
    nsDataFlat_after = []
    ewTitle_after = []
    nsTitle_after = []
    ewObjList_after = []
    nsObjList_after = []
    for imgEw_, imgNS_, imgEw_after, imgNS_after in zip(ewImgList, nsImgList, ewListAfter, nsListafter):
        print("------------", count, "/" + str(len(ewImgList)) + "-----------")
        imageAsArray_EW = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=imgEw_))
        imageAsArray_NS = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=imgNS_))
        imageAsArray_EW_after = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=imgEw_after))
        imageAsArray_NS_after = RT.ImageAsArray(RT.GetRasterInfo(inputRaster=imgNS_after))
        # dist_EW_obj = RT.Distribution(inputArray=imageAsArray_EW)
        # dist_NS_obj = RT.Distribution(inputArray=imageAsArray_NS)
        # dist_EW_obj_after = RT.Distribution(inputArray=imageAsArray_EW_after)
        # dist_NS_obj_after = RT.Distribution(inputArray=imageAsArray_NS_after)
        # ewObjList.append(dist_EW_obj)
        # nsObjList.append(dist_NS_obj)
        # nsObjList_after.append(dist_NS_obj_after)
        # nsObjList_after.append(dist_NS_obj_after)
        ewData.append(imageAsArray_EW)
        nsData.append(imageAsArray_NS)
        ewData_after.append(imageAsArray_EW_after)
        nsData_after.append(imageAsArray_NS_after)
        # ewDataFlat_after.append(imageAsArray_NS_after.flatten())
        # nsDataFlat_after.append(imageAsArray_NS_after.flatten())

        ew_txtTitle = "E/W_" + Path(imgEw_).stem.split("_")[1]
        ns_txtTitle = "N/S_" + Path(imgNS_).stem.split("_")[1]
        ew_txtTitle_after = "E/W_" + Path(imgEw_after).stem.split("_")[1] + "_filterd"
        ns_txtTitle_after = "N/S_" + Path(imgNS_after).stem.split("_")[1] + "_filtered"
        ewTitle.append(ew_txtTitle)
        nsTitle.append(ns_txtTitle)
        ewTitle_after.append(ew_txtTitle_after)
        nsTitle_after.append(ns_txtTitle_after)

        count += 1

    imEw = axEw.imshow(ewData[0], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
    imNS = axNs.imshow(nsData[0], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
    imEw_after = axEw_after.imshow(ewData_after[0], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
    imNS_after = axNs_after.imshow(nsData_after[0], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
    axEw.set_title(ewTitle[0])
    axNs.set_title(nsTitle[0])
    axEw.set_title(ewTitle[0])
    axNs.set_title(nsTitle[0])
    axNs.axis('off')
    axEw.axis('off')
    axNs_after.axis('off')
    axEw_after.axis('off')
    ColorBar(ax=axNs, mapobj=imNS, label="Displacement [m]")
    ColorBar(ax=axNs_after, mapobj=imNS, label="Displacement [m]")

    def updatefig(i):
        print("Frame =", i + 1)

        axEw.cla()
        axNs.cla()
        axNs_after.cla()
        axNs_after.cla()

        imEw = axEw.imshow(ewData[i], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
        imNS = axNs.imshow(nsData[i], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
        imEw_after = axEw_after.imshow(ewData_after[i], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
        imNS_after = axNs_after.imshow(nsData_after[i], cmap=my_cmap, vmin=-5, vmax=5, origin="upper")
        axEw.set_title(ewTitle[i])
        axNs.set_title(nsTitle[i])
        axEw_after.set_title(ewTitle_after[i])
        axNs_after.set_title(nsTitle_after[i])
        axNs.axis('off')
        axEw.axis('off')
        axNs_after.axis('off')
        axEw_after.axis('off')

        return imEw, imNS, imEw_after, imNS_after

    ani = animation.FuncAnimation(fig, updatefig, frames=len(ewData), interval=300, repeat_delay=3000)
    plt.show()

    #
    # ani = animation.ArtistAnimation(fig, ims, interval=300, blit=True, repeat_delay=1000)
    if saving:
        # Set up formatting for the movie files
        print(saving)
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=1, metadata=dict(artist='saif@caltech.edu'), bitrate=5800)
        # line_ani.save('lines.mp4', writer=writer)
        # ani.save(saving ,writer="imagemagick")#,for GIF
        ani.save(saving, writer)  # dpi=100, writer="imagemagick")
        # plt.show()


def Filtering(inputList, snrList, outputFolder):
    filter = MultScaleFilter.Filtering()
    filter.MultiScaleFilterWithSNR_(inputData=inputList, outputFolder=outputFolder, snrData=snrList,
                                    thresholdVelovity=15, thresholdMedian=1, step=5,
                                    multiScaleFactor=3,
                                    withSNR=True)

    return


if __name__ == '__main__':
    folderPath = "//home/cosicorr/0-WorkSpace/7day-increments/Original"
    folderPath ="/home/cosicorr/0-WorkSpace/7day-increments/results/1-filtered_1"

    # _2019-01-22-testf2_band02__2019-01-22-testf2_band02_filtered
    ewList = FileRT.GetFilesBasedOnExtension(path=folderPath, filter="*_band01_filtered.tif")
    nsList = FileRT.GetFilesBasedOnExtension(path=folderPath, filter="*_band02_filtered.tif")
    # snrList = FileRT.GetFilesBasedOnExtension(path=folderPath, filter="*_band03.tif")
    path = "/home/cosicorr/0-WorkSpace/7day-increments/results/2-filtered_2"
    ewListAfter = FileRT.GetFilesBasedOnExtension(path=path, filter="*_band01_filtered.tif")
    nsListAfter = FileRT.GetFilesBasedOnExtension(path=path, filter="*_band02_filtered.tif")
    # snrList = FileRT.GetFilesBasedOnExtension(path=folderPath, filter="*_band03.tif")
    plot_beforeandafter(ewList, nsList, ewListAfter, nsListAfter,
                        "//home/cosicorr/0-WorkSpace/7day-increments/results/Anim_filtered1_&_filtered2.mp4")
    # path = "/home/cosicorr/Desktop/0-WorkSpace/Temp/TempSandune/filtered"
    # ewList = FileRT.GetFilesBasedOnExtension(path)
    # plot_withDist(ewListAfter, nsListAfter, "/home/cosicorr/Desktop/0-WorkSpace/Temp/TempSandune/2-Anim_After.mp4")
    # plot_withDist(ewList, nsList, "/home/cosicorr/Desktop/0-WorkSpace/Temp/TempSandune/EW_anim_filtered.mp4")
    # plot(ewList,nsList)
    # Filtering(ewList, snrList, "/home/cosicorr/0-WorkSpace/7day-increments/results/2-filtered_2")
