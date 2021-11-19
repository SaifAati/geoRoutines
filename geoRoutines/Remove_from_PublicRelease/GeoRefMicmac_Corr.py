from geospatialroutine.Routine import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from geospatialroutine.Plotting.Plotting_Routine import ColorBar


def CorrMicMac_Batch(baseFile, targetFile, outputfolder=[]):
    # Using readlines()
    file1 = open(baseFile, 'r')
    baseImgs = file1.readlines()

    file2 = open(targetFile, 'r')
    targetImg = file2.readlines()
    # print(targetImg)
    for base, target in zip(baseImgs, targetImg):
        base = base.rstrip("\n")
        target = target.rstrip("\n")
        baseName = os.path.basename(base)
        targetName = os.path.basename(target)

        corrName = "CorrMicMac_" + os.path.basename(base)[0:8] + "_vs_" + os.path.basename(target)[0:8] + ".tif"
        rasterOutput = os.path.join(os.path.dirname(base), corrName)
        print(baseName, targetName, rasterOutput)
        CorrMicMac(folderPath=os.path.dirname(base), baseImage=baseName, targetImage=targetName,
                   rasterOutput=rasterOutput, display=False, sameRaster=False, rasterize=True)


    return


def CorrMicMac(folderPath, baseImage, targetImage, SzW=16, reg=0, rasterize=True, sameRaster=False, rasterOutput=[],
               display=False):
    """
    mm3d MM2DPosSism
    [Name=SzW] INT :: {Size of window (Def=4, mean 9x9)}
    [Name=Reg] REAL :: {Regularization (Def=0.3)}


    :param baseImage:
    :param targetImage:
    :return:
    """

    cmd = ["//home/cosicorr/Desktop/micmac/bin/mm3d"]
    cmd.extend(["MM2DPosSism"])
    cmd.extend([baseImage])
    cmd.extend([targetImage])
    cmd.extend(["SzW=" + str(SzW)])
    cmd.extend(["Reg=" + str(reg)])
    call = ""
    for cmd_ in cmd:
        call += cmd_ + " "
    print("call=", call)
    os.chdir(folderPath)
    os.system(call)
    if rasterize:
        px = os.path.join(os.path.join(folderPath, "MEC"), "Px1_Num6_DeZoom1_LeChantier.tif")
        py = os.path.join(os.path.join(folderPath, "MEC"), "Px2_Num6_DeZoom1_LeChantier.tif")
        if os.path.isfile(px) and os.path.isfile(py):
            # print(corrName,rasterOutput)
            GeoResMicMac_Corr(refGrid=os.path.join(folderPath, baseImage), px=px, py=py, outputPath=rasterOutput,
                              threshold=None, display=display, sameRaster=sameRaster)
        else:
            print("Px1_Num6_DeZoom1_LeChantier.tif or Px1_Num6_DeZoom1_LeChantier.tif does not exit!")
        return


def GeoResMicMac_Corr(refGrid, px, py, outputPath=None, threshold=None, display=False, sameRaster=True):
    refGridInfo = GetRasterInfo(refGrid, printInfo=True)
    pxInfo = GetRasterInfo(px, printInfo=True)
    pyInfo = GetRasterInfo(py, printInfo=True)
    ### ----- Verification: all rasters have to have the same diemension  -----###
    if refGridInfo["Dimension"] != pxInfo["Dimension"] or refGridInfo["Dimension"] != pyInfo["Dimension"]:
        print("Error!")
        return
    gsd = np.abs(refGridInfo["Resolution"][0])

    pxArray = gsd * ImageAsArray(imageInfo=pxInfo)
    print("---- loading PXarray")
    pyArray = gsd * ImageAsArray(imageInfo=pyInfo)
    print("---- loading PYarray")
    if threshold:
        pxArray = np.ma.masked_outside(pxArray, -gsd * threshold, gsd * threshold)
        pxArray = pxArray.filled(fill_value=-32654)
        pyArray = np.ma.masked_outside(pyArray, -gsd * threshold, gsd * threshold)
        pyArray = pyArray.filled(fill_value=-32654)

    print(pxArray, pxArray.shape)
    print(pyArray, pxArray.shape)
    if not outputPath:
        outputPath = os.path.join(os.path.dirname(px), "corr_micmac_geoRef.tif")
    if sameRaster:
        WriteRaster(refRasterPath=refGrid, newRasterPath=outputPath, Listarrays=[pxArray, pyArray], numberOfBands=2,
                    descriptions=["East-West", "North-South"], options=["COMPRESS=LZW"], noData=-32654)
    else:
        WriteRaster(refRasterPath=refGrid, newRasterPath=outputPath[0:-4] + "_EW.tif", Listarrays=[pxArray],
                    numberOfBands=1, descriptions=["East-West"], options=["COMPRESS=LZW"], noData=-32654)
        WriteRaster(refRasterPath=refGrid, newRasterPath=outputPath[0:-4] + "_NS.tif", Listarrays=[pyArray],
                    numberOfBands=1, descriptions=["North-South"], options=["COMPRESS=LZW"], noData=-32654)

    if display:
        fig = plt.figure(figsize=(10,10))
        gs = gridspec.GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])  # row,col
        ax2 = fig.add_subplot(gs[0, 1])  # row,col

        my_cmap = plt.get_cmap('rainbow')
        # my_cmap = plt.get_cmap('gist_earth')
        # my_cmap= plt.get_cmap("RdBu")
        # my_cmap = plt.get_cmap('seismic')
        # my_cmap = plt.get_cmap('gray')

        im1 = ax1.imshow(np.ma.masked_values(pxArray, -32654), cmap=my_cmap)  # , vmin=vMin, vmax=vMax, origin="upper")
        # cbar = plt.colorbar(im, label=label)
        ColorBar(ax=ax1, mapobj=im1)
        im2 = ax2.imshow(np.ma.masked_values(pyArray, -32654), cmap=my_cmap)  # , vmin=vMin, vmax=vMax, origin="upper")
        # cbar = plt.colorbar(im, label=label)
        ColorBar(ax=ax2, mapobj=im2)
        plt.show()

    return


if __name__ == '__main__':
    path = "//home/cosicorr/0-WorkSpace/3D-Correlation_project/DEMS_correlation/micmac/med9"
    path2= "//home/cosicorr/0-WorkSpace/3D-Correlation_project/DEMS_correlation/micmac/med9/szw8"
    refGrid = os.path.join(path, "1-Pre_dem_subset_HH_med9.tif")
    px = os.path.join(path2, "Px1_Num6_DeZoom1_LeChantier.tif")
    py = os.path.join(path2, "Px2_Num6_DeZoom1_LeChantier.tif")
    GeoResMicMac_Corr(refGrid, px, py, threshold=None, display=False)

    path = "/home/cosicorr/0-WorkSpace/Temp"
    baseImg = os.path.join(path, "20160221_180512_0c74_3B_AnalyticMS.tif")
    targetImg = os.path.join(path, "20180913_180725_1004_3B_AnalyticMS_SR.tif")
    # CorrMicMac(folderPath=path, baseImage="20160221_180512_0c74_3B_AnalyticMS.tif",
    #            targetImage="20180913_180725_1004_3B_AnalyticMS_SR.tif")
    path = "//home/cosicorr/0-WorkSpace/Landslide_project/Planet/PlanetSubset3/PS_B3_Subsets"
    # CorrMicMac_Batch(baseFile=os.path.join(path, "0-Master.txt"), targetFile=os.path.join(path, "1-Slave.txt"),
    #                  outputfolder=path)
