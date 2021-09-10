import os

from FilesCommandRoutine import GetFilesBasedOnExtension
import Main_routine as RT
import gdal
import numba


def Convert_NTF_2_TIFF(inputNTF, outputPath=None, bandsList=None):
    """

    :param inputNTF:
    :param bandsList: list of bands to be extracted, default all bands : int list
    :return:
    """
    ntfRasterInfo = RT.GetRasterInfo(inputRaster=inputNTF)
    raster = ntfRasterInfo["Raster"]

    metaData = ntfRasterInfo["MetaData"]
    rpc = ntfRasterInfo["RPC"]
    if not outputPath:
        newRasterPath = os.path.join(os.path.dirname(inputNTF), os.path.basename(inputNTF)[:-4] + ".tif")
        print("newRasterPath=", newRasterPath)
    else:
        newRasterPath = os.path.join(outputPath, os.path.basename(inputNTF)[:-4] + ".tif")
        print("newRasterPath=", newRasterPath)

    arrays = []
    if not bandsList:
        for i in range(ntfRasterInfo["NbBands"]):
            band = raster.GetRasterBand(i + 1)
            if band.GetOverviewCount() > 0:
                print("Band has {} overviews".format(band.GetOverviewCount()))
            if band.GetRasterColorTable():
                print("Band has a color table with {} entries".format(band.GetRasterColorTable().GetCount()))
            arrays.append(RT.ImageAsArray(imageInfo=ntfRasterInfo, bandNumber=i + 1))

        if raster.GetGCPs():
            RT.WriteRaster(refRasterPath=inputNTF,
                           newRasterPath=newRasterPath, Listarrays=arrays,
                           numberOfBands=ntfRasterInfo["NbBands"], metaData=metaData, RPC=rpc,
                           GCPs=[list(raster.GetGCPs()), raster.GetGCPProjection()])
        else:
            RT.WriteRaster(refRasterPath=inputNTF,
                           newRasterPath=newRasterPath, Listarrays=arrays,
                           numberOfBands=ntfRasterInfo["NbBands"], metaData=metaData, RPC=rpc)

    else:
        for index, band_ in enumerate(bandsList):
            print("Converting Band number: ", band_)

            band = raster.GetRasterBand(index + 1)
            if band.GetOverviewCount() > 0:
                print("Band has {} overviews".format(band.GetOverviewCount()))
                for over in range(band.GetOverviewCount()):
                    print(band.GetOverview(over))
            if band.GetRasterColorTable():
                print("Band has a color table with {} entries".format(band.GetRasterColorTable().GetCount()))

            arrays.append(RT.ImageAsArray(imageInfo=ntfRasterInfo, bandNumber=band_))

        # RT.WriteRaster(refRasterPath=inputNTF,
        #                    newRasterPath=newRasterPath, Listarrays=arrays,
        #                    numberOfBands=len(bandsList))

        if raster.GetGCPs():
            RT.WriteRaster(refRasterPath=inputNTF,
                           newRasterPath=newRasterPath, Listarrays=arrays,
                           numberOfBands=len(bandsList), metaData=metaData, RPC=rpc,
                           GCPs=[list(raster.GetGCPs()), raster.GetGCPProjection()])
        else:
            RT.WriteRaster(refRasterPath=inputNTF,
                           newRasterPath=newRasterPath, Listarrays=arrays,
                           numberOfBands=len(bandsList), metaData=metaData, RPC=rpc)

    # band.GetOverviewCount()
    return


def Convert_NTF_2_TIFF_V2(inputNTF, outputPath=None):
    if not outputPath:
        newRasterPath = os.path.join(os.path.dirname(inputNTF), os.path.basename(inputNTF)[:-4] + ".tif")
        print("newRasterPath=", newRasterPath)
    else:
        newRasterPath = os.path.join(outputPath, os.path.basename(inputNTF)[:-4] + ".tif")
        print("newRasterPath=", newRasterPath)
    srcDS = gdal.Open(inputNTF)
    gdal.Translate(newRasterPath, srcDS, format="GTiff",outputType=gdal.GDT_Float64)

    return


def ConvertTIF_2_NTF(inputTifFile, outputPath=None):
    if not outputPath:
        newRasterPath = os.path.join(os.path.dirname(inputTifFile), os.path.basename(inputTifFile)[:-4] + ".ntf")
        print("newRasterPath=", newRasterPath)
    else:
        newRasterPath = os.path.join(outputPath, os.path.basename(inputTifFile)[:-4] + ".ntf")
        print("newRasterPath=", newRasterPath)
        srcDS = gdal.Open(inputTifFile)
        gdal.Translate(newRasterPath, srcDS, format="nitf")

    return


# @numba.jit( parallel=True)
def BatchConvert(inputPath, filter="*.tif", outputPath=None, mode=1):
    """

    :param inputPath:
    :param filter:
    :param outputPath:
    :param mode: 1 Tif 2 NTF
                 2 NTF 2 tif
    :return:
    """

    if mode == 1:
        filesList = GetFilesBasedOnExtension(path=inputPath)
        for file in filesList:
            ConvertTIF_2_NTF(inputTifFile=file, outputPath=outputPath)
    if mode == 2:
        filesList = GetFilesBasedOnExtension(path=inputPath, filter="*.NTF")
        for i in range( len(filesList)):
            Convert_NTF_2_TIFF_V2(inputNTF=filesList[i], outputPath=outputPath)

    return


if __name__ == '__main__':
    # path = "//home/cosi/2-Data/4-Shisper/5-WV/Shisper/2017_111_stereo_M1_BS/"
    # ntif = path + "05MAY17WV020500017MAY05054910-M1BS_R3C1-500615038050_04_P001____GA_E0AAAAAAIAAO0.NTF"
    # # print(os.path.dirname(ntif))
    # # print(os.path.basename(ntif))
    # # Convert_NTF_2_TIFF(inputNTF=ntif,
    # #                    outputPath="/home/cosi/2-Data/4-Shisper/5-WV/Shisper/2017_111_stereo_M1_BS/Correlation",
    # #                    bandsList=[5])
    #
    # path = "//media/storage/Saif/1-Shishper/WV-Data/1-SortedData_WV/2019/2019-M1BS/NTF/"
    # files = []
    # # r=root, d=directories, f = files
    # for r, d, f in os.walk(path):
    #     for file in f:
    #         if '.NTF' in file:
    #             files.append(os.path.join(r, file))
    #
    # for imgPath in files:
    #     print("img:", imgPath)
    #     # Convert_NTF_2_TIFF(inputNTF=imgPath)
    #     Convert_NTF_2_TIFF_V2(inputNTF=imgPath,
    #                           outputPath="/media/storage/Saif/1-Shishper/WV-Data/1-SortedData_WV/2019/2019-M1BS/TIF/")

    path = "/home/cosi/2-Data/3-3DProject/Shisper_3dApproach"
    BatchConvert(inputPath=path, mode=2)
