import glob
import os
import zipfile
from pathlib import Path
from xml.dom import minidom

import dateutil.parser as dparser
import numpy as np
from osgeo import gdal

import Correlation_misc as CorrMisc
import FilesCommandRoutine as FileRT
import Main_routine as RT


class Sentinel:
    def __init__(self):
        self.inputDataPath_SL = ""
        self.savingPath = ""
        self.windDim = []
        self.band2Subset = "B04"

        self.pathMasterImg = ""
        self.pathSlaveImg = "",
        self.corrOption = ["frequency", "32", "32", "32", "32", "6", "6", "4", "0.9", "1", "1"]
        self.outputSavingPath = ""

        self.corrFolder = ""
        self.datesFile = ""

    def UnzipSentinel(self):
        """
        :param slRawDataDirectoryPAth:
        :param savingPath:
        :param folderName:
        :return:
        """
        ## Create the SL_UnzipedData_Folder
        unzipedDataPath_SL = self.savingPath

        folders = FileRT.FilesInDirectory(path=self.inputDataPath_SL)
        for folder_ in folders:
            if ".zip" in folder_:
                filePath = self.inputDataPath_SL + folder_  # Complete file path of each file in the directory
                print(filePath)
                zip_ref = zipfile.ZipFile(filePath, 'r')
                title = folder_[:-4]
                print(title)
                zip_ref.extractall(path=unzipedDataPath_SL + title + "/")
                zip_ref.close()
        return

    def GetSelectedBandInfo(self, unzippedS2Path, bandNumber="B04"):
        bandInfo = {}
        if os.path.basename(unzippedS2Path).startswith('S2') or os.path.basename(unzippedS2Path).startswith('L1C'):
            bandInfo['Satellite'] = "S2"
            print("unzippedS2Path:", unzippedS2Path)
            safePath = glob.glob(os.path.join(unzippedS2Path, "*.SAFE"))
            if not safePath:
                raise Exception('Could not find any SAFE file in ' + unzippedS2Path)
            else:
                safePath = safePath[0]
                # print("safePath:",safePath)
                bandInfo["Platform"] = self.OpenMTD_MSIL1C_xml(safePath=safePath)

            granulePath = glob.glob(os.path.join(safePath, 'GRANULE'))
            if not granulePath:
                raise Exception('Could not find any GRANULE file in ' + unzippedS2Path)
            else:
                granulePath = glob.glob(os.path.join(granulePath[0], "*"))[0]

            imgDataPath = glob.glob(os.path.join(granulePath, "IMG_DATA"))
            if not imgDataPath:
                raise Exception('Could not find any IMG_DATA file in ' + granulePath)
            else:
                imgDataPath = imgDataPath[0]

            bandsPaths = glob.glob(os.path.join(imgDataPath, "*.jp2"))
            bandsPaths.sort()
            for band_ in bandsPaths:
                if os.path.basename(band_).endswith(bandNumber + ".jp2"):
                    matchBandPath = band_
                    bandInfo["Number"] = bandNumber
                    bandInfo["Path"] = matchBandPath

            bandInfo["Date"], bandInfo["Time"] = self.OpenMTDxml(granulePath=granulePath)
            # get gml path
            qiPath = os.path.join(granulePath, 'QI_DATA')
            gmlPath = glob.glob(os.path.join(qiPath, '*DETFOO*' + bandNumber + '.gml'))
            if not gmlPath:
                print("gmlPath not found")
            else:
                bandInfo['Azimuth'] = gmlPath[0]


        else:
            raise Exception('Could not identify a valid S-2 product.')

        return bandInfo

    def OpenMTD_MSIL1C_xml(self, safePath):
        S2Platform = "S2A"
        xmlPath = glob.glob(os.path.join(safePath, '*.xml'))
        for xmlFilePath in xmlPath:
            if "MTD_MSIL1C" in os.path.basename(xmlFilePath):
                # parse an xml file by name
                mydoc = minidom.parse(xmlFilePath)
                items = mydoc.getElementsByTagName('SPACECRAFT_NAME')

                for elem in items:
                    S2Platform = elem.firstChild.data
        if S2Platform == "Sentinel-2B":
            S2Platform = "S2B"
        if S2Platform == "Sentinel-2A":
            S2Platform = "S2A"
        return S2Platform

    def OpenMTDxml(self, granulePath):
        xmlPath = glob.glob(os.path.join(granulePath, '*.xml'))[0]
        from xml.dom import minidom
        # parse an xml file by name
        mydoc = minidom.parse(xmlPath)
        items = mydoc.getElementsByTagName('SENSING_TIME')
        for elem in items:
            date = elem.firstChild.data

        acquisitionDate = date.split("T")[0]
        acquisitionTime = date.split("T")[1]

        # import xml.etree.ElementTree
        # e = xml.etree.ElementTree.parse(xmlPath).getroot()
        # acquisition_time = e[0][3].text
        # acquisition_time[0:10]
        # acquisition_time[12:]
        return acquisitionDate, acquisitionTime

    def SubsetSentinel(self, inputFolder, savingFolder):
        # "gdal_translate -projwin 449230.87415403477 4034940.424743298 472437.7298007704 4015846.1764263636 -of GTiff -co NBITS=16
        # gdalTranslateCommand = "gdal_translate -projwin %s %s %s %s -of GTiff -co NBITS=16 " % (
        #     str(self.windDim[0]), str(self.windDim[1]),
        #     str(self.windDim[2]), str(self.windDim[3]))
        ## Create folder to save SL_subsets
        subsetSavingFolder_S2 = FileRT.CreateDirectory(directoryPath=savingFolder,
                                                       folderName="S2_%s_Subsets" % (self.band2Subset))
        ## Create DateMetadata_Subset File
        dateMetadataSubset_S2 = open(subsetSavingFolder_S2 + "SubsetDates.txt", "w+")

        unzippedS2 = FileRT.FilesInDirectory(path=inputFolder)

        for unzippedS2Folder in unzippedS2:
            bandInfo = self.GetSelectedBandInfo(unzippedS2Path=os.path.join(inputFolder, unzippedS2Folder),
                                                bandNumber=self.band2Subset)
            print("Selected band info", bandInfo)
            bandPath = bandInfo["Path"]
            subsetTitle = bandInfo["Date"] + "_" + bandInfo["Platform"] + "_" + self.band2Subset + "_Subset.tif"

            metaData = []
            for key, value in bandInfo.items():
                temp = [key, value]
                metaData.append(temp)

            RT.SubsetImg(inputImg=bandPath, outputPath=os.path.join(subsetSavingFolder_S2, subsetTitle),
                         windDim=self.windDim, metaData=metaData)

            dateMetadataSubset_S2.write(
                "SubsetId = " + os.path.basename(bandPath) + " Date =" + bandInfo["Date"] + "\r\n")
        return

    def CreateBatchCorrelationFiles(self, corrStrategy=1):
        """
        :param pathMasterImg:
        :param pathSlaveImg:
        :param corrOption:
        :return:
        """
        CorrMisc.CreateBatchCorrelationFiles(pathMasterImgs=self.pathMasterImg, pathSlaveImgs=self.pathSlaveImg,
                                             savingPath=self.savingPath, corrOptions=self.corrOption,
                                             corrStrategy=corrStrategy)
        return

    def ExtractEWNSFromCorrelation(self, savingBandsPath):
        files = FileRT.FilesInDirectory(self.corrFolder)
        files.sort()
        corrFilesPath = []

        pathEW = FileRT.CreateDirectory(directoryPath=savingBandsPath, folderName="EW")
        pathNS = FileRT.CreateDirectory(directoryPath=savingBandsPath, folderName="NS")

        for file_ in files:
            if 'Corr' in file_ and '.hdr' not in file_ and ".txt" not in file_ and ".xml" not in file_:
                corrFilesPath.append(file_)

        for corrFile_ in corrFilesPath:
            ## Split the corr file name to retrive the name of master and slave img
            corrFile_Temp = corrFile_.replace("Corr_", "")
            corrFile_TempList = corrFile_Temp.split("_vs_")
            masterImg = corrFile_TempList[0].split(".tif")[0]
            slaveImg = corrFile_TempList[1].split(".tif")[0]
            ## open the date metadata to fecth the date of each image
            dateMasterString = masterImg.split("_")[0]
            dateMaster = dparser.parse(dateMasterString, fuzzy=True)
            print("Master=", dateMaster)

            dateSlaveString = slaveImg.split("_")[0]
            dateSlave = dparser.parse(dateSlaveString, fuzzy=True)
            print("Slave=", dateSlave)

            corr = gdal.Open(self.corrFolder + corrFile_)
            EW = corr.GetRasterBand(1)
            NS = corr.GetRasterBand(2)
            SNR = corr.GetRasterBand(3)
            ewArray = EW.ReadAsArray()
            nsArray = NS.ReadAsArray()
            snrArray = SNR.ReadAsArray()
            # velocityArray1 = np.sqrt(ewArray ** 2 + nsArray ** 2) / daysDiff.days

            titleEW = "Master_" + dateMasterString + "_vs_Slave_%s_SL_%s_EW.tif" % (
                dateSlaveString, masterImg.split("_")[2])
            titlePathEW = pathEW + titleEW

            RT.Array2MultiBandRaster(refRasterPath=self.corrFolder + corrFile_, newRasterPath=titlePathEW,
                                     arrays=[ewArray, snrArray], numberOfBands=2,
                                     descriptions=[EW.GetDescription(), SNR.GetDescription()])

            titleNS = "Master_" + dateMasterString + "_vs_Slave_%s_SL_%s_NS.tif" % (
                dateSlaveString, masterImg.split("_")[2])
            titlePathNS = pathNS + titleNS

            RT.Array2MultiBandRaster(refRasterPath=self.corrFolder + corrFile_, newRasterPath=titlePathNS,
                                     arrays=[nsArray, snrArray], numberOfBands=2,
                                     descriptions=[NS.GetDescription(), SNR.GetDescription()])

        return

    def ComputeEWNSVelocity(self, corrFolder, outputSavingPath, withSNR=True):
        """
        : corr file should be defined
        :param savingBandsPath:
        :return:
        """

        self.corrFolder = corrFolder
        files = FileRT.FilesInDirectory(self.corrFolder)
        files.sort()
        corrFilesPath = []

        pathEW = FileRT.CreateDirectory(directoryPath=outputSavingPath, folderName="EW_Velocity")
        pathNS = FileRT.CreateDirectory(directoryPath=outputSavingPath, folderName="NS_Velocity")

        for file_ in files:
            if 'Corr' in file_ and '.hdr' not in file_ and ".txt" not in file_ and ".xml" not in file_:
                corrFilesPath.append(file_)

        for corrFile_ in corrFilesPath:
            corrFilePath = os.path.join(corrFolder, corrFile_)
            velocityInfo = CorrMisc.ComputeEWNSVelocity(corrFilePath=corrFilePath)

            if withSNR:

                RT.WriteRaster(refRasterPath=corrFilePath,
                               newRasterPath=os.path.join(pathEW, velocityInfo["EW"]["FileName"]),
                               Listarrays=[velocityInfo["EW"]["Array"], velocityInfo["SNR"]["Array"]], numberOfBands=2,
                               descriptions=[velocityInfo["EW"]["Description"], velocityInfo["SNR"]["Description"]],
                               metaData=velocityInfo)

                RT.WriteRaster(refRasterPath=corrFilePath,
                               newRasterPath=os.path.join(pathNS, velocityInfo["NS"]["FileName"]),
                               Listarrays=[velocityInfo["NS"]["Array"], velocityInfo["SNR"]["Array"]], numberOfBands=2,
                               descriptions=[velocityInfo["NS"]["Description"], velocityInfo["SNR"]["Description"]],
                               metaData=velocityInfo)
            else:
                RT.WriteRaster(refRasterPath=corrFilePath,
                               newRasterPath=os.path.join(pathEW, velocityInfo["EW"]["FileName"]),
                               Listarrays=[velocityInfo["EW"]["Array"]], numberOfBands=1,
                               descriptions=[velocityInfo["EW"]["Description"]],
                               metaData=velocityInfo)

                RT.WriteRaster(refRasterPath=corrFilePath,
                               newRasterPath=os.path.join(pathEW, velocityInfo["NS"]["FileName"]),
                               Listarrays=[velocityInfo["NS"]["Array"]], numberOfBands=1,
                               descriptions=[velocityInfo["NS"]["Description"]],
                               metaData=velocityInfo)

        return

    def VelocityMapComputing(self, boolSNRband=True):
        """
        This function compute the velocity map between NS and EW of the correlation map
        :param corrPath: correlation folder path : string
        :param savingDirectoryPath: string
        :return:
        """
        files = FileRT.FilesInDirectory(self.corrFolder)
        files.sort()
        corrFilesPath = []
        for file_ in files:
            if 'Corr' in file_ and '.hdr' not in file_ and ".txt" not in file_ and ".xml" not in file_:
                corrFilesPath.append(file_)

        for corrFile_ in corrFilesPath:
            ## Split the corr file name to retrive the name of master and slave img
            corrFile_Temp = corrFile_.replace("Corr_", "")
            corrFile_TempList = corrFile_Temp.split("_vs_")
            masterImg = corrFile_TempList[0].split(".tif")[0]
            slaveImg = corrFile_TempList[1].split(".tif")[0]
            ## open the date metadata to fecth the date of each image
            dateMasterString = masterImg.split("_")[0]
            dateMaster = dparser.parse(dateMasterString, fuzzy=True)
            print("Master=", dateMaster)

            dateSlaveString = slaveImg.split("_")[0]
            dateSlave = dparser.parse(dateSlaveString, fuzzy=True)
            print("Slave=", dateSlave)

            ## Compute the delay between the tow images
            daysDiff = abs(dateMaster - dateSlave)
            print("daysDiff=", daysDiff.days)

            ## Create a velocity map
            title = "Master_" + dateMasterString + "_vs_Slave_%s_SL_%s_velocity.tif" % (
                dateSlaveString, masterImg.split("_")[2])
            print("title=", title)
            titlePath = self.savingPath + title

            # command = "gdal_calc.py --calc=""\"sqrt(A*A+B*B)/%s\" -A %s --A_band=1 -B %s --B_band=2 --outfile=%s" % (
            #     str(daysDiff.days), self.corrFolder + corrFile_, self.corrFolder + corrFile_, titlePath)

            # os.system(command)
            corr = gdal.Open(self.corrFolder + corrFile_)
            EW = corr.GetRasterBand(1)
            NS = corr.GetRasterBand(2)
            SNR = corr.GetRasterBand(3)
            ewArray = EW.ReadAsArray()
            ewArray = np.array(ewArray, dtype=np.float64)
            nsArray = NS.ReadAsArray()
            nsArray = np.array(nsArray, dtype=np.float64)
            snrArray = SNR.ReadAsArray()
            snrArray = np.array(snrArray, dtype=np.float64)
            velocityArray1 = np.sqrt(ewArray * ewArray + nsArray * nsArray) / daysDiff.days
            if boolSNRband == False:
                RT.Array2MultiBandRaster(refRasterPath=self.corrFolder + corrFile_, newRasterPath=titlePath,
                                         arrays=[velocityArray1], numberOfBands=1,
                                         descriptions=["Velocity"])
            if boolSNRband == True:
                RT.Array2MultiBandRaster(refRasterPath=self.corrFolder + corrFile_, newRasterPath=titlePath,
                                         arrays=[velocityArray1, snrArray], numberOfBands=2,
                                         descriptions=["Velocity", SNR.GetDescription()])

        return

    def ApplyROIMask(self, inputFolder, mask_Path, outputFolder="", setNaN=True):
        """
        Apply mask to map
        :param InputDirectory: images path : string
        :param outputDirectroy: results output path : string
        :param mask_Path: path of the mask file : string
        :return:
        """
        if not outputFolder:
            workspacePath = Path(inputFolder).parents[1]
            folders = FileRT.FilesInDirectory(path=workspacePath)
            if "MaskedData" in folders:
                tempPath = os.path.join(workspacePath, "MaskedData")
            else:
                tempPath = FileRT.CreateDirectory(directoryPath=workspacePath, folderName="MaskedData")

            outputFolder = FileRT.CreateDirectory(directoryPath=tempPath,
                                                  folderName=os.path.basename(
                                                      os.path.normpath(inputFolder)) + "_Masked")

        maskRaster = gdal.Open(mask_Path)
        maskArray = maskRaster.GetRasterBand(1).ReadAsArray()
        maskArray = maskArray.astype(float)
        if setNaN == True:
            maskArray[maskArray == 0] = np.nan

        files = FileRT.FilesInDirectory(inputFolder)

        for file_ in files:
            if ".enp" not in file_ and ".hdr" not in file_:
                print("Masking:", file_)
                arrays = []
                descriptions = []
                rasterInfo = RT.GetRasterInfo(inputFolder + file_)
                raster = rasterInfo["Raster"]
                for i in range(raster.RasterCount):
                    bandNumber = i + 1
                    bandTmp = raster.GetRasterBand(bandNumber)
                    bandTmpArr = RT.ImageAsArray(rasterInfo, bandNumber=bandNumber)
                    arrays.append(bandTmpArr * maskArray)
                    descriptions.append(bandTmp.GetDescription())

                title = outputFolder + file_[:-4] + "_Masked.tif"
                metaData = rasterInfo["MetaData"]
                metaData["MaskPath"] = mask_Path
                RT.WriteRaster(refRasterPath=inputFolder + file_, newRasterPath=title,
                               Listarrays=arrays, numberOfBands=raster.RasterCount, descriptions=descriptions,
                               metaData=metaData)

        return

    def ResizeSentinel2Landsat(self, landsatRefImg, inputPath, savingPath):
        # EW component
        files = FileRT.FilesInDirectory(path=inputPath)
        for file_ in files:
            RT.ResizeRaster(refImg=landsatRefImg, inputImg=inputPath + file_, outputPath=savingPath + file_)

    def DeleteSubsets(self, pathRefFile, subsetImgPath):
        ##### Version 0.01 ##### To be continued
        """
        this function will delete the subset that does not exit in the refFile
        :param pathRefFile:
        :return:
        """
        ## Read all imgs in subset folder
        files = FileRT.FilesInDirectory(path=subsetImgPath)

        for file_ in files:
            if ".txt" in file_:
                subsetDateFilePath = subsetImgPath + file_
        print(subsetDateFilePath)

        for file_ in files:
            if ".tif" in file_:
                ## fetch date of the img
                dateFile = open(subsetDateFilePath, "r")
                for line in dateFile:
                    if file_ in line:
                        dateStr = line.split("=")[2]

                        if dateStr not in open(pathRefFile).read():
                            print(False)
                            print(dateStr, file_)
                            os.remove(subsetImgPath + file_)

    def RenameWithDate(self, refTxtFile, imgPath, savingPath, prefix=""):
        srcFiles = FileRT.FilesInDirectory(path=imgPath)
        print(srcFiles)
        i = 1
        for file_ in srcFiles:
            if ".tif" in file_:
                ## fetch date of the img
                dateFile = open(refTxtFile, "r")
                for line in dateFile:
                    if file_ in line:
                        dateStr = line.split("=")[2]
                        print(dateStr[:-1])
                        dst = savingPath + dateStr[:-1] + prefix + "_SL_Subset.tif"

                        src = imgPath + file_
                        os.rename(src, dst)
                        print(i)
                        i = i + 1

        return

    def SubsetSentinelTCI(self, windDim, savingPath, inputDataPath_SL):
        # "gdal_translate -projwin 449230.87415403477 4034940.424743298 472437.7298007704 4015846.1764263636 -of GTiff -co NBITS=16
        gdalTranslateCommand = "gdal_translate -projwin %s %s %s %s -of GTiff -co NBITS=8 " % (
            str(windDim[0]), str(windDim[1]),
            str(windDim[2]), str(windDim[3]))
        ## Create folder to save SL_subsets
        subsetSavingFolder_SL = FileRT.CreateDirectory(directoryPath=savingPath,
                                                       folderName="SL_%s_Subset" % ("TCI"))
        ## Create DateMetadata_Subset File
        dateMetadataSubset_SL = open(subsetSavingFolder_SL + "SubsetDates.txt", "w+")

        folders = FileRT.FilesInDirectory(path=inputDataPath_SL)
        for folder_ in folders:
            pathfolder_ = inputDataPath_SL + folder_ + "/"
            safeFolders = FileRT.FilesInDirectory(path=pathfolder_)
            for safe_ in safeFolders:
                if ".SAFE" in safe_:
                    safeFolderPath = pathfolder_ + safe_
                    granulePath = safeFolderPath + "/" + "GRANULE/"  # Granule Folder PAth
                    temp_ = FileRT.FilesInDirectory(granulePath)  # The folder inside Granule folder

                    xmls = FileRT.FilesInDirectory(path=granulePath + temp_[0])
                    for xml_ in xmls:
                        if ".xml" in xml_:
                            mtdFilePath = granulePath + temp_[0] + "/" + xml_  # MTD_TL.xml file path
                            print(xml_)
                    print(mtdFilePath)

                    imgDataPath = granulePath + temp_[0] + "/IMG_DATA/"  # The IMG_Data folder
                    imgs = FileRT.FilesInDirectory(path=imgDataPath)

                    for img_ in imgs:  # Loop through all imgs inside IMG_DATA folder
                        if "_TCI" in img_:
                            img_path = imgDataPath + img_
                            title = img_[0:-4] + "_Subset.tif"
                            # print(title)
                            command_ = gdalTranslateCommand + img_path + " " + subsetSavingFolder_SL + title
                            # print(command_)
                            os.system(command_)
                            dateStr = self.OpenMTDxml(pathFile=mtdFilePath)
                            dateMetadataSubset_SL.write("SubsetId = " + title + " Date =" + dateStr + "\r\n")
        return

    def CreateRGB_SL(self):

        folders = FileRT.FilesInDirectory(path=self.inputDataPath_SL)
        for folder_ in folders:
            pathfolder_ = self.inputDataPath_SL + folder_ + "/"
            safeFolders = FileRT.FilesInDirectory(path=pathfolder_)
            for safe_ in safeFolders:
                if ".SAFE" in safe_:
                    safeFolderPath = pathfolder_ + safe_
                    granulePath = safeFolderPath + "/" + "GRANULE/"  # Granule Folder PAth
                    temp_ = FileRT.FilesInDirectory(granulePath)  # The folder inside Granule folder

                    mtdFilePath = granulePath + temp_[0] + "/MTD_TL.xml"  # MTD_TL.xml file path
                    imgDataPath = granulePath + temp_[0] + "/IMG_DATA/"  # The IMG_Data folder
                    xmls = FileRT.FilesInDirectory(path=granulePath + temp_[0])
                    for xml_ in xmls:
                        if ".xml" in xml_:
                            mtdFilePath = granulePath + temp_[0] + "/" + xml_  # MTD_TL.xml file path
                            print(xml_)
                    print(mtdFilePath)
                    dateStr = self.OpenMTDxml(pathFile=mtdFilePath)
                    imgs = FileRT.FilesInDirectory(path=imgDataPath)
                    for img_ in imgs:  # Loop through all imgs inside IMG_DATA folder
                        if "_B02.jp2" in img_:
                            blueBandPath = imgDataPath + img_

                        if "_B03.jp2" in img_:
                            greenBandPath = imgDataPath + img_
                        if "_B04.jp2" in img_:
                            redBandPath = imgDataPath + img_

                    vrtOutputFilePath = self.savingPath + "RGB_" + dateStr + "_SL_"
                    commandBuildvrt = "gdalbuildvrt -separate %s %s %s %s" % (
                        vrtOutputFilePath, redBandPath, greenBandPath, blueBandPath)
                    os.system(commandBuildvrt)
                    tifOutputFilePath = vrtOutputFilePath + ".tif"
                    commandTranslateVrt2Tiff = "gdal_translate %s %s" % (vrtOutputFilePath, tifOutputFilePath)
                    os.system(commandTranslateVrt2Tiff)

        return

    def SubsetRGB(self):
        gdalTranslateCommand = "gdal_translate -projwin %s %s %s %s -of GTiff " % (
            str(self.windDim[0]), str(self.windDim[1]),
            str(self.windDim[2]), str(self.windDim[3]))

        folders = FileRT.FilesInDirectory(path=self.inputDataPath_SL)
        for file_ in folders:

            if ".tif" in file_ and ".enp" not in file_:
                pathImg_SL = self.inputDataPath_SL + file_
                title = file_[:-4] + "Subset.tif"
                print(title)
                command_ = gdalTranslateCommand + pathImg_SL + " " + self.savingPath + title
                os.system(command_)

        return

    def ConvertJP2ToGeoTiff(self, inputFolder, outputFolder):
        files = FileRT.FilesInDirectory(path=inputFolder)
        for index, img in enumerate(files):
            imgRasterInfo = RT.GetRasterInfo(inputRaster=inputFolder + img)
            title = img[:-4] + ".tif"
            print("title=", title)
            RT.WriteRaster(refRasterPath=inputFolder + img, newRasterPath=outputFolder + title,
                           Listarrays=[RT.ImageAsArray(imgRasterInfo, 1)], numberOfBands=1)

        return


if __name__ == '__main__':
    sl = Sentinel()
    sl.ConvertJP2ToGeoTiff(inputFolder="/home/cosi/5-Temp_SA/Pamirselected/",
                           outputFolder="/home/cosi/5-Temp_SA/PamirTif/")
