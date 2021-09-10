import sys, os, tarfile
import dateutil.parser as dparser
import numpy as np
from osgeo import gdal, osr
from pathlib import Path

import Main_routine as RT
import FilesCommandRoutine as FileRT
import Correlation_misc as CorrMisc


class Landsat:
    def __init__(self):
        self.inputDataPath_LC = ""
        self.savingPath = ""
        self.windDim = []
        self.band2Subset = "B8"

        self.pathMasterImg = ""
        self.pathSlaveImg = "",
        self.corrOption = ["frequency", "32", "32", "32", "32", "4", "4", "4", "0.9", "1", "1"]
        self.outputSavingPath = ""

        self.corrFolder = ""
        self.datesFile = ""

    def UnzipLandsat(self):
        """
        :param lcRawDataDirectoryPath: Landsat directory path :string
        :param savingPath: where to create the folder to store unzipped files :string
        :return:
        """
        files = FileRT.FilesInDirectory(path=self.inputDataPath_LC)
        for file_ in files:
            if ".tar.gz" in file_:
                filePath_ = self.inputDataPath_LC + file_
                title = file_[0:-13]
                print("Unzipping:", title)
                tar = tarfile.open(filePath_, 'r:gz')
                tar.extractall(path=self.savingPath + title + "/")
                tar.close()
        return

    def SubsetLandsat(self):
        # "gdal_translate -projwin 449230.87415403477 4034940.424743298 472437.7298007704 4015846.1764263636 -of GTiff -co NBITS=16
        gdalTranslateCommand = "gdal_translate -projwin %s %s %s %s -of GTiff -co NBITS=16 " % (
            str(self.windDim[0]), str(self.windDim[1]),
            str(self.windDim[2]), str(self.windDim[3]))
        ## Create folder to save SL_subsets
        subsetSavingFolder_LC = FileRT.CreateDirectory(directoryPath=self.savingPath,
                                                       folderName="LC_%s_Subset" % (self.band2Subset))
        ## Create DateMetadata_Subset File
        dateMetadataSubset_LC = open(subsetSavingFolder_LC + "SubsetDates.txt", "w+")

        folders = FileRT.FilesInDirectory(path=self.inputDataPath_LC)
        for folder_ in folders:
            print(folder_)
            pathfolder_ = self.inputDataPath_LC + folder_ + "/"
            imgs = FileRT.FilesInDirectory(path=pathfolder_)
            for img_ in imgs:  # Loop through all imgs inside IMG_DATA folder
                if self.band2Subset in img_:
                    img_path = pathfolder_ + img_
                    title = img_[0:-4] + "_Subset.tif"
                    print(title)
                    command_ = gdalTranslateCommand + img_path + " " + subsetSavingFolder_LC + title
                    print(command_)
                    os.system(command_)

                # if "MTL.txt" in img_:
                #     mtlPath = pathfolder_ +img_
                #     dateStr = self.OpenMtlFile(mtlPath=mtlPath)
                #     dateMetadataSubset_LC.write("SubsetId = " + title + " Date =" + dateStr + "\r\n")
        return

    def OpenMtlFile(self, mtlPath=""):
        """
        Function to read MTL file associated with landsat images and fetch the acquisition date
        :param mtlPath: MTL file path : string
        :return: the acquisition date  : string
        """
        with open(mtlPath, "rt") as mtlFile:
            for line in mtlFile:  # For each line, read it to a string
                if "DATE_ACQUIRED" in line:
                    print(line[0:-1])
                    return (line.split('=')[1][0:-1])

    def CreateBatchCorrelationFiles(self):
        """
        :param pathMasterImg:
        :param pathSlaveImg:
        :param corrOption:
        :return:
        """

        CorrMisc.CreateBatchCorrelationFiles(pathMasterImgs=self.pathMasterImg, pathSlaveImgs=self.pathSlaveImg,
                                             savingPath=self.savingPath, corrOptions=self.corrOption)

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

            titleEW = "Master_" + dateMasterString + "_vs_Slave_%s_LC_149_%s_EW.tif" % (
                dateSlaveString, masterImg.split("_")[2])
            titlePathEW = pathEW + titleEW

            RT.Array2MultiBandRaster(refRasterPath=self.corrFolder + corrFile_, newRasterPath=titlePathEW,
                                     arrays=[ewArray, snrArray], numberOfBands=2,
                                     descriptions=[EW.GetDescription(), SNR.GetDescription()])

            titleNS = "Master_" + dateMasterString + "_vs_Slave_%s_LC_149_%s_NS.tif" % (
                dateSlaveString, masterImg.split("_")[2])
            titlePathNS = pathNS + titleNS

            RT.Array2MultiBandRaster(refRasterPath=self.corrFolder + corrFile_, newRasterPath=titlePathNS,
                                     arrays=[nsArray, snrArray], numberOfBands=2,
                                     descriptions=[NS.GetDescription(), SNR.GetDescription()])

        return

    def ComputeEWNSVelocity(self, corrFolder, outputSavingPath, withSNR=True):
        """
        : coor file should be defined
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
            with open(self.datesFile, "r") as dateFile:
                for line in dateFile:
                    if masterImg in line:
                        dateMasterString = line.split("=")[-1][1:-1]
                        dateMaster = dparser.parse(line.split("=")[-1], fuzzy=True)
                        print("Master=", dateMaster)
                    if slaveImg in line:
                        dateSlaveString = line.split("=")[-1][1:-1]
                        dateSlave = dparser.parse(line.split("=")[-1], fuzzy=True)
                        print("Slave=", dateSlave)
            ## Compute the delay between the tow images
            daysDiff = abs(dateMaster - dateSlave)
            print("daysDiff=", daysDiff.days)

            ## Create a velocity map
            title = "Master_" + dateMasterString + "_vs_Slave_%s_LC_%s_velocity.tif" % (
                dateSlaveString, masterImg.split("_")[2])
            print("title=", title)
            titlePath = self.savingPath + title

            corr = gdal.Open(self.corrFolder + corrFile_)
            EW = corr.GetRasterBand(1)
            NS = corr.GetRasterBand(2)
            SNR = corr.GetRasterBand(3)
            ewArray = EW.ReadAsArray()
            nsArray = NS.ReadAsArray()
            snrArray = SNR.ReadAsArray()
            velocityArray1 = np.sqrt(ewArray ** 2 + nsArray ** 2) / daysDiff.days
            if boolSNRband == False:
                RT.Array2MultiBandRaster(refRasterPath=self.corrFolder + corrFile_, newRasterPath=titlePath,
                                         arrays=[velocityArray1], numberOfBands=1,
                                         descriptions=["Velocity"])
            if boolSNRband == True:
                RT.Array2MultiBandRaster(refRasterPath=self.corrFolder + corrFile_, newRasterPath=titlePath,
                                         arrays=[velocityArray1, snrArray], numberOfBands=2,
                                         descriptions=["Velocity", SNR.GetDescription()])
        return

    def ApplyROIMask(self, inputDirectory, mask_Path, outputDirectroy="",setNaN=True):
        """
        Apply mask to map
        :param InputDirectory: images path : string
        :param outputDirectroy: results output path : string
        :param mask_Path: path of the mask file : string
        :return:
        """
        if not outputDirectroy:
            workspacePath = Path(inputDirectory).parents[1]
            folders = FileRT.FilesInDirectory(path=workspacePath)
            if "MaskedData" in folders:
                tempPath = os.path.join(workspacePath, "MaskedData")
            else:
                tempPath = FileRT.CreateDirectory(directoryPath=workspacePath, folderName="MaskedData")

            outputDirectroy = FileRT.CreateDirectory(directoryPath=tempPath,
                                                     folderName=os.path.basename(
                                                         os.path.normpath(inputDirectory)) + "_Masked")

        maskRaster = gdal.Open(mask_Path)
        maskArray = maskRaster.GetRasterBand(1).ReadAsArray()
        maskArray = maskArray.astype(float)
        if setNaN == True:
            maskArray[maskArray == 0] = np.nan

        FileRT.DirectoryAsEmpty(outputDirectroy)
        files = FileRT.FilesInDirectory(inputDirectory)

        for file_ in files:
            if ".enp" not in file_ and ".hdr" not in file_:
                print("Masking:", file_)
                arrays = []
                descriptions = []
                rasterInfo = RT.GetRasterInfo(inputDirectory + file_)
                raster = rasterInfo["Raster"]
                for i in range(raster.RasterCount):
                    bandNumber = i + 1
                    bandTmp = raster.GetRasterBand(bandNumber)
                    bandTmpArr = RT.ImageAsArray(rasterInfo, bandNumber=bandNumber)
                    arrays.append(bandTmpArr * maskArray)
                    descriptions.append(bandTmp.GetDescription())

                title = outputDirectroy + file_[:-4] + "_Masked.tif"
                metaData = rasterInfo["MetaData"]
                metaData["MaskPath"] = mask_Path
                RT.WriteRaster(refRasterPath=inputDirectory + file_, newRasterPath=title,
                               Listarrays=arrays, numberOfBands=raster.RasterCount, descriptions=descriptions,
                               metaData=metaData)

        return

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

    def RenameWithDate(self, refTxtFile, imgPath, savingPath, sufix):
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
                        print(dateStr[1:-1])
                        dst = savingPath + dateStr[1:-1] + sufix

                        src = imgPath + file_
                        os.rename(src, dst)
                        print(i)
                        i = i + 1

        return
