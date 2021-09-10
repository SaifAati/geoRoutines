import datetime
import os
from pathlib import Path
import dateutil.parser as dparser
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import geospatialroutine.FilesCommandRoutine as FileRT
import geospatialroutine.Routine as RT


###############################################  Main functions (.h) ###################################################
def CorrelationStrategy(baseImgList, targetImgList, baseFilePath, targetFilePath, corrOptions, corrOptionFile,
                        outputFile, savingPath, corrStrat=1):
    corrStrategy = corrStrat

    if corrStrategy == 1:
        baseImgList = baseImgList[:-1]
        targetImgList = targetImgList[1:]
        ## Create base File

        with open(baseFilePath, 'w') as f:
            for item in baseImgList:
                item_ = item
                f.write("%s\n" % item_)

        ## Create target File

        with open(targetFilePath, 'w') as f:
            for item in targetImgList:
                item_ = item
                f.write("%s\n" % item_)

        ## Create Correlation option file

        corrCommand = ""
        for i in range(len(corrOptions)):
            if i == 0:
                corrCommand = corrCommand + corrOptions[i]
            else:
                corrCommand = corrCommand + " " + corrOptions[i]
        print(corrCommand)

        with open(corrOptionFile, 'w') as f:
            for i in range(len(baseImgList)):
                f.write("%s\n" % corrCommand)

        ## Create Output File

        with open(outputFile, 'w') as f:
            for i in range(len(baseImgList)):
                item = os.path.join(savingPath,
                                    "Corr_" + Path(baseImgList[i]).stem + "_vs_" + Path(
                                        targetImgList[i]).stem + "_" + str(i))
                f.write("%s\n" % item)
    if corrStrategy == 2:
        newBaseList = []
        newTargetList = []
        outputList = []
        print(baseImgList)
        print("Correlation Strategy all")
        nbCorrelation = int(len(baseImgList) * (len(baseImgList) - 1) / 2)
        print("Number of correlation:", nbCorrelation)
        for i in range(len(baseImgList) - 1):
            newMList_ = (len(baseImgList) - 1 - i) * [baseImgList[i]]
            newBaseList.extend(newMList_)
            for j in range(i + 1, len(targetImgList)):
                newTargetList.append(targetImgList[j])

        print(newBaseList, len(newBaseList))
        print(newTargetList, len(newTargetList))

        with open(baseFilePath, 'w') as f:
            for item in newBaseList:
                item_ = item
                f.write("%s\n" % item_)

        ## Create target File

        with open(targetFilePath, 'w') as f:
            for item in newTargetList:
                item_ = item
                f.write("%s\n" % item_)

        ## Create Output File

        with open(outputFile, 'w') as f:
            for i in range(len(newBaseList)):
                item = os.path.join(savingPath,
                                    Path(newBaseList[i]).stem + "_VS_" + Path(newTargetList[i]).stem)
                f.write("%s\n" % item)

        ## Create Correlation option file

        corrCommand = ""
        for i in range(len(corrOptions)):
            if i == 0:
                corrCommand = corrCommand + corrOptions[i]
            else:
                corrCommand = corrCommand + " " + corrOptions[i]
        print(corrCommand)

        with open(corrOptionFile, 'w') as f:
            for i in range(len(newBaseList)):
                f.write("%s\n" % corrCommand)
    if corrStrategy == 3:
        newBaseList = []
        newTargetList = []
        outputList = []
        print(baseImgList)
        print("Correlation Strategy all (Back and Forward)")
        nbCorrelation = int(len(baseImgList) * (len(baseImgList) - 1))
        print("Number of correlation:", nbCorrelation)
        for base_ in baseImgList:
            newTargetList_ = []
            for target_ in targetImgList:
                if base_ != target_:
                    newBaseList.append(base_)
                    newTargetList_.append(target_)
            newTargetList.extend(newTargetList_)

        with open(baseFilePath, 'w') as f:
            for item in newBaseList:
                item_ = item
                f.write("%s\n" % item_)

        with open(targetFilePath, 'w') as f:
            for item in newTargetList:
                item_ = item
                f.write("%s\n" % item_)

        # ## Create Output File

        with open(outputFile, 'w') as f:
            for i in range(len(newBaseList)):
                item = os.path.join(savingPath,
                                    Path(newBaseList[i]).stem + "_VS_" + Path(newTargetList[i]).stem)
                f.write("%s\n" % item)
        #
        # ## Create Correlation option file

        corrCommand = ""
        for i in range(len(corrOptions)):
            if i == 0:
                corrCommand = corrCommand + corrOptions[i]
            else:
                corrCommand = corrCommand + " " + corrOptions[i]
        # print(corrCommand)
        #
        with open(corrOptionFile, 'w') as f:
            for i in range(len(newBaseList)):
                f.write("%s\n" % corrCommand)
    if corrStrategy == 4:

        from itertools import permutations
        newBaseList = []
        newTargetList = []
        outputList = []
        print(baseImgList)
        print("#############################################")
        print("--Correlation Strategy all (Back and Forward)")
        print("--Base list is the same as the target list---")
        print("#############################################")
        # nbCorrelation = int(len(baseImgList) * (len(baseImgList) - 1))

        perm = permutations(baseImgList,2)
        permList = []
        for i in list(perm):
            permList.append(i)
            newBaseList.append(i[0])
            newTargetList.append(i[1])
            outputList.append(os.path.join(savingPath,Path(i[0]).stem+"_VS_"+Path(i[1]).stem))
        # print(permList,len(permList))
        print("Number of correlation:", len(permList))

        with open(baseFilePath, 'w') as f:
            for item in newBaseList:
                item_ = item
                f.write("%s\n" % item_)

        with open(targetFilePath, 'w') as f:
            for item in newTargetList:
                item_ = item
                f.write("%s\n" % item_)
        with open(outputFile, 'w') as f:
            for item in outputList:
                item_ = item
                f.write("%s\n" % item_)
        corrCommand = ""
        for i in range(len(corrOptions)):
             if i == 0:
                 corrCommand = corrCommand + corrOptions[i]
             else:
                 corrCommand = corrCommand + " " + corrOptions[i]

        with open(corrOptionFile, 'w') as f:
            for i in range(len(newBaseList)):
                f.write("%s\n" % corrCommand)


    return


########################################################################################################################


################################ Warping Functions #####################################################################
def CreateBatchCorrelationFiles(pathbaseImgs, pathtargetImgs, savingPath, corrOptions, corrStrategy=1):
    """

    :param pathbaseImgs:
    :param pathtargetImgs:
    :param savingPath:
    :param corrOptions: ["frequency","64","64","32","32","8","8","4","0.9","1","1"]
    :param corrStrategy: 1 sequence correlation: e.g for N image we have N-1 correlation
                         2 all correlation :e.g for N image we have N(N-1)/2 correlation
                         3 all correlation back and forward :e.g for N image we have N(N-1) correlation

    :return:
    """

    filterList = [".hdr", ".geojson", ".anc", ".txt", ".svg", ".csv", ".qgz"]
    # list1 = FileRT.FilesInDirectory(path=pathbaseImgs)
    # baseImgList = [os.path.join(pathbaseImgs, item) for item in list1 if
    #                any(i in item for i in filterList) == False]
    # list1 = FileRT.FilesInDirectory(path=pathtargetImgs)
    # targetImgList = [os.path.join(pathtargetImgs, item) for item in list1 if
    #                  any(i in item for i in filterList) == False]
    baseFilePath = os.path.join(savingPath, "0-Base.txt")
    targetFilePath = os.path.join(savingPath, "1-Target.txt")
    corrOptionFile = os.path.join(savingPath, "2-CorrOptionFile.txt")
    outputFile = os.path.join(savingPath, "3-OutputFile.txt")

    baseImgList = [
        "/media/cosicorr/storage/Saif/1-Ridgecrest/1-WV/P1BS_sorted/Before_eq/Before_eq/2018-06-16-WV1/502260874070_01_P002_PAN/Ortho_150cm/18JUN16213928-P1BS-502260874070_01_P002_Ortho_150cm_LidarDEM",
        "/media/cosicorr/storage/Saif/1-Ridgecrest/1-WV/P1BS_sorted/After_eq/Post_eq/2019-07-14_wv2_stereo/WV220190714183850P00/19JUL14183850-P1BS-504254136040_01_P004/Ortho_150cm/19JUL14183850-P1BS-504254136040_01_P004_Ortho_150cm_LidarDEM",
        "/media/cosicorr/storage/Saif/1-Ridgecrest/1-Spot6_WorkSpace/Pre_SPOT6_P_20180915/Ortho_150cm_LidarDEM/IMG_SPOT6_P_201809151819247_SEN_4555226101_R1C1_Ortho_150cm_Lidar__.tif",
        "/media/cosicorr/storage/Saif/1-Ridgecrest/1-Spot6_WorkSpace/Post_Spot_RPC_SPOT6_20190724/Ortho_150cm_Lidar/IMG_SPOT6_P_201907241820211_SEN_4555225101_R1C1_Ortho_150cm_LidarDem_subset_.tif"]
    targetImgList = baseFilePath
    CorrelationStrategy(baseImgList=baseImgList,
                        targetImgList=targetImgList,
                        baseFilePath=baseFilePath,
                        targetFilePath=targetFilePath,
                        corrOptions=corrOptions,
                        corrOptionFile=corrOptionFile,
                        outputFile=outputFile,
                        savingPath=savingPath,
                        corrStrat=corrStrategy)

    return


########################################################################################################################


########################################################################################################################


def GetDiffDays(refDate, listDates):
    """
    Transform a list of dates in string format, to a list of differences in days compared to the reference date
    :param refDate: referencee date: string
    :param listDate: list of dates : List[string]
    return : list of dates in days
    """
    refDate_ = dparser.parse(refDate, fuzzy=True)
    listDates_ = []
    for dateStr in listDates:
        daysDiff = dparser.parse(dateStr, fuzzy=True) - refDate_
        listDates_.append(daysDiff.days)
    return listDates_


def GetDatesFromDiffDays(refDate, diffDaysList):
    refDate_ = dparser.parse(refDate, fuzzy=True)
    listDates = []
    for dateDays in diffDaysList:
        date = refDate_ + datetime.timedelta(days=int(dateDays))
        listDates.append(date.strftime("%Y-%m-%d"))
    return listDates


def GetMidDates(datebase, datetarget):
    datebase = datetime.datetime.strptime(datebase, '%Y-%m-%d')
    # print("datebase=", datebase)
    datetarget = datetime.datetime.strptime(datetarget, '%Y-%m-%d')
    # print("datetarget=", datetarget, )
    deltaT = int((datetarget - datebase).days)
    # print(deltaT,"\n")
    midDate = datebase + (datetarget - datebase) / 2
    # print(midDate)

    return midDate, deltaT


def VelocityFileInfo(velocityFileName, option=1):
    """
    :param corrFile: name of the correlationn file raster: string
    :param option: 1 the name contain the base and salve date
            option 2: the name contain the target date
            option 3: form the metadata of the raster
            option 4: only one date (new Time step)
    :return: doctrinaire of correlation information based on the name of the correlation raster: dic
    """
    if option == 1:
        vInfo = {"base": {}, "target": {}}
        infoList = os.path.basename(velocityFileName).split("_vs_")
        infobase = infoList[0]
        infotarget = infoList[1]
        infobase = infobase.split("_")
        infotarget = infotarget.split("_")
        vInfo["base"]["Date"] = infobase[1]
        vInfo["base"]["Satellite"] = infobase[2]
        vInfo["base"]["Band"] = infobase[3]

        vInfo["target"]["Date"] = infotarget[1]
        vInfo["target"]["Satellite"] = infotarget[2]
        vInfo["target"]["Band"] = infotarget[3]

        vInfo["Component"] = infotarget[4]

        datebase = dparser.parse(vInfo["base"]["Date"], fuzzy=True)
        # print("base Date=", datebase)
        datetarget = dparser.parse(vInfo["target"]["Date"], fuzzy=True)
        # print("target Date=", datetarget)

        ## Compute the delay between the tow images
        daysDiff = abs(datebase - datetarget)
        vInfo["TimeSpan"] = daysDiff.days
        vInfo["MidDate"], _ = GetMidDates(datebase, datebase)
        # print("daysDiff=", daysDiff.days)
        print("velocityInfoFile:", vInfo)
        return vInfo
    if option == 2:
        vInfo = {"target": {}}
        infoList = os.path.basename(velocityFileName).split("_")
        vInfo["target"]["Date"] = infoList[1]
        vInfo["target"]["Satellite"] = infoList[2]
        vInfo["target"]["Band"] = infoList[3]
        vInfo["Component"] = infoList[4]
        return vInfo
    if option == 3:
        rasterInfo = RT.GetRasterInfo(inputRaster=velocityFileName)
        metaData = rasterInfo["Raster"].GetMetadata()
        print(metaData)
        vInfo = {"base": {}, "target": {}}
        vInfo["base"]["Date"] = metaData["MasterDate"]
        vInfo["base"]["Satellite"] = metaData["MasterPlatform"]
        vInfo["base"]["Band"] = metaData["MasterBand"]

        vInfo["target"]["Date"] = metaData["SlaveDate"]
        vInfo["target"]["Satellite"] = metaData["SlavePlatform"]
        vInfo["target"]["Band"] = metaData["SlaveBand"]

        vInfo["MidDate"], vInfo["TimeSpan"] = GetMidDates(vInfo["base"]["Date"], vInfo["target"]["Date"])

        print("velocityInfoFile:", vInfo)
        return vInfo
    if option == 4:
        ## Return all the meta data
        rasterInfo = RT.GetRasterInfo(inputRaster=velocityFileName)
        metaData = rasterInfo["Raster"].GetMetadata()
        return metaData


def CorrFileInfo(corrFileName, oneDic=False):
    """
    :param corrFile: name of the correlationn file raster: string
    :return: doctrinaire of correlation information based on the name of the correlation raster: dic
    """

    if not oneDic:
        corrInfo = {"base": {}, "target": {}}
        infoList = os.path.basename(corrFileName).split("_")
        corrInfo["base"]["Date"] = infoList[1]
        corrInfo["base"]["Satellite"] = infoList[2]
        corrInfo["base"]["Band"] = infoList[3]

        infoList = os.path.basename(corrFileName).split("_vs_")[1].split("_")
        corrInfo["target"]["Date"] = infoList[0]
        corrInfo["target"]["Satellite"] = infoList[1]
        corrInfo["target"]["Band"] = infoList[2]

        datebase = dparser.parse(corrInfo["base"]["Date"], fuzzy=True)
        # print("base Date=", datebase)
        datetarget = dparser.parse(corrInfo["target"]["Date"], fuzzy=True)
        # print("target Date=", datetarget)

        ## Compute the delay between the tow images
        daysDiff = abs(datebase - datetarget)
        corrInfo["TimeSpan"] = daysDiff.days
        # print("daysDiff=", daysDiff.days)
        print("corrInfoFile:", corrInfo)
    if oneDic:
        corrInfo = {}
        infoList = os.path.basename(corrFileName).split("_")
        corrInfo["baseDate"] = infoList[1]
        corrInfo["baseSatellite"] = infoList[2]
        corrInfo["baseBand"] = infoList[3]

        infoList = os.path.basename(corrFileName).split("_vs_")[1].split("_")
        corrInfo["targetDate"] = infoList[0]
        corrInfo["targetSatellite"] = infoList[1]
        corrInfo["targetBand"] = infoList[2]

        datebase = dparser.parse(corrInfo["baseDate"], fuzzy=True)
        # print("base Date=", datebase)
        datetarget = dparser.parse(corrInfo["targetDate"], fuzzy=True)
        # print("target Date=", datetarget)

        ## Compute the delay between the tow images
        daysDiff = abs(datebase - datetarget)
        corrInfo["TimeSpan"] = daysDiff.days
        # print("daysDiff=", daysDiff.days)
        print("corrInfoFile:", corrInfo)

    return corrInfo


def ComputeEWNSVelocity(corrFilePath):
    """
    :param corrFilePath: path of the correlation raster, we assume the correlation raster contain 3 maps EW, NS and SNR
    :return:
    """
    velocities = {"EW": {}, "NS": {}, "SNR": {}}
    rasterInfo = RT.GetRasterInfo(inputRaster=corrFilePath)
    if rasterInfo["Error"]:
        raise Exception("Could not open raster with Gdal:" + corrFilePath)

    print("Correlation raster Info:", rasterInfo)
    if not rasterInfo["NbBands"] == 3:
        raise Exception('Could not find 3 bands in the correlation raster  ' + rasterInfo["ImageName"])

    ewArray = RT.ImageAsArray(imageInfo=rasterInfo, bandNumber=1)
    nsArray = RT.ImageAsArray(imageInfo=rasterInfo, bandNumber=2)
    snrArray = RT.ImageAsArray(imageInfo=rasterInfo, bandNumber=3)

    corrFileInfo_ = CorrFileInfo(corrFileName=rasterInfo["ImageName"])

    prefix = "base_" + corrFileInfo_["base"]["Date"] + "_" + corrFileInfo_["base"]["Satellite"] + "_" + \
             corrFileInfo_["base"]["Band"] + "_vs_target_" + corrFileInfo_["target"]["Date"] + "_" + \
             corrFileInfo_["target"]["Satellite"] + "_" + corrFileInfo_["target"]["Band"]

    velocities["baseDate"] = corrFileInfo_["base"]["Date"]
    velocities["targetDate"] = corrFileInfo_["target"]["Date"]
    velocities["baseBand"] = corrFileInfo_["base"]["Band"]
    velocities["targetBand"] = corrFileInfo_["target"]["Band"]
    velocities["basePlatform"] = corrFileInfo_["base"]["Satellite"]
    velocities["targetPlatform"] = corrFileInfo_["target"]["Satellite"]
    velocities["TimeSpan"] = corrFileInfo_["TimeSpan"]

    velocities["EW"]["FileName"] = prefix + "_EW_Velocity.tif"
    velocities["EW"]["Array"] = ewArray / corrFileInfo_["TimeSpan"]
    velocities["EW"]["Description"] = "E/W:" + prefix + "_EW_Velocity.tif"

    velocities["NS"]["FileName"] = prefix + "_NS_Velocity.tif"
    velocities["NS"]["Array"] = nsArray / corrFileInfo_["TimeSpan"]
    velocities["NS"]["Description"] = "N/S:" + velocities["NS"]["FileName"]

    velocities["SNR"]["Array"] = snrArray
    velocities["SNR"]["Description"] = "SNR"

    print("EW-NS velocities maps:", velocities)
    return velocities


def ComputeVelocityMap(ewVelocityPath, nsVelocityPath, savingPath):
    """

    """
    ## Check if we have the same number of images in both folders
    ewImgs = FileRT.FilesInDirectory(path=ewVelocityPath)
    nsImgs = FileRT.FilesInDirectory(path=nsVelocityPath)
    if len(ewImgs) != len(nsImgs):
        print("Error, EW and NS folder don't have the same number of images ")
        return
    for index, ewImg_ in enumerate(ewImgs):
        ewImgPath = ewVelocityPath + ewImg_
        nsImgPath = nsVelocityPath + nsImgs[index]

        ewRasterInfo = RT.GetRasterInfo(inputRaster=ewImgPath)
        nsRasterInfo = RT.GetRasterInfo(inputRaster=nsImgPath)

        ewArray = RT.ImageAsArray(imageInfo=ewRasterInfo)
        nsArray = RT.ImageAsArray(imageInfo=nsRasterInfo)
        velocityArray = np.sqrt(ewArray ** 2 + nsArray ** 2)
        name = ewImg_.split(("EW"))[0] + "velocity.tif"
        metaData = ewRasterInfo["MetaData"]
        metaData["ewVelocityMap"] = ewVelocityPath
        metaData["nsVelocityMap"] = nsVelocityPath
        RT.WriteRaster(refRasterPath=ewImgPath, newRasterPath=savingPath + name, Listarrays=[velocityArray],
                       numberOfBands=1, descriptions=["VelocityMap computed based on EW and NS  "], metaData=metaData)

    return


def FromtargetDateGetbaseDate(targetDate, velocityFolder):
    ## Open the folder of data to be displayed
    files = FileRT.FilesInDirectory(velocityFolder)
    files.sort()
    fileList = []  # list of raster to be displayed
    dateList = []  # list of date
    satType = []

    for file_ in files:
        if ".tif" in file_ and ".xml" not in file_ and ".enp" not in file_:
            if "target_" + targetDate in file_:
                vInfo = VelocityFileInfo(velocityFileName=file_, option=1)
                return vInfo["base"]["Date"]


def CorrelationReport(ewPath, nsPath, outputFolder, reportFileName, saveReport=True, decimal=4, latex=False):
    """

    :param ewPath:
    :param nsPath:
    :param outputFolder:
    :param saveReport:
    :param decimal:
    :return:
    """
    from pandas import DataFrame
    ewImgs = FileRT.FilesInDirectory(path=ewPath)
    nsImgs = FileRT.FilesInDirectory(path=nsPath)
    data = {"Pre-Ortho [date/platform]": [], "Post-Ortho [date/platform]": [], "Time span [d]": [],
            "EW mean/std/min/max [m]": [], "NS mean/std/min/max [m]": [], "SNR mean/std [m]": []}

    for index, img_ in enumerate(ewImgs):
        ewRasterPath = os.path.join(ewPath, img_)
        ewRasterInfo = RT.GetRasterInfo(inputRaster=ewRasterPath)
        data["Pre-Ortho [date/platform]"].append(
            (ewRasterInfo["MetaData"]["baseDate"] + " (" + ewRasterInfo["MetaData"]["basePlatform"] + ")"))
        data["Post-Ortho [date/platform]"].append(
            (ewRasterInfo["MetaData"]["targetDate"] + " (" + ewRasterInfo["MetaData"]["targetPlatform"] + ")"))
        data["Time span [d]"].append(ewRasterInfo["MetaData"]["TimeSpan"])

        # compute stats
        ewMean, ewStd, ewMax, ewMin = CorrStat(rasterPath=ewRasterPath, bandNb=1)
        data["EW mean/std/min/max [m]"].append(
            str(np.round(ewMean, decimal)) + "/" + str(np.round(ewStd, decimal)) + "/" + str(
                np.round(ewMin, decimal)) + "/" + str(
                np.round(ewMax, decimal)))

        snrMean, snrStd, snrMax, snrMin = CorrStat(rasterPath=ewRasterPath, bandNb=2)
        data["SNR mean/std [m]"].append(
            str(np.round(snrMean, decimal)) + "/" + str(np.round(snrStd, decimal)))

        nsRasterPath = os.path.join(nsPath, nsImgs[index])
        nsMean, nsStd, nsMax, nsMin = CorrStat(rasterPath=nsRasterPath, bandNb=1)
        data["NS mean/std/min/max [m]"].append(
            str(np.round(nsMean, decimal)) + "/" + str(np.round(nsStd, decimal)) + "/" + str(
                np.round(nsMin, decimal)) + "/" + str(
                np.round(nsMax, decimal)))

    df = DataFrame(data, columns=["Pre-Ortho [date/platform]", "Post-Ortho [date/platform]", "Time span [d]",
                                  "EW mean/std/min/max [m]", "NS mean/std/min/max [m]", "SNR mean/std [m]"])
    print(df)

    if saveReport:
        base_filename = reportFileName
        with open(os.path.join(outputPath, base_filename), 'w') as outfile:
            df.to_string(outfile, index=False)

    return


def CorrStat(rasterPath, bandNb=1):
    """
    :param rasterPath:
    :param bandNb:
    :return:
    """
    rasterInfo = RT.GetRasterInfo(inputRaster=rasterPath)
    rasterArray = RT.ImageAsArray(rasterInfo, bandNb)

    mean = np.nanmean(rasterArray)
    std = np.nanstd(rasterArray)
    max = np.nanmax(rasterArray)
    min = np.nanmin(rasterArray)
    return mean, std, max, min


def VelocityProfile(inputPath):
    vImgs = FileRT.FilesInDirectory(path=velocityPath)
    meanList = []
    stdList = []
    maxList = []
    dateList = []
    for index, img_ in enumerate(vImgs):
        mean, std, max, min = CorrStat(rasterPath=inputPath + img_)
        rasterInfo = RT.GetRasterInfo(inputRaster=inputPath + img_)
        date_str = rasterInfo["MetaData"]["targetDate"]  #
        # print(date_str)
        date = datetime.datetime.strptime(date_str, '%Y-%m-%d').date()
        dateList.append(date)
        meanList.append(mean)
        maxList.append(max)

    start = datetime.datetime.strptime("2013-01-01", "%Y-%m-%d")
    end = dateList[-1]  # datetime.datetime.strptime("07-07-2014", "%d-%m-%Y")
    # date_generated = [start + datetime.timedelta(days=x * 366) for x in range(1, 7)]
    date_generated = [datetime.datetime.strptime("2014-01-01", "%Y-%m-%d"),
                      datetime.datetime.strptime("2015-01-01", "%Y-%m-%d"),
                      datetime.datetime.strptime("2016-01-01", "%Y-%m-%d"),
                      datetime.datetime.strptime("2017-01-01", "%Y-%m-%d"),
                      datetime.datetime.strptime("2018-01-01", "%Y-%m-%d"),
                      datetime.datetime.strptime("2019-01-01", "%Y-%m-%d")]
    dates = ["2013", "2014", "2015", "2016", "2017", "2018", "2019"]
    print(date_generated)
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 1)
    ax1 = plt.subplot(gs[0, 0])  # row 0, col 0
    # plt.imshow(midValMatrix,origin="lower")
    for index, year in enumerate(date_generated):
        plt.axvline(x=year, color="green")
        ax1.text(year - datetime.timedelta(days=180), 1.65, dates[index])
    ax1.text(year + datetime.timedelta(days=180), 1.65, dates[-1])
    color = 'tab:red'

    ax1.plot(dateList, meanList, color=color, marker="o", label="vMean")
    ax1.tick_params(axis='y', labelcolor=color)
    # im = plt.imshow(midValMatrix)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color2 = 'tab:blue'
    ax2.plot(dateList, maxList, color=color2, marker="o", linestyle="dashed", label="vMax")
    ax2.tick_params(axis='y', labelcolor=color2)

    ax1.set_title("Mean Max velocity of Shisper Galcier ")

    ax1.set_ylabel("Mean velocity [m/day]", color=color)
    ax2.set_ylabel('Max velocity [m/day]', color=color2)  # we already handled the x-label with ax1

    ## X-axis
    ax1.set_xlabel("Date")
    from matplotlib.dates import DateFormatter
    date_form = DateFormatter("%y/%m")
    ax1.xaxis.set_major_formatter(date_form)
    # fig.autofmt_xdate()
    # # plt.yticks(yData)
    xStart, xEnd = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(xStart, xEnd, 180))  ## 180; equivalent to 6 months

    ax1.legend(loc=(0.01, 0.80), frameon=False)
    ax2.legend(loc=(0.01, 0.77), frameon=False)
    plt.show()

    return


def RenameCorr_SalveDate(inputFolder, outputFolder="", sameFolder=True):
    files = FileRT.FilesInDirectory(path=inputFolder)
    files.sort()
    for file_ in files:
        # vInfo = VelocityFileInfo(velocityFileName=file_)
        # print(vInfo)
        # rasterInfo = RT.GetRasterInfo(inputRaster=os.path.join(inputFolder,file_))
        # print(rasterInfo["MetaData"]["targetDate"])
        # newName = vInfo["target"]["Date"]
        targetName = file_.split("_vs_")[1]
        # targetName = "target_"+rasterInfo["MetaData"]["targetDate"]+"_"+rasterInfo["MetaData"]["targetPlatform"]+".tif"
        print(targetName)
        src = os.path.join(inputFolder, file_)
        if sameFolder:
            dst = os.path.join(inputFolder, targetName)
        else:
            dst = os.path.join(outputFolder, targetName)
        os.rename(src, dst)

    return


def ComputeMean_Orientation(pathEW, pathNS, outputPath):
    def angle(ew, ns):
        if ew > 0 and ns > 0:
            return np.pi / 2 - np.abs(np.arctan(ns / ew))
        elif ew > 0 and ns < 0:
            return np.pi / 2 + np.abs(np.arctan(ns / ew))
        elif ew < 0 and ns < 0:
            return (3 * np.pi / 2) - np.abs(np.arctan(ns / ew))
        elif ew < 0 and ns > 0:
            return (3 * np.pi / 2) + np.abs(np.arctan(ns / ew))
        # elif ew==0 and ns!=0:
        #     return np.nan

    ewFiles = FileRT.GetFilesBasedOnExtension(pathEW)
    nsFiles = FileRT.GetFilesBasedOnExtension(pathNS)
    print(len(ewFiles), len(nsFiles))
    imgTemp = RT.ImageAsArray(RT.GetRasterInfo(ewFiles[0]))

    ewDataCube = np.empty((len(ewFiles), imgTemp.size))
    nsDataCube = np.empty((len(nsFiles), imgTemp.size))
    for i in range(len(ewFiles)):
        imgTempEW = RT.ImageAsArray(RT.GetRasterInfo(ewFiles[i])).flatten()
        imgTempNS = RT.ImageAsArray(RT.GetRasterInfo(nsFiles[i])).flatten()
        ewDataCube[i, :] = imgTempEW
        nsDataCube[i, :] = imgTempNS
    ewVcum = np.nansum(ewDataCube, axis=0) / len(ewFiles)
    nsVcum = np.nansum(nsDataCube, axis=0) / len(nsFiles)
    RT.WriteRaster(refRasterPath=ewFiles[0], newRasterPath=os.path.join(outputPath, "EW_Cum_norm.tif"),
                   Listarrays=[ewVcum.reshape(imgTemp.shape)], numberOfBands=1)
    RT.WriteRaster(refRasterPath=ewFiles[0], newRasterPath=os.path.join(outputPath, "NS_Cum_norm.tif"),
                   Listarrays=[nsVcum.reshape(imgTemp.shape)], numberOfBands=1)
    anglesFlatten = np.copy(nsVcum)
    for i in range(len(ewVcum)):
        anglesFlatten[i] = angle(ns=nsVcum[i], ew=ewVcum[i])

    anglesRad = anglesFlatten.reshape(imgTemp.shape)

    RT.WriteRaster(refRasterPath=ewFiles[0], newRasterPath=os.path.join(outputPath, "OrientationMap.tif"),
                   Listarrays=[anglesRad], numberOfBands=1)
    anglesRaster = RT.ImageAsArray(RT.GetRasterInfo(os.path.join(outputPath, "OrientationMap.tif")))

    print(anglesRaster)
    angles = np.degrees(anglesRaster)
    print(np.nanmax(angles), np.nanmean(angles), np.nanmin(angles))

    return


def ConvertDisplacementMap2CSv(displacementMapPath, bandNumber, outputPath=None,
                               csvHeader=["xMap(UTM)", "yMap(UTM)", "velocity(m/year)"]):
    """
    Convert displacement map to csv file
    :param displacementMapPath:
    :param bandNumber:
    :param outputPath:
    :return:
    """

    rasterInfo = RT.GetRasterInfo(displacementMapPath)
    print(rasterInfo)

    dataDic = {csvHeader[0]: [], csvHeader[1]: [], csvHeader[2]: []}

    for xPix in np.arange(rasterInfo.get("Dimension")[0]):
        for yPix in np.arange(rasterInfo.get("Dimension")[1]):  # width

            xMap, yMap = RT.Pixel2Map(imageInfo=rasterInfo, x=xPix, y=yPix)
            dataDic[csvHeader[0]].append(xMap)
            dataDic[csvHeader[1]].append(yMap)

            vel = RT.ImageAsArray(imageInfo=rasterInfo, bandNumber=bandNumber)[yPix, xPix]
            dataDic[csvHeader[2]].append(vel)

    df = pd.DataFrame(data=dataDic)
    if outputPath:
        # df.to_csv("/home/cosicorr/0-WorkSpace/meanEW_Velocity_UTM11N_corr.csv")
        df.to_csv(outputPath)
    print(df)
    return


if __name__ == '__main__':
    ewPath = "//home/cosi/2-Data/4-Shisper/2-PostProcessing_S2/1-Shishper_Side/1-EW-NS-velocity/EW_Velocity_Masked/"
    nsPath = "//home/cosi/2-Data/4-Shisper/2-PostProcessing_S2/1-Shishper_Side/1-EW-NS-velocity/NS_Velocity_Masked/"
    outputPath = "//home/cosi/2-Data/4-Shisper/2-PostProcessing_S2/1-Shishper_Side/1-EW-NS-velocity/"
    # CorrelationReport(ewPath=ewPath, nsPath=nsPath, outputFolder=outputPath, reportFileName="CorrReport_WithMask.txt",saveReport=True)

    velocityPath = "//home/cosi/2-Data/4-Shisper/4-Results/0-Final-Results_onlineData/Shisper/2-Velocities-Data-Cube_PC10/"
    # VelocityProfile(inputPath=velocityPath)

    path = "H:\From_CosiWorkStation\\2-Data\\4-Shisper\\3-PostProcessing_Joint_LC_SL\MochowarSide\\2-NS_EW_S2_L8\\NS_targetDate"
    # RenameCorr_SalveDate(inputFolder=path)

    # ComputeMean_Orientation(
    #     pathNS="H:\From_CosiWorkStation\\2-Data\\4-Shisper\\3-PostProcessing_Joint_LC_SL\ShisperSide\\2-S2_and_L8\\0-NS_Masked_Gaussian_Masked_targetName",
    #     pathEW="H:\From_CosiWorkStation\\2-Data\\4-Shisper\\3-PostProcessing_Joint_LC_SL\ShisperSide\\2-S2_and_L8\\0-EW_Masked_Gaussian_Masked_targetName",
    #     outputPath="H:\From_CosiWorkStation\\2-Data\\4-Shisper\\3-PostProcessing_Joint_LC_SL\ShisperSide\\2-S2_and_L8")

    orthoFolder = "/home/cosicorr/0-WorkSpace/Test_Data_ElMayor/Correlation_Set2_test/Test3/OrthoTemp"
    outputfolder = "/home/cosicorr/0-WorkSpace/3D-Correlation_project/Ridgecrest/7p1_3DDA/Sets/WV_Spot_Sets/WV_Spot_sets"
    CreateBatchCorrelationFiles(pathbaseImgs=orthoFolder, pathtargetImgs=orthoFolder, savingPath=outputfolder,
                                corrOptions=["frequency", "128", "128", "32", "32", "8", "8", "2", "0.9", "0", "1"],
                                corrStrategy=4)
