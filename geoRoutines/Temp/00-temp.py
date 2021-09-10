# importing os module
# importing shutil module
import shutil

from geospatialroutine.FilesCommandRoutine import *
from geospatialroutine.Overlap import OverlappingClass
from geospatialroutine.Routine_Lyr import Intersection


class Img:
    def __init__(self, fpPath):
        self.fpPath = fpPath
        self.imgFolder = os.path.basename(fpPath)
        self.preIntersect = []
        self.postIntersect = []


def SelectIntesctImagesWithRefShp(imgsFolder, refShp, destinationFolder=None):
    fpList = ExtractSubfiles(inputdirectory=imgsFolder, fileExtension=[".geojson"])

    fpListIntersect = []
    for fp_ in fpList:
        if Intersection(fp1=refShp, fp2=fp_) != 0:
            fpListIntersect.append(fp_)
    print("------ From {} images only {} images intersect with the input shapefile".format(len(fpList),
                                                                                           len(fpListIntersect)))
    # print(fpListIntersect)
    if destinationFolder:
        for index, fp_ in enumerate(fpListIntersect):
            print("--- Copying:", index + 1, "/", len(fpListIntersect))
            srcPath = os.path.dirname(fp_)
            # List files and directories
            # print("Before copying file:")
            # print(os.listdir(srcPath))
            des = os.path.join(destinationFolder, os.path.basename(srcPath))
            # # Copy the content of
            # # source to destination
            # # using shutil.copy() as parameter
            destination = shutil.copytree(srcPath, des, copy_function=shutil.copy)
            return destination


def Temp(fpPath, fpPreList, fpPostList, thPre=40, thPost=50):
    fpObj = Img(fpPath)

    fpPreList.remove(fpObj.fpPath)
    if len(fpPreList) > 0:
        for fp_ in fpPreList:
            over = OverlappingClass(imgfp1=fpObj.fpPath, imgfp2=fp_, display=False)
            if float(over.overlapPerc) > thPre:
                fpObj.preIntersect.append([fp_, over.overlapPerc])

    for fpPost_ in fpPostList:
        over2 = OverlappingClass(imgfp1=fpObj.fpPath, imgfp2=fpPost_, display=False)
        if float(over2.overlapPerc) > thPost and len(fpObj.preIntersect) > 0:
            for index_, temp_ in enumerate(fpObj.preIntersect):
                over3 = OverlappingClass(imgfp1=temp_[0], imgfp2=fpPost_, display=False)
                if float(over3.overlapPerc) > thPost:
                    fpObj.postIntersect.append([fpPost_, over2.overlapPerc, index_, over3.overlapPerc])
    return fpObj


if __name__ == '__main__':
    preFolderPath = "G:\Ridgecrest\\0-Raw DATA\Post_eq"
    ruptureShp = "G:\Ridgecrest\Main_Rupture_7.1\PolylineFault.shp"
    destinationFolder = "G:\Ridgecrest\\0-Raw DATA\PostIntersectWithRupture"
    # SelectIntesctImagesWithRefShp(imgsFolder=preFolderPath,refShp=ruptureShp)

    fpPreList = ExtractSubfiles(inputdirectory="G:\Ridgecrest\\0-Raw DATA\PreIntersectWithRupture",
                                fileExtension=[".geojson"])
    fpPostList = ExtractSubfiles(inputdirectory="G:\Ridgecrest\\0-Raw DATA\PostIntersectWithRupture",
                                 fileExtension=[".geojson"])
    imgObj = Temp(fpPath=fpPreList[10], fpPostList=fpPostList, fpPreList=fpPreList)
    print("--- ", imgObj.fpPath, "-----")
    print(imgObj.preIntersect)
    print(imgObj.postIntersect)
    #     print(fpObj_.postIntersect)

    # preObjList = []
    # for fp_ in fpPreList:
    #     temp = Img(fp_)
    #     preObjList.append(temp)
    #
    #
    # for index,fpObj_ in enumerate (preObjList):
    #     print(" ---- img:", index +1, " ---------")
    #     fpPreList.remove(fpObj_.fpPath)
    #     if len(fpPreList)>0:
    #         for fp_ in fpPreList:
    #             over = OverlappingClass(imgfp1=fpObj_.fpPath,imgfp2=fp_,display=False)
    #             if float(over.overlapPerc) >40:
    #                 fpObj_.preIntersect.append([fp_,over.overlapPerc])
    #
    #     for fpPost_ in fpPostList:
    #         over2 = OverlappingClass(imgfp1=fpObj_.fpPath, imgfp2=fpPost_, display=False)
    #         if float(over2.overlapPerc) > 50 and len(fpObj_.preIntersect)>0 :
    #             for index_,temp_ in enumerate(fpObj_.preIntersect):
    #                 over3 = OverlappingClass(imgfp1=temp_[0], imgfp2=fpPost_, display=False)
    #                 if float(over3.overlapPerc) > 50:
    #                     fpObj_.postIntersect.append([fpPost_,over2.overlapPerc,index_,over3.overlapPerc])
    #
    #     # print(fpObj_.preIntersect)
    #     # print(fpObj_.postIntersect)
    #     # break
    #
    # # for index, fpObj_ in enumerate(preObjList):
    # #     overlapList = []
    # #     for preOver_ in fpObj_.preIntersect:
    # #         overlapList.append(float(preOver_[1]))
    # #     preMaxIndex = overlapList.index(max(overlapList))
    # #     print(preMaxIndex)
    # #     postOverlapList = []
    # #     for postOver_ in fpObj_.postIntersect:
    # #         if float(postOver_[2]) == preMaxIndex:
    # #             print(postOver_)
    # #             postOverlapList
    # #     break
    # for index, fpObj_ in enumerate(preObjList):
    #     print("--- ",fpObj_.fpPath, "-----")
    #     print(fpObj_.preIntersect)
    #     print(fpObj_.postIntersect)
