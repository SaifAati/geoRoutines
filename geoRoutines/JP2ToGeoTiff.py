# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2020

import argparse
import glob
import os

from osgeo import gdal


def GetFilesBasedOnExtension(directoryPath, filter="*.jp2"):
    """
    :param directoryPath: input directory : string
    :param filter: filter to select files: string
    :return: list of files: list
    """
    os.chdir(directoryPath)
    filesList = []
    for file in glob.glob(filter):
        filesList.append(os.path.join(directoryPath, file))
    filesList.sort()
    print("The list of ", filter, " files are:", filesList, " Ntot=", len(filesList))

    return filesList


def RunConv(args):
    inputFolder = args.inputPath
    outputFolder = args.outputPath

    files = GetFilesBasedOnExtension(directoryPath=inputFolder)
    for index, inputRasterPath in enumerate(files):
        # imgRasterInfo = GetRasterInfo(inputRaster=imgPath)
        if not outputFolder:
            newRasterPath = os.path.join(os.path.dirname(inputRasterPath),
                                         os.path.basename(inputRasterPath)[:-4] + ".tif")
            print("newRasterPath=", newRasterPath)
        else:
            newRasterPath = os.path.join(outputFolder,
                                         os.path.basename(inputRasterPath)[:-4] + ".tif")
            print("newRasterPath=", newRasterPath)

        srcDS = gdal.Open(inputRasterPath)
        gdal.Translate(newRasterPath, srcDS, format="GTiff", outputType=gdal.GDT_Float64)
    return


def main():
    parser = argparse.ArgumentParser(description="Convert JP2 raster file to GeoTiff ")
    parser.add_argument("-inputPath", help="Directory path containing all JP2 files", type=str,
                        required=True)
    parser.add_argument("-outputPath",
                        help="Directory path where all converted data will be stored, if this argument kept empty the images will be stored in the same input directory ",
                        type=str, required=False)
    parser.set_defaults(func=RunConv)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
