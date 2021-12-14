# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2020

from datetime import datetime, timedelta
import os
import geoRoutines.FilesCommandRoutine as FileRT

import time

_start_time = time.time()
def GetTime(seconds):
    sec = timedelta(seconds=seconds)
    d = datetime(1, 1, 1) + sec

    # print("DAYS:HOURS:MIN:SEC")
    print("%d:%d:%d:%d" % (d.day - 1, d.hour, d.minute, d.second))

def tic():
    global _start_time
    _start_time = time.time()


def tac():
    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec, 60)
    (t_hour, t_min) = divmod(t_min, 60)
    print('Time passed: {}hour:{}min:{}sec'.format(t_hour, t_min, t_sec))


def Unzip(outputFolder, zipPath):
    """
    We need to add th .rar and tar cases
    """
    import zipfile
    from pathlib import Path

    if ".zip" in zipPath:
        zip_ref = zipfile.ZipFile(zipPath, 'r')
        folderName = Path(zipPath).stem
        zip_ref.extractall(path=os.path.join(outputFolder, folderName))
        zip_ref.close()
    return


def UnzipBatch(inputFolder, outputFolder=None, compressType="*.zip"):
    from pathlib import Path

    zipFiles = FileRT.GetFilesBasedOnExtension(inputFolder, filter=compressType)
    if outputFolder == None:
        outputFolder = inputFolder
    for index, file_ in enumerate(zipFiles):
        print("*********** Unzipping:" + Path(file_).stem + ": " + str(index + 1) + "/" + str(len(zipFiles))+"*******\n")
        Unzip(outputFolder=outputFolder, zipPath=file_)
        break

    return
