import fnmatch
import os
from shutil import rmtree, copyfile
import shutil


# import Correlation_misc as CorrMis
def GetFilesBasedOnExtension(path, filter="*.tif", disp=False):
    import glob
    os.chdir(path)
    filesList = []

    for file in glob.glob(filter):
        filesList.append(os.path.join(path, file))
    filesList.sort()
    if disp:
        print("The list of ", filter, " files are:", filesList, " Ntot=", len(filesList))

    return filesList


def GetFilesBasedOnExtensions(path, filterList=["*.tif", "*.vrt"], disp=False):
    import glob
    os.chdir(path)
    filesList = []

    for filter in filterList:
        for file in glob.glob(filter):
            filesList.append(os.path.join(path, file))
    filesList.sort()
    if disp:
        print("The list of ", filterList, " files are:", filesList, " Ntot=", len(filesList))

    return filesList


def FilesInDirectory(path, exclusionFilter=[], displayFile=False):
    files = os.listdir(path)
    files.sort()
    if exclusionFilter:
        oldFiles = files
        files = []
        for file_ in oldFiles:
            if len(exclusionFilter) > 0:
                check = any(item in file_ for item in exclusionFilter)
                # print(check)
                if not check:
                    # print(file_)
                    files.append(file_)

    if displayFile == True:
        print("The list of  files are:", files, " Ntot=", len(files))
    else:
        print(" Ntot=", len(files))
    files_ =[]
    for file_ in files:
        files_.append(os.path.join(path,file_))
    return files_


def GetOnlyTifImages(path):
    imgs = FilesInDirectory(path)
    imgsPath = []
    for index, img_ in enumerate(imgs):
        if ".tif" in img_ and ".hdr" not in img_ and ".aux.xml" not in img_:
            imgsPath.append(os.path.join(path, img_))
    return imgsPath


def DirectoryAsEmpty(directoryPath):
    files = FilesInDirectory(directoryPath)
    for i in range(len(files)):
        filePath = directoryPath + files[i]
        os.remove(filePath)
    return


def CreateDirectory(directoryPath, folderName, cal=None):
    """
    Create a Folder in the directory.
    Before that this function verify if the name of the folder exist in the directory
    if the folder exist the user will chose either crete new one or delete it
    :param directoryPath: the path where to create the folder : string
    :param folderName: the name of the folder : string
    :return:
    """
    # define the name of the directory to be created

    path = os.path.join(directoryPath, folderName)

    if os.path.exists(path):
        if cal == None:
            print("<< %s >> folder already exist type y to delete and recreate new one or type n :  " % folderName)
            cal = input("Type y/n:")
            if cal == "y":
                try:
                    rmtree(path)
                    os.makedirs(path, exist_ok=True)
                except OSError:
                    print("Creation of the directory %s failed " % path)
                else:
                    print("Successfully created the directory <<%s>> " % path)
                    return (path )

            if cal == "n":
                return (path )
        if cal == "y":
            try:
                rmtree(path)
                os.makedirs(path, exist_ok=True)
            except OSError:
                print("Creation of the directory %s failed " % path)
            else:
                print("Successfully created the directory <<%s>> " % path)
                return (path )
        if cal == "n":
            return (path )
    else:
        try:

            os.makedirs(path, exist_ok=True)

        except OSError:
            print("Creation of the directory %s failed " % path)

        else:
            print("Successfully created the directory << %s >> " % path)
            return (path )


def CreateTxtOfFiles(inputFolder, outputFileName="ListofImgs.txt"):
    files = FilesInDirectory(path=inputFolder)
    with open(os.path.join(inputFolder, outputFileName), 'w') as f:
        for file_ in files:
            f.write("%s\n" % file_)

    return


def DeleteSubsets(inputFilePath, refTxtFile):
    ##### Version 0.01 ##### To be continued
    """
    this function will delete the subset that does not exit in the refFile
    :param pathRefFile:
    :return:
    """
    ## Read all imgs in subset folder
    files = FilesInDirectory(path=inputFilePath)

    for file_ in files:
        if ".tif" in file_ and file_ not in open(refTxtFile).read():
            print(False)
            print(file_)
            os.remove(inputFilePath + file_)


def CompareTwoFiles(file1, file2):
    lineList1 = [line.rstrip('\n') for line in open(file1)]
    lineList2 = [line.rstrip('\n') for line in open(file2)]
    for f in lineList1:
        if f not in lineList2:
            print(f)
    return


def UncompressFile(compressedFilePath, output=None):
    import tarfile
    import zipfile
    if tarfile.is_tarfile(compressedFilePath):
        file = tarfile.open(compressedFilePath)

        file.list()
        print(file.getmembers())
        print(file.getnames())
        file.close()

    if zipfile.is_zipfile(compressedFilePath):
        with zipfile.ZipFile(compressedFilePath, 'r') as zip_ref:
            directory_to_extract_to = output
            zip_ref.extractall(directory_to_extract_to)
    return


def UncompressBatch(directoryInput, directoryOutput):
    filesList = GetFilesBasedOnExtension(path=directoryInput, filter="*.zip")
    nbTot = len(filesList)
    for index, file_ in enumerate(filesList):
        baseName = os.path.basename(file_)[0:-4]
        print("status :", index + 1, "/", nbTot)
        UncompressFile(compressedFilePath=file_, output=os.path.join(directoryOutput, baseName))

    return


def CopyFile(inputFilePath, outputFolder, overWrite=True):
    outputFilePath = os.path.join(outputFolder, os.path.basename(inputFilePath))
    files = os.listdir(outputFolder)
    if os.path.basename(inputFilePath) not in files:
        copyfile(src=inputFilePath, dst=outputFilePath)
    else:
        print(os.path.basename(inputFilePath), ": exist in destination folder!!")
        if overWrite:
            print("-- replacing with the new file !! --")
            copyfile(src=inputFilePath, dst=outputFilePath)
        else:
            print("-- Keeping the old one !!--")

    return outputFilePath


def Copyfiles(inputdirectory, destinationDirectory, filter=".NTF"):
    for root, dirs, files in os.walk(inputdirectory):
        for name in files:
            if name.endswith((filter)):
                ntfFile = os.path.join(root, name)
                print(ntfFile)
                copyfile(src=ntfFile, dst=os.path.join(destinationDirectory, name))

    return


def LocateFile(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''

    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def ExtractSubfiles(inputdirectory, fileExtension=[".NTF"], disp=False):
    """

    """
    filesList = []
    for root, dirs, files in os.walk(inputdirectory):
        for name in files:
            if any(name.endswith(ele) for ele in fileExtension):
                file = os.path.join(root, name)
                filesList.append(file)
    if disp:
        print("Subfiles:", filesList, "NbTot=", len(filesList))
    return filesList


def ContentFolderDelete(folder, exception=None):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        if filename != exception:
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))
    return


def ContentOfFolder(folderPath):
    # https://thispointer.com/python-how-to-get-list-of-files-in-directory-and-sub-directories/
    fileListNames = []
    dirPathList = []
    dirNamesList = []

    for (dirpath, dirnames, filenames) in os.walk(folderPath):
        fileListNames.append(filenames)
        dirPathList.append(dirpath)
        # print(dirnames,dirpath,filenames)
        # print(dirnames)

    fileList = []
    for dirPath_, list_ in zip(dirPathList, fileListNames):
        for name_ in list_:
            fileList.append(os.path.join(dirPath_, name_))
    #
    print("===> #files:",len(fileList))

    return fileList,dirPathList


if __name__ == '__main__':
    tmpPath = "/home/cosi/2-Data/4-Shishper/3-PostProcessing_Joint_LC_SL/ShishperSide/2-S2_and_L8/"
    file1 = tmpPath + "log_not_filtered.txt"
    file2 = tmpPath + "log_filtered.txt"

    # CompareTwoFiles(file1, file2)

    # UncompressFile(compressedFilePath="C:\\temp\\13JUL19WV031300019JUL13191105-P1BS-503388194040_01_P001_________AAE_0AAAAABPAIW0.tar")

    # UncompressBatch(directoryInput="G:\Data\Ridgecrest\WV_Data_ridgecrest_Chris\cont_",
    #                 directoryOutput="G:\Data\Ridgecrest\WV_Data_ridgecrest_Chris\cont_\\Unzipped")

    # Copyfiles(inputdirectory="G:\Data\Ridgecrest\WV_Data_ridgecrest_Chris\cont_\\Unzipped\Post_eq\\",
    #           destinationDirectory="G:\Data\Ridgecrest\WV_Data_ridgecrest_Chris\cont_\\Unzipped\All_NTF_Post_pan\\")
    rsmList = ["/media/cosicorr/storage/Saif/1-Ridgecrest/Data/1-WV/P1BS_sorted/Before_eq/Before_eq/2017-11-30-WV2/503916437010_01_P004_PAN"]
    ContentOfFolder(rsmList[0])