import os
import pandas
from pathlib import Path
import requests
import numpy as np
# from bs4 import BeautifulSoup
from geoRoutines.FilesCommandRoutine import GetFilesBasedOnExtension
from bs4 import BeautifulSoup
import geoRoutines.FilesCommandRoutine as FileRT


def download_to_local():
    import logging
    import os
    from google.cloud import storage
    global table_id
    global bucket_name

    bucket_name = 'your-bucket-name'
    prefix = 'your-bucket-directory/'
    dl_dir = 'your-local-directory/'

    # storage_client = storage.Client()
    client = storage.Client("https://storage.googleapis.com")
    # storage_client = storage.Client.from_service_account_json('service-292354190926@gs-project-accounts.iam.gserviceaccount.com')
    print(client)
    # bucket = storage_client.get_bucket(bucket_name=bucket_name)
    # blobs = bucket.list_blobs(prefix=prefix)  # Get list of files
    # for blob in blobs:
    #     filename = blob.name.replace('/', '_')
    #     blob.download_to_filename(dl_dir + filename)  # Download

    # logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.DEBUG)
    # bucket_name = "mybucket"
    # table_id = "shakespeare"
    # storage_client = storage.Client.from_service_account_json("/google-cloud/keyfile/service_account.json")
    # # The “folder” where the files you want to download are
    # folder="/google-cloud/download/{}".format(table_id)
    # delimiter="/"
    # bucket=storage_client.get_bucket(bucket_name)
    # blobs=bucket.list_blobs(prefix=table_id, delimiter=delimiter) #List all objects that satisfy the filter.
    # # Download the file to a destination
    # def download_to_local():
    #  logging.info("File download Started…. Wait for the job to complete.")
    #  # Create this folder locally if not exists
    #  if not os.path.exists(folder):
    #     os.makedirs(folder)
    #  # Iterating through for loop one by one using API call
    #  for blob in blobs:
    #     logging.info("Blobs: {}".format(blob.name))
    #     destination_uri = "{}/{}".format(folder, blob.name)
    #     blob.download_to_filename(destination_uri)
    #     logging.info("Exported {} to {}".format(
    #     blob.name, destination_uri))


def DownloadPlanetData(csvFile=None):
    csvList = GetFilesBasedOnExtension(path=csvFile, filter="*.csv")

    ## Create folders
    for csv_ in csvList:
        dirName = os.path.join(os.path.dirname(csv_), Path(csv_).stem)
        if not os.path.exists(dirName):
            os.makedirs(dirName)

    ## Start the download
    for csv_ in csvList:
        dirName = os.path.join(os.path.dirname(csv_), Path(csv_).stem)
        print("\n----- ", Path(csv_).stem, "\n")
        df = pandas.read_csv(csv_)
        files2Download = df["Link to Tar file"]
        nbTot = len(files2Download)
        sizeList = [float(i[:-2]) for i in df["Size of Tar file"]]
        print("Toal size=", np.sum(sizeList), " GB")

        print("Files to download:\n", files2Download)
        downloadedSize = 0
        for index, file in enumerate(files2Download):
            basename = os.path.basename(file).split("?")[0]
            print(basename)
            urlFile = file
            r = requests.get(urlFile, stream=True)
            savingPath = os.path.join(dirName, basename)
            print("status :", index + 1, "/", nbTot, "size:", sizeList[index], " GB")

            print("Remaining Size :", np.sum(sizeList) - downloadedSize, " GB")
            print(savingPath)

            with open(savingPath, 'wb') as f:
                f.write(r.content)

            downloadedSize += sizeList[index]
    return


def listFD(url, userName, pwd, ext=''):
    page = requests.get(url, auth=(userName, pwd)).text
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]


def main(savingDirectory):
    url = "https://cad4nasa.gsfc.nasa.gov/restricted_data/Avouac/Req892/"
    url = "https://cad4nasa.gsfc.nasa.gov/restricted_data/Avouac/Req892/_cont/"
    url = "https://cad4nasa.gsfc.nasa.gov/restricted_data/Avouac/Req927/"
    url = "https://cad4nasa.gsfc.nasa.gov/restricted_data/Avouac/Req928/"
    url = "https://cad4nasa.gsfc.nasa.gov/restricted_data/Avouac/Req971/"
    url = "https://cad4nasa.gsfc.nasa.gov/restricted_data/Avouac/Req982/"
    url = "https://cad4nasa.gsfc.nasa.gov/restricted_data/Avouac/Req992/"
    # url = "https://cad4nasa.gsfc.nasa.gov/restricted_data/Avouac/Req982/"
    # userName = "2o2o_user"
    userName = "user_2o2l"
    # pwd = "DG_@cce55"
    pwd = "@cce55_DG"
    filesInUrl = listFD(url=url, userName=userName, pwd=pwd, ext="")
    files2Download = []
    alreadydownloaded = FileRT.FilesInDirectory(savingDirectory)
    print(alreadydownloaded)
    # filesinSavingFiles = GetFilesBasedOnExtension(savingDirectory,filter= "*.zip")
    for item in filesInUrl:
        if ".NTF" in item or ".tar" in item or ".zip" in item or ".ntf" in item:
            # if "WV2" in item and "P" in item:

            basename = os.path.basename(item)
            if basename not in alreadydownloaded:
                if "-M1BS" not in basename:
                    temp = os.path.join(savingDirectory, basename)
                    # if  temp not in filesinSavingFiles:
                    files2Download.append(item)

    nbTot = len(files2Download)
    print("Files to download", files2Download, nbTot)

    for index, file in enumerate(files2Download):
        basename = os.path.basename(file)
        urlFile = file
        r = requests.get(urlFile, auth=(userName, pwd), stream=True)
        savingPath = os.path.join(savingDirectory, basename)
        print("status :", index + 1, "/", nbTot)
        print(savingPath)

        with open(savingPath, 'wb') as f:
            f.write(r.content)

    return


if __name__ == '__main__':
    # savingPath = "G:\Kumamoto Eq _2016_Mw7.1"
    # savingPath = "G:\Kumamoto Eq _2016_Mw7.1\WV3"
    # savingPath = "H:\Kumamoto Eq _2016_Mw7.1\WV1"
    # savingPath = "H:\Kumamoto Eq _2016_Mw7.1\WV2"
    # savingPath = "/media/cosicorr/storage/Saif2/NewWV_data"
    savingPath = "//media/cosicorr/storage/Saif/NewData/#982"
    # savingPath = "/media/cosicorr/storage/Saif/NewData/981"
    savingPath = "/media/cosicorr/storage/Saif/NewData/#922_Kamamutu"
    main(savingDirectory=savingPath)
    # DownloadPlanetData(csvFile="/media/storage/Saif/Planet_autoCalibration_project/Ridgecrest")

    # download_to_local()
