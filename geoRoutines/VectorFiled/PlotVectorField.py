import os
import matplotlib.pyplot as plt
import numpy as np

from geospatialroutine.VectorField.VectorDataLoad import VectorDataLoad_FromRaster


def ReduceVectorField(ewVal, nsVal, windowSize=[25, 25], step=10):
    """
    The idea is to make a copy of EW and NS and fill it with zero,
    Then change the value for each step and window by the mean value
    :param ewVsl:
    :param nsVal:
    :param windowSize:
    :param step:
    :return:
    """

    def crop_center(imgArray, index, windowSize):
        x, y = index
        cropx, cropy = windowSize
        startx = x - cropx // 2
        starty = y - cropy // 2
        return imgArray[startx:x + cropx // 2 + 1, starty:y + cropy // 2 + 1]

    ewValRed = np.zeros(ewVal.shape)
    nsValRed = np.zeros(nsVal.shape)
    nbIndexX = (np.shape(ewVal)[1] - windowSize[0]) // step
    nbIndexY = (np.shape(ewVal)[0] - windowSize[1]) // step
    indexList = []
    for i in range(0, nbIndexY):
        for j in range(0, nbIndexX):
            indexX = windowSize[0] // 2 + i * step
            indexY = windowSize[1] // 2 + j * step
            indexList.append((indexX, indexY))

    for index_ in indexList:
        windValEw = crop_center(imgArray=ewVal, index=index_, windowSize=windowSize)
        ewValRed[index_[0], index_[1]] = np.nanmean(windValEw)
        windValNs = crop_center(imgArray=nsVal, index=index_, windowSize=windowSize)
        nsValRed[index_[0], index_[1]] = np.nanmean(windValNs)

    return ewValRed, nsValRed

def PlotVectorField(ewVal,nsVal,scaleFactor=10,threshold= None):
    M = np.sqrt(ewVal ** 2 + nsVal ** 2)
    print(np.nanmean(M), np.nanmax(M),np.nanstd(M))
    if threshold:
        M = np.ma.masked_where((threshold < M )& (M < 0.01) , M)
        print(np.nanmean(M), np.nanmax(M), np.nanstd(M))
    ##u,v = np.meshgrid(xCoord,yCoord)

    # X, Y = np.arange(0, np.shape(M)[1], 1), np.arange(0, np.shape(M)[0], 1)
    # X, Y = np.meshgrid(X, Y)
    ## M = sqrt(u * u + v * v)  # magnitude
    ## x, y = u, v
    ## im = plt.imshow(M)
    qq = plt.quiver(ewVal, nsVal, M ,scale=scaleFactor)#, cmap=plt.cm.jet,)
    ## plt.colorbar(qq, cmap=plt.cm.jet)
    plt.show()

    return


if __name__ == '__main__':

    path = "C:\\temp"
    ewPath = os.path.join(path, "EW_Cum_norm.tif")
    nsPath = os.path.join(path, "NS_Cum_norm.tif")
    ewVal, nsVal, xCoord, yCoord = VectorDataLoad_FromRaster(ewFile=ewPath, nsFile=nsPath)
    ewValRed, nsValRed = ReduceVectorField(ewVal, nsVal)

    M = np.sqrt(ewValRed**2+nsValRed**2)
    print(np.nanmean(M),np.nanmax(M))
    # u,v = np.meshgrid(xCoord,yCoord)

    X, Y = np.arange(0,np.shape(M)[1] , 1), np.arange(0, np.shape(M)[0], 1)
    X, Y = np.meshgrid(X, Y)
    # M = sqrt(u * u + v * v)  # magnitude
    # x, y = u, v
    im = plt.imshow(M)
    qq = plt.quiver(ewValRed,nsValRed, M, cmap=plt.cm.jet,scale=10)
    plt.colorbar(qq, cmap=plt.cm.jet)
    plt.show()
