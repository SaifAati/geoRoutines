
import rasterio
from rasterio.plot import show
from matplotlib import pyplot


def PlotWithRasterio(rasterPath):
    """
    https://rasterio.readthedocs.io/en/latest/topics/plotting.html
    :param rasterPath:
    :return:
    """
    src = rasterio.open(rasterPath)
    show(src.read(), transform=src.transform)
    pyplot.show()

def PlotHistogram(rasterPath):
    from rasterio.plot import show_hist
    src = rasterio.open(rasterPath)
    show_hist( src, bins = 200, lw = 0.0, stacked = False, alpha = 0.3,
    histtype = 'stepfilled', title = "Histogram")

if __name__ == '__main__':
    PlotHistogram(rasterPath="/home/cosicorr/0-WorkSpace/AutoCal_Planet_workSpace/SkySat_Corr_Ridgecrest/20190121_182150_ssc2_u0001_pansharpened.tif_VS_20190706_212238_ssc8_u0001_pansharpened.tif")