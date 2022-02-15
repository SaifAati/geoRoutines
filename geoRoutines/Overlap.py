# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2020
import geojson
import geopandas
import matplotlib.pyplot as plt
from geoRoutines.RandomColors import GenerateColors
# import geoRoutines.Remove_from_PublicRelease.Routine as RT
from geoRoutines.georoutines import ComputeEpsg
class cOverlapping:
    ## Note need to use the intersection in routine lyer
    def __init__(self, imgfp1, imgfp2, display=False):
        self.imgfp1 = imgfp1
        self.imgfp2 = imgfp2
        self.intersection = 0
        self.overlapPerc = 0

        self.Intersection()
        self.Area()
        if display:
            self.Visualise_overlapping()

    def Visualise_overlapping(self):


        fig, ax = plt.subplots()
        colors = GenerateColors(N=len(self.fpDF_List), pastel_factor=0.5)
        if self.intersection != 0:
            self.interDF.plot(ax=ax, alpha=0.5, cmap='tab10')
        self.fpDF_List[0].plot(ax=ax, facecolor='none', edgecolor=colors[0],lw=4)
        self.fpDF_List[1].plot(ax=ax, facecolor='none', edgecolor=colors[1],lw=4)
        plt.show()

        return

    def Intersection(self):

        fpList = [self.imgfp1, self.imgfp2]
        self.fpDF_List = []
        for fp_ in fpList:
            fpDataFrame = geopandas.read_file(fp_)
            self.fpDF_List.append(fpDataFrame)

        self.interDF = geopandas.overlay(self.fpDF_List[0], self.fpDF_List[1], how='intersection')
        if self.interDF.index.start == self.interDF.index.stop:
            intersection = 0
        else:
            intersection = geojson.Feature(geometry=self.interDF.loc[0, "geometry"])

        self.intersection = intersection
        return self.intersection, self.interDF

    def Area(self):
        if self.intersection != 0:
            dfCopy = self.interDF.copy()
            coords = self.intersection['geometry']['coordinates']
            epsg = ComputeEpsg(lon=coords[0][0][0], lat=coords[0][0][1])
            epsg_ = "epsg:" + str(epsg)
            # dfCopy = dfCopy.to_crs({'init': epsg_})
            dfCopy = dfCopy.to_crs("epsg:"+str(epsg))
            self.inter_area = (dfCopy['geometry'].area / 10 ** 6)[0]

            # fpDFCopy = self.fpDF_List[0].copy().to_crs({'init': epsg_})
            fpDFCopy = self.fpDF_List[0].copy().to_crs("epsg:"+str(epsg))
            fpArea = (fpDFCopy["geometry"].area / 10 ** 6)[0]
            self.overlapPerc = float('%1.1f' % ((self.inter_area / fpArea) * 100))
            # print("Intersection area (Km^2): ", self.inter_area, " Overlapping %: ", self.overlapPerc)
        else:
            self.overlapPerc = 0
        return self.overlapPerc