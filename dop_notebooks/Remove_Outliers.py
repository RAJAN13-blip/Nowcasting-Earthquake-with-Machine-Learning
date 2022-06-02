#!/opt/local/bin python
#import sys
#sys.path.reverse()

    #   Earthquake Methods library of methods and functions
    #   
    #   This code base collects the methods and functions used to make
    #   plots and maps of earthquake data and activity
    #
    ######################################################################

import sys
import matplotlib
import matplotlib.mlab as mlab
#from matplotlib.pyplot import figure, show
import numpy as np
from numpy import *
from mpl_toolkits.basemap import Basemap
from array import array
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as patches

import datetime
import dateutil.parser

import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import os

import math
from math import exp

import MCUtilities
from MCUtilities import *
import MCCalcMethods

from matplotlib import cm

import http.client
from urllib.error import HTTPError

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


#   -------------------------------------------------------

def reject_outliers(data_lat, data_lng, sd_factor):

    #   This method removes lat-lng points that are far from the median of the burst
    
    median_lat = np.median(data_lat)
    median_lng = np.median(data_lng)
    
    great_circle_distance   =   []
    data_lat_removed        =   []
    data_lng_removed        =   []
    
    for i in range(len(data_lat)):
        lat1 = median_lat
        lng1 = median_lng
        lat2 = data_lat[i]
        lng2 = data_lng[i]
        dist = compute_great_circle_distance(lat1, lng1, lat2, lng2)
        great_circle_distance.append(dist)
        
    median_gcdist = np.median(great_circle_distance)
    
    for i in range(len(data_lat)):
        if great_circle_distance[i] < sd_factor * median_gcdist:
            data_lat_removed.append(data_lat[i])
            data_lng_removed.append(data_lng[i])

    return data_lat_removed, data_lng_removed
    
def compute_great_circle_distance(lat1, lng1, lat2, lng2):

    # Build an array of x-y values in degrees defining a circle of the required radius

    pic = 3.1415926535/180.0
    Radius = 6371.0

    lat1 = float(lat1) * pic
    lng1 = float(lng1) * pic
    lat2 = float(lat2) * pic
    lng2 = float(lng2) * pic
    
    delta_lng = lng1 - lng2
    
    delta_radians = math.sin(lat1)*math.sin(lat2) + math.cos(lat1)*math.cos(lat2)*cos(delta_lng)
    delta_radians = math.acos(delta_radians)
    
    great_circle_distance = delta_radians * Radius

    return great_circle_distance

    #.................................................................

    
if __name__ == '__main__':

    data_lat = [33.9253333, 33.9271667, 33.9325, 33.92, 33.9553333, 33.9256667, 33.9075, 33.92, 33.913, 33.9161667, 33.9136667, 33.9383333,
             33.9248333, 33.9403333, 33.9376667, 33.9225, 33.9036667, 33.8976667, 33.9191667, 33.9046667, 33.9243333, 33.9093333,
             33.9056667, 33.912, 33.7608333, 33.9301667, 33.8988333, 33.9011667, 33.9031667, 33.9038333, 33.9038333, 33.9011667, 
             33.9613333, 33.9571667, 33.9618333, 34.1188333, 33.9118333, 33.911, 33.9456667, 33.9006667, 33.9326667, 33.903, 33.908,
             33.9093333, 33.9256667, 33.9063333, 33.9565, 33.901, 33.9, 33.9573333, 33.9243333, 33.908, 33.9205, 33.9168333, 33.919, 33.896]

    data_lng = [-117.9341667, -117.9188333, -117.9158333, -117.9316667, -117.8911667, -117.9255, -117.954, -117.9303333, -117.9273333,
             -117.9313333, -117.9326667, -117.898, -117.925, -117.9068333, -117.9128333, -117.9003333, -117.9625, -117.9576667,
             -117.9356667, -117.9383333, -117.9178333, -117.9451667, -117.9543333, -117.9488333, -117.5196667, -117.9181667,
             -117.9626667, -117.9536667, -117.968, -117.9551667, -117.9418333, -117.943, -117.8923333, -117.9035, -117.8945,
             -118.4888333, -117.9258333, -117.9358333, -117.9065, -117.9546667, -117.9141667, -117.9495, -117.9475, -117.9411667,
             -117.9288333, -117.9521667, -117.8876667, -117.9578333, -117.9583333, -117.899, -117.9115, -117.9441667, -117.9333333,
             -117.9393333, -117.9375, -117.9556667]

    median_lat = np.median(data_lat)
    median_lng = np.median(data_lng)

    print('')
    print('Length data_lat, Median Latitude: ', len(data_lat), median_lat)
    print('')
    print('Original Latitude Data: ', data_lat)
    
#    reject_outliers(data_lat)

    print('')
    print('Length data_lng, Median Longitude: ', len(data_lat), median_lng)
    print('')
    print('Original Longitude Data: ', data_lng)
    print('')

    
#     great_circle_distance   =   []
#     data_lat_removed        =   []
#     data_lng_removed        =   []
#     
# 
#     
#     for i in range(len(data_lat)):
#         lat1 = median_lat
#         lng1 = median_lng
#         lat2 = data_lat[i]
#         lng2 = data_lng[i]
#         dist = compute_great_circle_distance(lat1, lng1, lat2, lng2)
#         great_circle_distance.append(dist)
#         
#     median_gcdist = np.median(great_circle_distance)
#     
#     for i in range(len(data_lat)):
#         if great_circle_distance[i] < 5.* median_gcdist:
#             data_lat_removed.append(data_lat[i])
#             data_lng_removed.append(data_lng[i])
#             

    sd_factor = 5.
    data_lat_removed, data_lng_removed = reject_outliers(data_lat, data_lng, sd_factor)
    
    print('')
    print('Length data_lat_removed: ', len(data_lat_removed))
    print('data_lat_removed: ', data_lat_removed)
    print('')
    
    print('')
    print('Length data_lng_removed: ', len(data_lng_removed))
    print('data_lng_removed: ', data_lng_removed)
    print('')
    
    print('Removed ', len(data_lat)-len(data_lat_removed), ' Points')
    
    
    
    
    
    
    
    
    
    
    
    
    
    