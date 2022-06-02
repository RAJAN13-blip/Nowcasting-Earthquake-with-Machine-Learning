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

import scipy.stats
from scipy.stats import poisson

import datetime
import dateutil.parser

import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import os

import math
from math import exp

import BSTUtilities
from BSTUtilities import *
import BSTCalcMethods

from matplotlib import cm

import http.client
from urllib.error import HTTPError

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from tabulate import tabulate

import time

    ###########################################################################
    ###########################################################################

    # Get the catalog parameters
    
def get_current_date_time():

    now = datetime.datetime.now()

    print("Current year: %d" % now.year)
    print("Current month: %d" % now.month)
    print("Current day: %d" % now.day)
    print("Current hour: %d" % now.hour)
    print("Current minute: %d" % now.minute)
    print("Current second: %d" % now.second)
    
    slash = '/'
    colon = ':'
    
    year  = str(now.year)
    month = str(now.month)
    day   = str(now.day)
    
    if now.month < 10:
        month = '0'+ str(now.month)
        
    if now.day < 10:
        day = '0'+ str(now.day)
    
    current_date = year + slash + month + slash + day
    current_time = str(now.hour) + colon + str(now.minute) + colon + str(now.second)
    
    return (current_date, current_time, year, month, day)

def write_to_file(output_file_name,catalog):

    #   This function writes the earthquake data to file
    
#   output_file_name = "USGS_WorldWide.catalog"
    output_file = open(output_file_name, "a")
    
    #############################################################

    #   Write output record file

    # Format the data

    date_string = ' '
    time_string = ' '
    ts = 0.
    
    i=-1
    for line in catalog:
        i += 1
        if i > 0:
            line_decode = line.decode('UTF-8')  #   Convert from Byte literal to String
            items = line_decode.strip().split(',')
#           items = line.strip().split(',')
        #   print items

            date = dateutil.parser.parse(items[0].split('T')[0])     # gives year-month-day (space) hours:minutes:seconds
    #       date = dateutil.parser.parse(items[0])    <-- gives year-month-day only
            ts = date.year + float(date.strftime("%j"))/366

            lat                 = items[1]
            lon                 = items[2]
            dep                 = items[3]
            mag                 = items[4]
            dsb                 = items[0].split('T')[0].split('-')
            date_string         = dsb[0]+'/'+ dsb[1] + '/' + dsb[2]
            time_string         = items[0].split('T')[1]
            time_string         = time_string[:-1]
            
    #   These next checks are in case depths are missing, or depths & mags are listed as 'Unk'

            if dep == 'Unk':
                dep = '0.0'

            if mag == 'Unk':
                mag = '0.0'

            if mag == 'ML':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mb':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mw':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mc':
                mag = items[4]
                dep = '0.0'

            if mag == 'Md':
                mag = items[4]
                dep = '0.0'

            if mag == 'Ms':
                mag = items[4]
                dep = '0.0'
                
            output_file.write("%s\t%s\t%f\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))
                
    output_file.close()
                
    return None

def download_worldwide_catalog(NELat, NELng, SWLat, SWLng,Magnitude,begin_date,end_date,output_file_name):

    #   For instructions, refer to:  https://earthquake.usgs.gov/fdsnws/event/1/

    data = {
    
    # Adjust these  -   CA-NV
        "minmagnitude": Magnitude,
        "minlatitude": SWLat,
        "maxlatitude": NELat,
        "minlongitude": SWLng,
        "maxlongitude": NELng,
        "mindepth": -10.0,         # Leaving these depth params in leads to missing some earthquake events
        "maxdepth": 1000.0,


    # Date parameters
#        "starttime": "2070/01/01",
        "starttime": begin_date,
        "endtime": end_date,


    # Leave these
        "eventtype": "earthquake",
        "format": "csv",
        "orderby": "time-asc",
    }
    
    block_size = 10000
    event_offset = 1
#     event_count = [0 for i in range(1)]
    
   
    #   First do a count of events
    
    url = "https://earthquake.usgs.gov/fdsnws/event/1/count?"
    params = urllib.parse.urlencode(data)
    query_string = url + str(params)
    
    error_code = 0
    
    try:
        response_count = urllib.request.urlopen(query_string)
        event_count = response_count.readlines()
#         print('event_count: ', event_count)
        event_count = event_count[0].decode('UTF-8')
#         print('event_count: ', event_count)
        number_events = int(event_count)
        
    except:
        error_code = 1
        number_events = 0
        print('')
        print('1 Download barfed, Error Code: ', error_code)
        pass

    print('')
    print('Number of Events: ', number_events)
    n_blocks = int(event_count)/block_size        #   Number of complete blocks
    n_blocks = int(n_blocks)
    print('Number of complete blocks of size ' + str(block_size) + ' =', n_blocks)
    print('')
    
    for i_block in range(0,n_blocks):
    
    #   -------------------------------------------------------------
    
        event_offset = i_block * block_size + 1
    
        data.update({'offset':event_offset})
        data.update({'limit':block_size})
    
        url = "https://earthquake.usgs.gov/fdsnws/event/1/query?"
        params = urllib.parse.urlencode(data)
        query_string = url + str(params)
    
        error_code = 0
        
        try:
            response = urllib.request.urlopen(query_string)
            catalog = response.readlines()
        except:
            error_code = 1
            print('')
            print('2 Download barfed, Error Code: ', error_code)
            pass
        
        print('Block Number: ', i_block)
        
    #   -------------------------------------------------------------
    
        write_to_file(output_file_name,catalog)
                
    residual_events = number_events
    
    if number_events > block_size:
        residual_events = number_events%(n_blocks*block_size)
        
    if residual_events > 0:
        if n_blocks == 0:
            event_offset = 1
        if n_blocks > 0:
            event_offset = n_blocks * block_size + 1
        data.update({'offset':event_offset})
        data.update({'limit':block_size})
    
        url = "https://earthquake.usgs.gov/fdsnws/event/1/query?"
        params = urllib.parse.urlencode(data)
        query_string = url + str(params)
    
        error_code = 0
        try:
            response = urllib.request.urlopen(query_string)
            catalog = response.readlines()
        except:
            error_code = 1
            print('')
            print('3 Download barfed, Error Code: ', error_code)
            pass
        
        if error_code == 0:
            write_to_file(output_file_name,catalog)
    
    return None
    
def get_worldwide_catalog(maxlat, minlat, maxlng, minlng, completeness_mag, start_date):

    time_seconds = 0.0
    restart = 'NO'          #   We want to output a new file
    output_file_name = "USGS_WorldWide.catalog"
    
    if restart == 'NO':     #   Erase the prior contents of the file
        output_file = open(output_file_name, "w")   #   This statement dumps any previous file
        output_file.close()
        
    current_date, current_time, current_year, current_month, current_day = get_current_date_time()
    print('current_date, current_time: ', current_date, current_time)
    
    start_year = int(start_date.split('/')[0])
    end_year = int(current_year)    
    number_years = end_year - start_year + 1
#
#   -------------------------------------------------------
#     
#   Write the WorldWide catalog by appending data year by year

    print('')
    print('------------------------------------------')
    print('')
    print('Building WorldWide Catalog for M>' + str(completeness_mag))
    print('')
    print('------------------------------------------')
    print('')
    
    for i in range (0,number_years):
        begin_year = str(int(i) + start_year)
        begin_date = str(begin_year)+ '/01/01'
        end_date   = str(int(begin_year)+1) + '/01/01'   
        last_date  = str(begin_year) + '/12/31'
        
        print('')
        print('')
        print('Begin Date to End Date: ', begin_date, ' to ', last_date) 

#       download_worldwide_catalog(85.0,179.99,-85.0,-179.99,completeness_mag, begin_date,end_date,output_file_name)
        download_worldwide_catalog(maxlat,maxlng,minlat,minlng,completeness_mag, begin_date,end_date,output_file_name)
        
        time.sleep(time_seconds)
        
    return

    #################################################################
    
def get_master_catalog(NELat, NELng, SWLat, SWLng, Magnitude):

#   This code will use shapely.py to filter the worldwide catalog
#       into a rectangular region

    settings_params = BSTUtilities.get_settings()
    region_type = settings_params[0]
    location = settings_params[3]

    data_string = 'Building catalog for region around ' + location + '...'
    
#    print ''
    print('------------------------------------------')
    print('')
    print(data_string)
    print('')
    print('------------------------------------------')
#    print ''

    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])
        settings_params[1] = completeness_mag
        settings_params[2] = earthquake_depth
        save_settings(settings_params)

    output_file_name = "USGS_master.catalog"
    output_file = open(output_file_name, "w")
    output_file.close()

    number_polygon_vertices = 4

    #   Construct the string of polygon vertices.  Note that the order is lat, long pairs
    
    vertex_lat = []
    vertex_lng = []
    
    #   Order of vertices of large rectangular region:  NW, NE, SE, SW
    
    vertex_lat.append(NELat)
    vertex_lat.append(NELat)
    vertex_lat.append(SWLat)
    vertex_lat.append(SWLat)
    
    vertex_lng.append(SWLng)
    vertex_lng.append(NELng)
    vertex_lng.append(NELng)
    vertex_lng.append(SWLng)
    
    point_list = []

    for i in range(number_polygon_vertices):
        point_list.append((float(vertex_lat[i]),float(vertex_lng[i])))
    
    polygon = Polygon(point_list)
#       

#     input_file_name     = "USGS_WorldWide.catalog"
#     input_file          =  open(input_file_name, "r")
# 
#     output_file_name = "USGS_master.catalog"
#     output_file = open(output_file_name, "w")
#     
#     for line in input_file:
#         items = line.strip().split()
#         dep    = items[6]
#         mag    = items[5]
#         eq_lat = items[4]
#         eq_lng = items[3]
#         
#         point = Point((float(eq_lat),float(eq_lng)))
#         
#         if (float(dep) <= float(earthquake_depth) and mag >= float(completeness_mag) and polygon.contains(point) == True):
# #           print items[0],items[1],items[2],items[3],items[4],items[5],items[6]
#             print >> output_file, items[0],items[1],items[2],items[3],items[4],items[5],items[6]
#         
#     output_file.close()
#     input_file.close()

    return
    
def get_circle_master_catalog(Latitude,Longitude,Radius,Magnitude):

#   This code will use shapely.py to filter the worldwide catalog
#       into a rectangular region

    settings_params = BSTUtilities.get_settings()
    region_type = settings_params[0]
    location = settings_params[3]

    data_string = 'Building circle catalog for region around ' + location + '...'
    
    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])
        settings_params[1] = completeness_mag
        settings_params[2] = earthquake_depth
        save_settings(settings_params)

    output_file_name = "USGS_circle.master.catalog"
    output_file = open(output_file_name, "w")
    output_file.close()

    #   Compute vertices that define the circle around the city
    
    lng_circle_dg, lat_circle_dg = BSTUtilities.createCircleAroundWithRadius(Latitude, Longitude, Radius)
    
    number_polygon_vertices = len(lng_circle_dg)
    
    point_list = []

    for i in range(number_polygon_vertices):
        point_list.append((float(lat_circle_dg[i]),float(lng_circle_dg[i])))
    
    polygon = Polygon(point_list)
    
#   print point_list
#       

    input_file_name     = "USGS_WorldWide.catalog"
    input_file          =  open(input_file_name, "r")

    output_file_name = "USGS_circle.master.catalog"
    output_file = open(output_file_name, "w")
    
    for line in input_file:
        items = line.strip().split()
        dep    = items[6]
        mag    = items[5]
        eq_lat = items[4]
        eq_lng = items[3]
        
        point = Point((float(eq_lat),float(eq_lng)))
        
        if (float(dep) <= float(earthquake_depth) and float(mag) >= float(completeness_mag) and polygon.contains(point) == True):
#           print items[0],items[1],items[2],items[3],items[4],items[5],items[6]
            print(items[0],items[1],items[2],items[3],items[4],items[5],items[6], file=output_file)
        
    output_file.close()
    input_file.close()

    return None
    
    #############################################################

    # Set the catalog parameters for circular region

def get_circle_catalog(Circle_Lat, Circle_Lng, Radius_float, Rebuild_flag, \
        circle_catalog_date_start, circle_catalog_date_end):

    settings_params = BSTUtilities.get_settings()
    region_type = settings_params[0]
    location = settings_params[3]

    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])
        
#     print 'completeness_mag: ', completeness_mag
#     print ''

    if Rebuild_flag == 'OFF':
        data_string = 'Building catalog for  M>' + str(completeness_mag) + ' for small circle around ' + location + '...'
        
    if Rebuild_flag == 'ON':
            data_string = 'Rebuilding catalog for M>' + str(completeness_mag) + ' for small circle around ' + location 
    
    print('------------------------------------------')
    print('')
    print(data_string)
    print('')
    print('------------------------------------------')
#     print ''

    #############################################################
    
    get_circle_master_catalog(Circle_Lat, Circle_Lng, Radius_float, completeness_mag)
    
    # Open the master output file to read
    input_circle_master = open("USGS_circle.master.catalog", "r")

    #   Write output record file

    # Open the output file
    output = open("USGS_circle.catalog", "w")

    # Format the data

    mag_array   =   []
    date_array  =   []
    time_array  =   []
    year_array  =   []
    depth_array =   []
    lat_array   =   []
    lng_array   =   []

    i=-1
    for line in input_circle_master:
        i+=1
        items = line.strip().split()

        try:
        
            lat                 = items[4]
            lon                 = items[3]
            dep                 = items[6]
            mag                 = items[5]
            date_string         = items[0]
            time_string         = items[1]
            ts                  = items[2]

            mag_array.append(mag)           #   List of magnitudes
            date_array.append(date_string)
            time_array.append(time_string)
            year_array.append(ts)
            depth_array.append(dep)
            lat_array.append(lat)
            lng_array.append(lon)
            
            if float(mag) >= 6.49:
                print('')
                print('Data for events >= 6.5: ')
                print(items)
                print('')

            if float(dep) <= float(earthquake_depth) and float(ts) >= circle_catalog_date_start and \
                    float(ts) <= circle_catalog_date_end:
                output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))

        except:
            pass

    # Finalize the output file
    input_circle_master.close()
    output.close()

    return mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array
    
def read_circle_catalog():

    mag_array   =   []
    date_array  =   []
    time_array  =   []
    year_array  =   []
    depth_array =   []
    lat_array   =   []
    lng_array   =   []

    data_file = open("USGS_circle.catalog","r")
    
    for line in data_file:
        items = line.strip().split()
        
        try:
        
            lat                 = items[4]
            lon                 = items[3]
            dep                 = items[6]
            mag                 = items[5]
            date_string         = items[0]
            time_string         = items[1]
            ts                  = items[2]

            mag_array.append(mag)           #   List of magnitudes
            date_array.append(date_string)
            time_array.append(time_string)
            year_array.append(ts)
            depth_array.append(dep)
            lat_array.append(lat)
            lng_array.append(lon)
            
        except:
            pass
    
    data_file.close()  

    return mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array

    ######################################################################
    
def get_circle_data(Location):

#   Read Circle_Location and data from "Settings_File.txt"

    settings_params = BSTUtilities.get_settings()
    default_settings_params = settings_params

    region_type = settings_params[0]
    completeness_mag    =   float(settings_params[1])
    earthquake_depth    =   float(settings_params[2])

    if region_type == 'Circle':     #   If still a circle, proceed
        Circle_Location     =   settings_params[3]
        small_radius        =   float(settings_params[6])

    if region_type != 'Circle':     #   If not a circle, reset settings_params to null
        Circle_Location     =   'None'
        Circle_Lat          =   0.0
        Circle_Lng          =   0.0
        Radius_float        =   0.0

        settings_params     =   []
        settings_params.append(region_type)
        settings_params.append(completeness_mag)
        settings_params.append(earthquake_depth)
        settings_params.append(Circle_Location)
        settings_params.append(Circle_Lat)
        settings_params.append(Circle_Lng)
        settings_params.append(Radius_float)

   #
   #    Get the data for the circle: Read pre-defined locations file
   #

    input_file = open("city_locations.txt", "r")
    i=0
    for line in input_file:
        i +=  1
    input_file.close()  # Put the file pointer back at top

    number_circle_locations = i

    # Create arrays of length i filled with zeros

    Circle_Location_file   = ["" for x in range(i)]

    CircleLat       = np.zeros(i)
    CircleLng       = np.zeros(i)
    Radius          = np.zeros(i)

    input_file = open("city_locations.txt", "r")

    i=-1
    for line in input_file:
        i+=1
        line    = line.strip()
        items   = line.split(',')

    #   Fill arrays

        items_array = np.asarray(items)

        Circle_Location_file[i]   = items_array[0]
        CircleLat[i]       = float(items_array[1])
        CircleLng[i]       = float(items_array[2])
        Radius[i]          = small_radius

    input_file.close()  # Put the file pointer back at top


    Current_Circle_Location = Circle_Location

#    if new_location_resp == 'y':
#        print ' Enter location (Case sensitive: Overrides previous parameter set):'
#        Circle_Location = raw_input()
#        region_type = 'Circle'

#    if new_location_resp != 'y':
#        Circle_Location = Current_Circle_Location
#        region_type = 'Circle'

    circle_location_flag = 0
    counter = 0

    while (circle_location_flag == 0):
        for j in range(0,number_circle_locations):

            counter += 1

            if Circle_Location == Circle_Location_file[j]:

                CircleCenterLat = CircleLat[j]
                CircleCenterLng = CircleLng[j]
                CircleRadius    = Radius[j]

                Circle_Lat         = float(CircleCenterLat)
                Circle_Lng         = float(CircleCenterLng)
                Radius_float       = float(CircleRadius)
 
                settings_params[0]  =   region_type
                settings_params[1]  =   completeness_mag
                settings_params[2]  =   earthquake_depth
                settings_params[3]  =   Circle_Location
                settings_params[4]  =   Circle_Lat
                settings_params[5]  =   Circle_Lng
                settings_params[6]  =   Radius_float
          
                circle_location_flag = 1

#   Write new value of Circle_Location to "Settings_File.txt"

    save_settings(settings_params)

    return (Circle_Location, Circle_Lat, Circle_Lng, Radius_float)

    #.................................................................


    
