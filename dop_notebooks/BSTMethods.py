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

# from file_read_backwards import FileReadBackwards

    ######################################################################
    
    ###########################################################################
    ###########################################################################
    
    #   Miscellaneous Methods
    
    ###########################################################################
    ###########################################################################

    # Set the catalog parameters
    
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
            items = line.strip().split(',')
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
    
    block_size = 20000
    event_offset = 1
    
   
    #   First do a count of events
    
    url = "https://earthquake.usgs.gov/fdsnws/event/1/count?"
    params = urllib.parse.urlencode(data)
    query_string = url + str(params)
    
    error_code = 0
    try:
        response_count = urllib.request.urlopen(query_string)
        event_count = response_count.readlines()
        number_events = int(event_count[0])
    except:
        error_code = 1
        number_events = 0
        print('')
        print('Download barfed, Error Code: ', error_code)
        pass
    
    print('')
    print('Number of Events: ', number_events)
    n_blocks = int(event_count[0])/block_size        #   Number of complete blocks
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
            print('Download barfed, Error Code: ', error_code)
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
            print('Download barfed, Error Code: ', error_code)
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
        
    return

    #################################################################
    
def get_master_catalog(NELat, NELng, SWLat, SWLng, Magnitude):

#   This code will use shapely.py to filter the worldwide catalog
#       into a rectangular region

    settings_params = get_settings()
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

    input_file_name     = "USGS_WorldWide.catalog"
    input_file          =  open(input_file_name, "r")

    output_file_name = "USGS_master.catalog"
    output_file = open(output_file_name, "w")
    
    for line in input_file:
        items = line.strip().split()
        dep    = items[6]
        mag    = items[5]
        eq_lat = items[4]
        eq_lng = items[3]
        
        point = Point((float(eq_lat),float(eq_lng)))
        
        if (float(dep) <= float(earthquake_depth) and mag >= float(completeness_mag) and polygon.contains(point) == True):
#           print items[0],items[1],items[2],items[3],items[4],items[5],items[6]
            print(items[0],items[1],items[2],items[3],items[4],items[5],items[6], file=output_file)
        
    output_file.close()
    input_file.close()

    return
    
def get_circle_master_catalog(Latitude,Longitude,Radius,Magnitude):

#   This code will use shapely.py to filter the worldwide catalog
#       into a rectangular region

    settings_params = get_settings()
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

def get_catalog(NELat, NELng, SWLat, SWLng, MagLo):

    settings_params = get_settings()
    region_type = settings_params[0]
    location = settings_params[3]

    data_string = 'Building regional catalog for large region centered on ' + location + '...'
    
    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])
        settings_params[1] = completeness_mag
        settings_params[2] = earthquake_depth
        save_settings(settings_params)

    #############################################################
    
        #   This method now calls the get_master_catalog() method which writes
    #       the master catalog.  All the other get_catalog-type methods
    #       will open and filter the master

    #    completeness_mag = 2.99
    
    get_master_catalog(NELat, NELng, SWLat, SWLng, completeness_mag)

    #   Write output record file

    # Open the master output file to read
    output_master = open("USGS_master.catalog", "r")

    #   Now write the output file.  So first, open the output file for the last events
#    output = open("USGS_master.catalog", "w")

    #   Now write the output_last file.  So first, open the output file for the last events
    output_last = open("USGS_last.catalog", "w")

    # Format the data

    i=-1
    for line in output_master:
        i+=1
        items = line.strip().split()
        
        lat                 = items[4]
        lon                 = items[3]
        dep                 = items[6]
        mag                 = items[5]
        date_string         = items[0]
        time_string         = items[1]
        ts                  = items[2]
        
    #   These next checks are in case depths are missing, or depths & mags are listed as 'Unk'

    #   print lat,lon,dep,mag,date_string,time_string

#       if (float(dep) <= float(earthquake_depth)):
#            output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))

        if (float(dep) <= float(earthquake_depth)):
            output_last.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))

    #   We use only the small earthquakes after the last large earthquake. So if this event is not the last
    #       large event, close the file, move the pointer to the top, and re-open the file for writing

        if (float(mag) >= float(MagLo)):
            output_last.close()
            output_last = open("USGS_last.catalog", "w")           

    # Finalize the output file
    output_master.close()
#    output.close()
    output_last.close()

    print(' ')
    print('     Read in Catalog Data, now Building Working File')

    #############################################################

    #   Now Write the Output Working File

    #############################################################

    #   First count number of lines in master data file that we just downloaded above   

    data_file = open("USGS_master.catalog", "r")

    #   Now open the working file where we will re-write the data

    j=0
    for line in data_file:
        j+=1
    h=j

    data_file.close()
    data_file = open("USGS_master.catalog", "r") # Put the file pointer back at the top

    working_file = open("EQ_Working_Catalog.txt", "w")
    print('')
    
    i=0
    j=0
    
    for line in data_file:
        j+=1
        items = line.strip().split()
    #
        percent_complete = 100.0*float(j)/float(h) 
        print('     Percent File Completed: ', "{0:.2f}".format(percent_complete),'%', "                                 \r", end=' ')   

    #
    #   Remember that at this point, all the below variables are string variables, not floats
    #

        date_string         = items[0]
        time_string         = items[1]
        ts                  = items[2]
        lon                 = items[3]
        lat                 = items[4]
        mag                 = items[5]
        dep                 = items[6]

        lat     = 0.0001*float(int(10000.0*float(lat + '000000')))

        if float(mag) >= MagLo:
            i+=1
            event = str(i)
            if i <10:
                event = '0'+str(i)

            working_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (event, date_string, time_string, ts, lon, lat, mag, dep))
            sys.stdout.flush()
            
    # Finalize the output file
    
    print('')

    working_file.close()
    data_file.close()

    return None
    
    ######################################################################

    # Set the catalog parameters for circular region

def get_circle_catalog(Circle_Lat, Circle_Lng, Radius_float, MagLo, Rebuild_flag, \
        circle_catalog_date_start, circle_catalog_date_end):

    settings_params = get_settings()
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
    output_circle_master = open("USGS_circle.master.catalog", "r")

    #   Write output record file

    # Open the output file
    output = open("USGS_circle.catalog", "w")
    output_last_circle = open("USGS_last.circle.catalog", "w")      

    # Format the data

    mag_array   =   []
    date_array  =   []
    time_array  =   []
    year_array  =   []
    depth_array =   []
    lat_array   =   []
    lng_array   =   []

    i=-1
    for line in output_circle_master:
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

    #       print lat,lon,dep,mag,date_string,time_string

            if float(dep) <= float(earthquake_depth) and float(ts) >= circle_catalog_date:
                output_last_circle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))

    #   We use only the small earthquakes after the last large earthquake. So if this event is not the last
    #       large event, close the file, move the pointer to the top, and re-open the file for writing

            if (float(mag) >= float(MagLo)):
                output_last_circle.close()
                output_last_circle = open("USGS_last.circle.catalog", "w")      

        except:
           pass

    # Finalize the output file
    output_circle_master.close()
    output.close()
    output_last_circle.close()

    #############################################################

    #   Now Write the Output Working File
    
    print(' ')
    print('     Read in Catalog Data, now Building Working File')

    #############################################################

    #   First count number of lines in master data file that we just downloaded above

    data_file = open("USGS_circle.catalog", "r")

    #   Now open the working file where we will re-write the data

    j = 0
    for line in data_file:
        j += 1
    h = j

    data_file.close()
    data_file = open("USGS_circle.catalog",
                     "r")  # Put the file pointer back at the top

    os.system("rm EQ_Working_Circle_Catalog.txt")

    working_file = open("EQ_Working_Circle_Catalog.txt", "w")
    
    i = 0
    j = 0
    print('')
    for line in data_file:
        j += 1
        items = line.strip().split()


        percent_complete = 100.0 * float(j) / float(h)
        print('     Percent File Completed: ', "{0:.2f}".format(
            percent_complete), '%', "                                 \r", end=' ')

        #
        #   Remember that at this point, all the below variables are string variables, not floats
        #

        date_string = items[0]
        time_string = items[1]
        ts = items[2]
        lon = items[3]
        lat = items[4]
        mag = items[5]
        dep = items[6]

        lat = 0.0001 * float(int(10000.0 * float(lat + '000000')))

        if float(mag) >= MagLo:
            i += 1
            event = str(i)
            if i < 10:
                event = '0' + str(i)

            working_file.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (event, date_string, time_string, ts, lon, lat, mag, dep))
            sys.stdout.flush()
            
    print('')

    # Finalize the output file
    
    working_file.close()
    data_file.close()

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

    settings_params = get_settings()
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

def get_polygon_data(Location):

#   Read Polygon_Location and data from "Settings_File.txt"

    settings_params = get_settings()
    default_settings_params = settings_params

    region_type = settings_params[0]
    completeness_mag    =   float(settings_params[1])
    earthquake_depth    =   float(settings_params[2])

    polygon_data = {}

    if region_type == 'Polygon':     #   If still a Polygon, proceed
        Polygon_Location     =   settings_params[3]

    if region_type != 'Polygon':     #   If not a Polygon, reset settings_params to null
        Polygon_Location    =   'None'

    polygon_data[0] = ['None', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0',]    # Dummy List = value of dict for index[0]
    default_vertex_points = polygon_data[0]


   #
   #    Get the data for the Polygon: Read pre-defined locations file
   #

    input_file = open("polygonlocations.txt", "r")
    i=0
    for line in input_file:
        i +=  1
    input_file.close()  # Put the file pointer back at top

    number_polygon_locations = i

    # Create arrays of length i filled with zeros

    Polygon_Location_List   = ["" for x in range(i)]

    input_file = open("polygonlocations.txt", "r")

    polygon_data = {}

    i=-1
    for line in input_file:
        i+=1
        line    = line.strip()
        items   = line.split(',')
        Polygon_Location_List[i] = items[0]

#   print items for testing purposes

        vertex_points = []

        for j in range(0,len(items)):
            vertex_points.append(items[j])

        polygon_data[i] = vertex_points

    input_file.close()  # Put the file pointer back at top

    print(' ')
    print('     Current pre-defined Locations are: ')
    print(' ')

    for i in range(0,number_polygon_locations):
        print('        ', Polygon_Location_List[i])

    print(' ')
    print('     Current region is: ')
    print('        ', Location)
    print(' ')
    print('     Current polygon location is: ', Polygon_Location)
    print(' ')
    print('     Pick a new polygon location? (y/n):')
    new_location_resp   =   input()

    Current_Polygon_Location = Polygon_Location

    if new_location_resp == 'y':
        print(' Enter location (Case sensitive: Overrides previous parameter set):')
        Polygon_Location = input()
        region_type = 'Polygon'
        print('Polygon_Location: ', Polygon_Location)

    if new_location_resp != 'y':
        region_type = 'Polygon'
        print(' ')
        print(' Then we use the current region (which could be the dummy region)')
        print(' ')

    polygon_location_flag = 0
    counter = 0
    settings_params = []

    while (polygon_location_flag == 0):

        for i in range(0,number_polygon_locations):

            counter += 1

            if Polygon_Location == polygon_data[i][0]:
                polygon_location_flag = 1

                settings_params.append(region_type)
                settings_params.append(completeness_mag)
                settings_params.append(earthquake_depth)

                vertex_points = polygon_data[i]

                for j in range(0,len(vertex_points)):
                    settings_params.append(vertex_points[j])

        if polygon_location_flag == 0 and counter >= number_polygon_locations:
            polygon_location_flag = 1
            settings_params = default_settings_params
            vertex_points = default_vertex_points
            print(' ')
            print('     Invalid location, try again...')
            print(' ')
            print('     (Press any key to continue)')
            resp = input()
            print(' ')

#   Write new value of Polygon_Location to "Settings_File.txt"

    save_settings(settings_params)

    return (vertex_points)

    #.................................................................

def freqmagRegion(NELat, NELng, SWLat, SWLng, MagLo, Location, max_mag):

    print_data = 'TRUE'

    settings_params = get_settings()
    region_type = settings_params[0]
    Circle_Location = settings_params[3]

    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])

    #   Open the complete data file and collect the earthquakes for plotting

    data_file = open("USGS_master.catalog", "r")     #  This is all the data and all the small events

    #   Find the number of lines in the file

    i=0
    for line in data_file:
        i       +=  1

    data_file.close()  # Put the file pointer back at top

    number_eqs = i

    # Create arrays of length i filled with zeros

    mag                 =   zeros(number_eqs)

    # Bins for Number-magnitude plot

    min_mag = completeness_mag
    bin_diff= 0.1
    
    number_mag_bins     =   (max_mag - min_mag) / bin_diff + 5      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    data_file = open("USGS_master.catalog", "r")     #  This is all the data and all the small events

    i=-1

    for line in data_file:

        i+= 1
        items = line.strip().split()
    #
    #   Remember that at this point, all the below variables are string variables, not floats
    #
        mag[i]          = float(items[5])
        bin_number      = int((mag[i]-min_mag)/bin_diff + 1 )
        if mag[i] >= completeness_mag:
            freq_mag_bins_pdf[bin_number]    +=  1

    data_file.close()                                               # Close the data file
    
    #.................................................................

    #   Tabulate the number of events in each bin

    for i in range(0,range_mag_bins):
        for j in range(i,range_mag_bins):                           # Loop over all bins to compute the GR cumulative SDF
            freq_mag_bins_sdf[i] += freq_mag_bins_pdf[j]
    print('')
    
#   print freq_mag_bins_sdf

    number_nonzero_bins=0
    for i in range(0,range_mag_bins):                               # Find the number of nonzero bins
        if freq_mag_bins_sdf[i] > 0:
            number_nonzero_bins+=1

    range_bins = int(number_nonzero_bins)                         # Find number of nonzero bins

    log_freq_mag_bins   =   zeros(number_nonzero_bins)              # Define the log freq-mag array

    for i in range(0,range_bins):
        if freq_mag_bins_sdf[i] > 0.0:
            log_freq_mag_bins[i] = -100.0                               # To get rid of it.
            log_freq_mag_bins[i] = math.log10(freq_mag_bins_sdf[i])     # Take the log-base-10 of the survivor function

    mag_bins  =   zeros(number_nonzero_bins)
    for i in range(0,range_bins):
        mag_bins[i]=min_mag + float(i)*bin_diff

    for i in range(0,number_mag_bins):
        mag_array[i]=min_mag + float(i)*bin_diff

    #.................................................................

#    if print_data == 'TRUE':
#        print ''
#        print mag_bins
#        print ''
#        print freq_mag_bins_pdf
#        print ''
#        print freq_mag_bins_sdf
#        print ''
#        print log_freq_mag_bins
#        print ''

    #.................................................................

    return mag_bins, log_freq_mag_bins, number_nonzero_bins

    ######################################################################

def freqmagCircle(NELat, NELng, SWLat, SWLng, MagLo, Location, max_mag):

    print_data = 'TRUE'

    settings_params = get_settings()
    region_type = settings_params[0]
    Circle_Location = settings_params[3]

    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])

    #   Open the complete data file and collect the earthquakes for plotting

    data_file = open("USGS_circle.catalog", "r")     #  This is all the data and all the small events

    #   Find the number of lines in the file

    i=0
    for line in data_file:
        i       +=  1

    data_file.close()  # Put the file pointer back at top

    number_eqs = i

    # Create arrays of length i filled with zeros

    mag                 =   zeros(number_eqs)

    # Bins for Number-magnitude plot

    min_mag = completeness_mag
    bin_diff= 0.1

    number_mag_bins     =   (max_mag - min_mag) / bin_diff + 5      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    data_file = open("USGS_circle.catalog", "r")     #  This is all the data and all the small events

    i=-1

    for line in data_file:

        i+= 1
        items = line.strip().split()
    #
    #   Remember that at this point, all the below variables are string variables, not floats
    #
        mag[i]          = float(items[5])
        bin_number      = int((mag[i]-min_mag)/bin_diff + 1 )
        if mag[i] >= completeness_mag:
            freq_mag_bins_pdf[bin_number]    +=  1

    data_file.close()                                               # Close the data file

    #.................................................................

    #   Tabulate the number of events in each bin

    for i in range(0,range_mag_bins):
        for j in range(i,range_mag_bins):                           # Loop over all bins to compute the GR cumulative SDF
            freq_mag_bins_sdf[i] += freq_mag_bins_pdf[j]
    print('')

    number_nonzero_bins=0
    for i in range(0,range_mag_bins):                               # Find the number of nonzero bins
        if freq_mag_bins_sdf[i] > 0:
            number_nonzero_bins+=1

    range_bins = int(number_nonzero_bins)                         # Find number of nonzero bins

    log_freq_mag_bins   =   zeros(number_nonzero_bins)              # Define the log freq-mag array

    for i in range(0,range_bins):
        if freq_mag_bins_sdf[i] > 0.0:
            log_freq_mag_bins[i] = -100.0                               # To get rid of it.
            log_freq_mag_bins[i] = math.log10(freq_mag_bins_sdf[i])     # Take the log-base-10 of the survivor function

    mag_bins  =   zeros(number_nonzero_bins)
    for i in range(0,range_bins):
        mag_bins[i]=min_mag + float(i)*bin_diff

    for i in range(0,number_mag_bins):
        mag_array[i]=min_mag + float(i)*bin_diff

 
    #.................................................................

    return mag_bins, log_freq_mag_bins, number_nonzero_bins

    ######################################################################

def freqmagLastCircle(NELat, NELng, SWLat, SWLng, MagLo, Location, max_mag):

    print_data = 'TRUE'

    settings_params = get_settings()
    region_type = settings_params[0]
    Circle_Location = settings_params[3]

    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])

    #   Open the complete data file and collect the earthquakes for plotting

    data_file = open("USGS_last.circle.catalog", "r")     #  This is all the data and all the small events

    #   Find the number of lines in the file

    i=0
    for line in data_file:
        i       +=  1

    data_file.close()  # Put the file pointer back at top

    number_eqs = i

    # Create arrays of length i filled with zeros

    mag                 =   zeros(number_eqs)

    # Bins for Number-magnitude plot

    min_mag = completeness_mag
    bin_diff= 0.1

    number_mag_bins     =   (max_mag - min_mag) / bin_diff + 5      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    data_file = open("USGS_last.circle.catalog", "r")     #  This is all the data and all the small events

    i=-1

    for line in data_file:

        i+= 1
        items = line.strip().split()
    #
    #   Remember that at this point, all the below variables are string variables, not floats
    #
        mag[i]          = float(items[5])
        bin_number      = int((mag[i]-min_mag)/bin_diff + 1 )
        if mag[i] >= completeness_mag:
            freq_mag_bins_pdf[bin_number]    +=  1

    data_file.close()                                               # Close the data file

    #.................................................................

    #   Tabulate the number of events in each bin

    for i in range(0,range_mag_bins):
        for j in range(i,range_mag_bins):                           # Loop over all bins to compute the GR cumulative SDF
            freq_mag_bins_sdf[i] += freq_mag_bins_pdf[j]
    print('')

    number_nonzero_bins=0
    for i in range(0,range_mag_bins):                               # Find the number of nonzero bins
        if freq_mag_bins_sdf[i] > 0:
            number_nonzero_bins+=1

    range_bins = int(number_nonzero_bins)                         # Find number of nonzero bins

    log_freq_mag_bins   =   zeros(number_nonzero_bins)              # Define the log freq-mag array

    for i in range(0,range_bins):
        if freq_mag_bins_sdf[i] > 0.0:
            log_freq_mag_bins[i] = -100.0                               # To get rid of it.
            log_freq_mag_bins[i] = math.log10(freq_mag_bins_sdf[i])     # Take the log-base-10 of the survivor function

    mag_bins  =   zeros(number_nonzero_bins)
    for i in range(0,range_bins):
        mag_bins[i]=min_mag + float(i)*bin_diff

    for i in range(0,number_mag_bins):
        mag_array[i]=min_mag + float(i)*bin_diff

    #.................................................................

#    if print_data == 'TRUE':
#        print ''
#        print mag_bins
#        print ''
#        print freq_mag_bins_pdf
#        print ''
#        print freq_mag_bins_sdf
#        print ''
#        print log_freq_mag_bins
#        print ''

    #.................................................................

    return mag_bins, log_freq_mag_bins, number_nonzero_bins

    ######################################################################

def freqmagLast(NELat, NELng, SWLat, SWLng, MagLo, Location, max_mag):

    print_data = 'TRUE'

    settings_params = get_settings()
    region_type = settings_params[0]
    Circle_Location = settings_params[3]

    if region_type == 'Circle':
        completeness_mag = float(settings_params[1])
        earthquake_depth = float(settings_params[2])

    #   Open the complete data file and collect the earthquakes for plotting

    data_file = open("USGS_last.catalog", "r")     #  This is all the data and all the small events

    #   Find the number of lines in the file

    i=0
    for line in data_file:
        i       +=  1

    data_file.close()  # Put the file pointer back at top

    number_eqs = i

    # Create arrays of length i filled with zeros

    mag                 =   zeros(number_eqs)

    # Bins for Number-magnitude plot

    min_mag = completeness_mag
    bin_diff= 0.1

    number_mag_bins     =   (max_mag - min_mag) / bin_diff + 5      #   Bin size = 0.1.  Assume min mag of interest is 3.0
    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   zeros(number_mag_bins)
    freq_mag_bins_sdf       =   zeros(number_mag_bins)
    freq_mag_pdf_working    =   zeros(number_mag_bins)
    mag_array               =   zeros(number_mag_bins)  

    data_file = open("USGS_last.catalog", "r")     #  This is all the data and all the small events

    i=-1

    for line in data_file:

        i+= 1
        items = line.strip().split()
    #
    #   Remember that at this point, all the below variables are string variables, not floats
    #
        mag[i]          = float(items[5])
        bin_number      = int((mag[i]-min_mag)/bin_diff + 1 )
        if mag[i] >= completeness_mag:
            freq_mag_bins_pdf[bin_number]    +=  1

    data_file.close()                                               # Close the data file

    #.................................................................

    #   Tabulate the number of events in each bin

    for i in range(0,range_mag_bins):
        for j in range(i,range_mag_bins):                           # Loop over all bins to compute the GR cumulative SDF
            freq_mag_bins_sdf[i] += freq_mag_bins_pdf[j]
    print('')

    number_nonzero_bins=0
    for i in range(0,range_mag_bins):                               # Find the number of nonzero bins
        if freq_mag_bins_sdf[i] > 0:
            number_nonzero_bins+=1

    range_bins = int(number_nonzero_bins)                         # Find number of nonzero bins

    log_freq_mag_bins   =   zeros(number_nonzero_bins)              # Define the log freq-mag array

    for i in range(0,range_bins):
        if freq_mag_bins_sdf[i] > 0.0:
            log_freq_mag_bins[i] = -100.0                               # To get rid of it.
            log_freq_mag_bins[i] = math.log10(freq_mag_bins_sdf[i])     # Take the log-base-10 of the survivor function

    mag_bins  =   zeros(number_nonzero_bins)
    for i in range(0,range_bins):
        mag_bins[i]=min_mag + float(i)*bin_diff

    for i in range(0,number_mag_bins):
        mag_array[i]=min_mag + float(i)*bin_diff

    #.................................................................

    return mag_bins, log_freq_mag_bins, number_nonzero_bins
    
    ###########################################################################
    
def write_correlation_circle_catalog(City_Latitude, City_Longitude, City_Radius, EQ_Latitude, EQ_Longitude,\
            earthquake_depth, MagLo, completeness_mag):
    
    #   First estimate the correlation length in the region from the current circle catalog
    
    Radius_Correlation = City_Radius 
    Radius_Region = 2.0*City_Radius
    Correlation_Scale = 1.0*(Radius_Region - City_Radius)

    #   Compute vertices that define the circle around the city
    
    lng_circle_dg, lat_circle_dg = BSTUtilities.createCircleAroundWithRadius(City_Latitude, City_Longitude, Radius_Region)
    
    number_polygon_vertices = len(lng_circle_dg)
    
    point_list = []

    for i in range(number_polygon_vertices):
        point_list.append((float(lat_circle_dg[i]),float(lng_circle_dg[i])))
    
    polygon = Polygon(point_list)
    
    stop_date = '1999/10/10'
    stop_flag = 0
    
    input_file = open('USGS_master.catalog','r') 

    output_file_name = "USGS_circle.correlation_temp.catalog"
    output_file = open(output_file_name, "w")
    
    for line in input_file:
        items = line.strip().split()
        eq_date = items[0]
        dep     = float(items[6])
        mag     = float(items[5])
        eq_lat  = float(items[4])
        eq_lng  = float(items[3])
        
        point = Point((float(eq_lat),float(eq_lng)))
        
        if eq_date == stop_date:
            stop_flag = 1
        
        if (float(dep) <= float(earthquake_depth) and float(mag) >= float(completeness_mag) and polygon.contains(point) == True\
                    and stop_flag == 0):

    #   In the next bit, we print to file only the EQs since the last large EQ (including the last large EQ)
    
            if float(mag) >= float(MagLo):

                great_circle_radius = BSTUtilities.compute_great_circle_distance(City_Latitude, City_Longitude,eq_lat, eq_lng)
                if great_circle_radius < City_Radius:
                    output_file.close()
                    output_file =  open(output_file_name, "w")  #   Re-initialize file
                
                print(items[0],items[1],items[2],items[3],items[4],items[5],items[6], file=output_file)
            else:
                great_circle_radius = BSTUtilities.compute_great_circle_distance(City_Latitude, City_Longitude,eq_lat, eq_lng)
                if great_circle_radius <= City_Radius:
                    print(items[0],items[1],items[2],items[3],items[4],items[5],items[6], file=output_file)
                else:
                    if np.random.uniform() < math.exp(- (great_circle_radius - City_Radius)/(Correlation_Scale)):
                        print(items[0],items[1],items[2],items[3],items[4],items[5],items[6], file=output_file)
            
    output_file.close()
    input_file.close()
    
    data_file = open('USGS_circle.correlation_temp.catalog','r')
    items_temp = []

    for line in data_file:
        items = line.strip()
        items_temp.append(items)
        
    data_file.close()
    
    os.system("rm USGS_circle.correlation_temp.catalog")
    
    new_output_file = open('USGS_circle.correlation.catalog','w')
    
    items_reversed = list(reversed(items_temp))
    
    for i in range(len(items_temp)):
        print(items_reversed[i], file=new_output_file)
        
    new_output_file.close()
    
    return Radius_Correlation, Radius_Region
    
    ###########################################################################
    
def get_regional_rate(completeness_mag, earthquake_depth, cutoff_start_year):

    number_earthquakes = 0

    input_file = open("USGS_WorldWide.catalog", "r")

    for line in input_file:
        items = line.strip().split()

        date                = items[0]
        time                = items[1]
        year                = items[2]
        dep                 = items[6]
        mag                 = items[5]
            
        if float(mag) >= float(completeness_mag) and float(dep) <= earthquake_depth and float(year) >= cutoff_start_year:
            number_earthquakes += 1
                
    input_file.close()

    time_interval = (float(year) - cutoff_start_year)*365.0
    
    regional_rate = float(number_earthquakes)/time_interval
    
    return regional_rate

    ###########################################################################
    ###########################################################################
    
    #   Plotting Methods
    
    ###########################################################################
    ###########################################################################

def map_epicenters(NELat, NELng, SWLat, SWLng, MagLo, Location):
   
    lat_0_center = 0.5*(NELat + SWLat)
    lon_0_center = 0.5*(NELng + SWLng)


    llcrnrlat = SWLat
    llcrnrlon  = SWLng
    urcrnrlat  = NELat
    urcrnrlon  = NELng

    m = Basemap(projection='cyl',llcrnrlat=SWLat, urcrnrlat=NELat,
            llcrnrlon=SWLng, urcrnrlon=NELng, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,8),labels=[0,0,0,1])
    m.drawparallels(np.arange(-90,90,4),labels=[1,0,0,0])

    #   m.shadedrelief()
    #   m.bluemarble()
    #   m.etopo()
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    #m.drawrivers()

    m.etopo()

    settings_params = get_settings()

    region_type = settings_params[0]

    #   -------------------------------------------

    if region_type == 'Circle':
        Polygon_Location = 'None'

        Circle_Location = settings_params[3]
        earthquake_depth = float(settings_params[2])

        if Circle_Location != 'None':
            earthquake_depth = float(settings_params[2])
            CircleCenterLat = float(settings_params[4])
            CircleCenterLng = float(settings_params[5])
            CircleRadius    = float(settings_params[6])

         #   -------------------------------------------

    if region_type == 'Polygon':
        Circle_Location = 'None'

        number_polygon_vertices = (len(settings_params)-4)/2

        Polygon_Location = settings_params[3]
        earthquake_depth = float(settings_params[2])

        if Polygon_Location != 'None':
            earthquake_depth = float(settings_params[2])

        x_poly = []
        y_poly = []

        for j in range(0,number_polygon_vertices):
            y_poly.append(settings_params[4+2*j])
            x_poly.append(settings_params[5+2*j])

    #   Close the polygon

        y_poly.append(settings_params[4])
        x_poly.append(settings_params[5])

        for i in range(0,len(x_poly)):
            x_poly[i] = float(x_poly[i])
            y_poly[i] = float(y_poly[i])

  
    #   -------------------------------------------


    #   Open input file and find the number of lines in the file

    input_file = open("EQ_Working_Catalog.txt", "r")

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs=i

    #   Open input file again

    input_file = open("EQ_Working_Catalog.txt", "r")

    # Loop over lines and extract variables of interest

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(i)]
    yrs                 =   zeros(i)
    lng                 =   zeros(i)
    lat                 =   zeros(i)
    mag                 =   zeros(i)
    dep                 =   zeros(i)
    time_string         =   ["" for x in range(i)]
    date_string         =   ["" for x in range(i)]

    # Loop over lines and extract variables of interest

    i=-1
    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    input_file.close()  # Put the file pointer back at top

    z_limit = len(lat)
    range_limit = z_limit - 1

    print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' around ' + Location
    print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
    print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
    print_text_4 = '     at depths < ' + str(settings_params[2]) + ' km'
    #print ' '
    #print print_text_1
    #print ' '
    #print print_text_2
    #print print_text_3
    #print ' '

    #.................................................................

    for z in range(z_limit):
        x,y = m(lng,lat)

    #   -----------------------
        if mag[z] >= 5.0:
            mark_size = 3
        if mag[z] >= 5.5:
            mark_size = 4
        if mag[z] >= 6:
            mark_size = 6
        if mag[z] >= 6.5:
            mark_size = 8
        if mag[z] >= 7:
            mark_size = 10
        if mag[z] >= 7.5:
            mark_size = 12
        if mag[z] >= 8.0:
            mark_size = 14
        if mag[z] >= 8.5:
            mark_size = 16
    #   ----------------------
            
    #   m.plot(x[z], y[z], "ro", mfc='none', mec='r', markeredgewidth=2.0, ms=mark_size[z])
        m.plot(x[z], y[z], "ro", ms=mark_size)
        m.plot(x[z], y[z], "o", ms=mark_size, fillstyle='none', color='k', lw=0.4)
        

    if (Circle_Location != 'None'):
#        x_circle_dg, y_circle_dg = draw_big_circle(CircleCenterLat, CircleCenterLng, CircleRadius)
        x_circle_dg, y_circle_dg = BSTUtilities.createCircleAroundWithRadius(CircleCenterLat, CircleCenterLng, CircleRadius)
        
        m.plot(x_circle_dg, y_circle_dg, "b-", lw=1.3)

    if (Polygon_Location != 'None'):
        m.plot(x_poly, y_poly, "b-", lw=1.6)


    #   ------------------------------------------------------
    #   Plot 150 km circle at latest event

    lat_latest = lat[i]
    lng_latest = lng[i]
    x,y = m(lng_latest,lat_latest)

    #   Draw a circle of radius 200 km

    Radius = 200.0  # km
    mark_size = int(155.0 * Radius/200.0)   #  approx at mid latitudes

    #   mark_size = 155
    #   m.plot(x,y, "o", mfc='none', ms=mark_size)

    #   ------------------------------------------------------

    sfv     =   (NELat-SWLat)/12.0
    sfh     =   (NELng-SWLng)/16.0

    sfh     =   sfh * 0.9
    sfht    =   sfh * 1.1  # move the text over a bit more

    lat=[llcrnrlat - 0.30*sfv, llcrnrlat + 0.25*sfv,llcrnrlat + 0.80*sfv, llcrnrlat + 1.30*sfv, llcrnrlat + 1.80*sfv, llcrnrlat + 2.40*sfv, llcrnrlat + 3.1*sfv, llcrnrlat + 3.95*sfv]
    lng=[urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht, urcrnrlon + 1.2*sfht]
    x,y = m(lng,lat)

    latc=[llcrnrlat - 0.25*sfv, llcrnrlat + 0.30*sfv,llcrnrlat + 0.85*sfv, llcrnrlat + 1.375*sfv, llcrnrlat + 1.90*sfv, llcrnrlat + 2.50*sfv, llcrnrlat + 3.2*sfv, llcrnrlat + 4.0*sfv]
    lngc=[urcrnrlon + 0.72*sfh, urcrnrlon + 0.72*sfh, urcrnrlon + 0.73*sfh, urcrnrlon + 0.74*sfh, urcrnrlon + 0.76*sfh, urcrnrlon + 0.78*sfh, urcrnrlon + 0.80*sfh, urcrnrlon + 0.80*sfh]
    xc,yc = m(lngc,latc)

    mag_value=['5.0+','5.5+','6+','6.5+','7+','7.5+', '8+', '8.5+']

    mark_size = [3, 4, 6, 8, 10, 12, 14, 16]

    for z in range(8):
        plt.text(x[z],y[z],mag_value[z], fontsize=10,)
        m.plot(xc[z], yc[z], "ro", ms=mark_size[z], clip_on=False)
        m.plot(xc[z], yc[z], "ro", ms=mark_size[z], clip_on=False, fillstyle='none', color='k', lw=0.4)

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]

    if region_type == 'Circle':

        if Circle_Location == 'None':
            SupTitle_text = 'Earthquake Epicenters for M>' + str(MagLo) + ' near ' + Location

        if Circle_Location != 'None':
            Circle_Location_actual = Circle_Location.split('-')
            SupTitle_text = 'Earthquakes M>' + str(MagLo) + ' near ' + City_Location

    if region_type == 'Polygon':

        if Polygon_Location == 'None':
            SupTitle_text = 'Earthquake Epicenters for M>' + str(MagLo) + ' in ' + Location + ' at Depth < ' + str(earthquake_depth) +  'km'

        if Polygon_Location != 'None':
            Polygon_Location_actual = Polygon_Location.split('-')
            SupTitle_text = 'Earthquakes M>' + str(MagLo) + ' in ' + Location + ' at Depth < ' + str(earthquake_depth) + ' km'

    plt.suptitle(SupTitle_text, fontsize=14)

    Title_text = 'From: ' + date_string[0] + ' To: ' + date_string[range_limit] + ' R<' + str(int(CircleRadius)) \
            + ' km, D<'+ str(int(earthquake_depth)) +' km' 
    plt.title(Title_text, fontsize=12)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')

    figure_name = './Data/' + Circle_Location + '_Seismicity.png'

    matplotlib.pyplot.savefig(figure_name,dpi=600)
    #matplotlib.pyplot.savefig('./Data/Seismicity-Map.jpg')
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

    #plt.show()

    return None

    ######################################################################
    
def map_epicenters_local(NELat_local, NELng_local, SWLat_local, SWLng_local, completeness_mag, MagLo, Location):
   
    lat_0_center = 0.5*(NELat_local + SWLat_local)
    lon_0_center = 0.5*(NELng_local + SWLng_local)


    llcrnrlat  = SWLat_local
    llcrnrlon  = SWLng_local
    urcrnrlat  = NELat_local
    urcrnrlon  = NELng_local

    m = Basemap(projection='cyl',llcrnrlat=SWLat_local, urcrnrlat=NELat_local,
            llcrnrlon=SWLng_local, urcrnrlon=NELng_local, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,4),labels=[0,0,0,1], color='k', textcolor='k', linewidth=0.4)
    m.drawparallels(np.arange(-90,90,2),labels=[1,0,0,0],color='k', textcolor='k', linewidth=0.4)

    m.fillcontinents(color='tan', lake_color='aqua')
    m.drawcoastlines(linewidth=0.4, linestyle='solid', color='k',)
    m.drawcountries()
    m.drawstates()
    m.drawrivers()

    settings_params = get_settings()

    region_type = settings_params[0]

    #   -------------------------------------------

    if region_type == 'Circle':
        Polygon_Location = 'None'

        Circle_Location = settings_params[3]
        earthquake_depth = float(settings_params[2])

        if Circle_Location != 'None':
            earthquake_depth = float(settings_params[2])
            CircleCenterLat = float(settings_params[4])
            CircleCenterLng = float(settings_params[5])
            CircleRadius    = float(settings_params[6])

         #   -------------------------------------------

    if region_type == 'Polygon':
        Circle_Location = 'None'

        number_polygon_vertices = (len(settings_params)-4)/2

        Polygon_Location = settings_params[3]
        earthquake_depth = float(settings_params[2])

        if Polygon_Location != 'None':
            earthquake_depth = float(settings_params[2])

        x_poly = []
        y_poly = []

        for j in range(0,number_polygon_vertices):
            y_poly.append(settings_params[4+2*j])
            x_poly.append(settings_params[5+2*j])

    #   Close the polygon

        y_poly.append(settings_params[4])
        x_poly.append(settings_params[5])

        for i in range(0,len(x_poly)):
            x_poly[i] = float(x_poly[i])
            y_poly[i] = float(y_poly[i])

  
   #   -------------------------------------------
   #
   #    Re-download the data in the circular region for the local seismicity plot

    settings_params = get_settings()

    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Radius_float        =   float(settings_params[6])

    mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = read_circle_catalog()

   #    Reverse the magnitude array for the magnitudes of events in the circular region

    mag_array_reversed  =   list(reversed(mag_array))
    date_array_reversed =   list(reversed(date_array))
    time_array_reversed =   list(reversed(time_array))
    year_array_reversed =   list(reversed(year_array))
    depth_array_reversed=   list(reversed(depth_array))
    lat_array_reversed  =   list(reversed(lat_array))
    lng_array_reversed  =   list(reversed(lng_array))

#    print lat_array_reversed[:]
#    print lng_array_reversed[:]
    

   #    Now find the number of small events within the circle since last large event

    todays_count = 0    

    last_eq_mag     = '(No Data)'
    last_eq_date    = '(No Data)'
    last_eq_time    = '(No Data)'
    last_eq_year    = '(No Data)'
    last_eq_lat     = '(No Data)'
    last_eq_lng     = '(No Data)'
    last_eq_depth   = '(No Data)'

   #.................................................................
   #
   #   This next bit finds the data on the last big earthquake

    mag_flag = 0
    for i in range(0,len(mag_array_reversed)):
        if float(mag_array_reversed[i]) >= MagLo and (mag_flag == 0) and (float(depth_array_reversed[i]) <= earthquake_depth):

            last_eq_mag         = str(mag_array_reversed[i])
            last_eq_date        = date_array_reversed[i]
            last_eq_time        = time_array_reversed[i]  
            last_eq_year        = year_array_reversed[i]
            last_eq_lat         = lat_array_reversed[i]
            last_eq_lng         = lng_array_reversed[i]
            last_eq_depth       = depth_array_reversed[i]

            mag_flag = 1
        if float(mag_array_reversed[i]) < MagLo and mag_flag == 0:
            todays_count += 1

   #.................................................................
   #
   #    Write catalog of most recent small earthquakes in the circle

    working_file_local = open("EQ_Working_Catalog_Local.txt", "w")    

    for i in range(0,todays_count):
        working_file_local.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" % \
                (i+1, date_array_reversed[i], time_array_reversed[i], year_array_reversed[i], lng_array_reversed[i], \
                lat_array_reversed[i], mag_array_reversed[i], depth_array_reversed[i]))

    if (last_eq_mag != '(No Data)'):
        working_file_local.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" % (i+2, last_eq_date, last_eq_time, last_eq_year, \
                last_eq_lng, last_eq_lat, last_eq_mag, last_eq_depth))

    working_file_local.close()


   #.................................................................



    #   Open input file and find the number of lines in the file

    input_file = open("EQ_Working_Catalog_Local.txt", "r")

    i=0
    for line in input_file:
        i       +=  1

    input_file.close()  # Put the file pointer back at top

    number_eqs=i

    #   Open input file again

    input_file = open("EQ_Working_Catalog_Local.txt", "r")

    # Loop over lines and extract variables of interest

    # Create arrays of length i filled with zeros

    indx_string         =   ["" for x in range(i)]
    yrs                 =   zeros(i)
    lng                 =   zeros(i)
    lat                 =   zeros(i)
    mag                 =   zeros(i)
    dep                 =   zeros(i)
    time_string         =   ["" for x in range(i)]
    date_string         =   ["" for x in range(i)]

    # Loop over lines and extract variables of interest

    i=-1
    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        i       +=  1

        indx_string[i] =   data_array[0]
        date_string[i] =   data_array[1]
        time_string[i] =   data_array[2]

        yrs[i]     =    float(data_array[3])
        lng[i]     =    float(data_array[4])
        lat[i]     =    float(data_array[5])
        mag[i]     =    float(data_array[6])
        dep[i]     =    float(data_array[7])

    input_file.close()  # Put the file pointer back at top

    z_limit = len(lat)
    range_limit = z_limit - 1

    if z_limit > 0:
        print_text_1 = '     Found ' + str(number_eqs) + ' earthquakes having M>' + str(MagLo) + ' around ' + Location
        print_text_2 = '     Occurring from: ' + date_string[0] + ' ' + time_string[0] 
        print_text_3 = '               thru: ' + date_string[number_eqs-1] + ' ' + time_string[number_eqs - 1]
        print_text_4 = '     at depths < ' + str(settings_params[2]) + ' km'


    #.................................................................

    color_start =   0.0
    color_stop  =   1.0

    #   z_limit is the number of points

    cm_subsection = linspace(color_start, color_stop, z_limit) 
    colors = [ cm.rainbow(x) for x in cm_subsection ]

    mag_diff = 0.1*(float(MagLo) - float(completeness_mag))

    mag1 = completeness_mag
    mag2 = completeness_mag + 1.0*mag_diff
    mag3 = completeness_mag + 2.0*mag_diff
    mag4 = completeness_mag + 3.0*mag_diff
    mag5 = completeness_mag + 4.0*mag_diff
    mag6 = completeness_mag + 5.0*mag_diff
    mag7 = completeness_mag + 6.0*mag_diff
    mag8 = completeness_mag + 7.0*mag_diff
    mag9 = completeness_mag + 8.0*mag_diff
    mag10= completeness_mag + 9.0*mag_diff
    mag11= MagLo

    if z_limit >0:

        for z in range(z_limit):
            x,y = m(lng,lat)

    #   -----------------------

            w = z_limit-z-1

            if mag[w] >= mag1:
                mark_size = 4
            if mag[w] >= mag2:
                mark_size = 5
            if mag[w] >= mag3:
                mark_size = 6
            if mag[w] >= mag4:
                mark_size = 7
            if mag[w] >= mag5:
                mark_size = 8
            if mag[w] >= mag6:
                mark_size = 9
            if mag[w] >= mag7:
                mark_size = 10
            if mag[w] >= mag8:
                mark_size = 11
            if mag[w] >= mag9:
                mark_size = 12
            if mag[w] >= mag10:
                mark_size = 13
            if mag[w] >= mag11:
                mark_size = 20
                
    #   ----------------------
     
            m.plot(x[w], y[w], "o", ms=mark_size, color=colors[z])
            m.plot(x[w], y[w], "o", ms=mark_size, fillstyle='none', color='k', lw=0.4)
            

    if (Circle_Location != 'None'):
        x_circle_dg, y_circle_dg = draw_big_circle(CircleCenterLat, CircleCenterLng, CircleRadius)
        x_circle_dg, y_circle_dg = BSTUtilities.createCircleAroundWithRadius(CircleCenterLat, CircleCenterLng, CircleRadius)
        
        m.plot(x_circle_dg, y_circle_dg, "b-", lw=1.3)

    if (Polygon_Location != 'None'):
        m.plot(x_poly, y_poly, "b-", lw=1.6)

    #   ------------------------------------------------------

    sfv     =   (NELat_local-SWLat_local)/12.0
    sfh     =   (NELng_local-SWLng_local)/16.0

    sfh     =   sfh * 0.9
    sfht    =   sfh * 1.35  # move the text over a bit more

    lat=[llcrnrlat, llcrnrlat + 0.775*sfv]
    lng=[urcrnrlon + 1.25*sfht, urcrnrlon + 1.25*sfht]

    x,y = m(lng,lat)

    latc=[llcrnrlat + 0.15*sfv, llcrnrlat + 0.9*sfv]
    lngc=[urcrnrlon + 0.85*sfh, urcrnrlon + 0.85*sfh]
    
    xc,yc = m(lngc,latc)

    mag_value=[str(completeness_mag),str(MagLo)+'+']

    mark_size = [4, 20]

    for z in range(2):
        plt.text(x[z],y[z],mag_value[z], fontsize=10)
        m.plot(xc[z], yc[z], "ro", ms=mark_size[z], clip_on=False)
        m.plot(xc[z], yc[z], "ro", ms=mark_size[z], clip_on=False, fillstyle='none', color='k', lw=0.6)

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]

    if region_type == 'Circle':

        if Circle_Location == 'None':
            SupTitle_text = 'Earthquake Epicenters for M>' + str(MagLo) + ' near ' + Location

        if Circle_Location != 'None':
            Circle_Location_actual = Circle_Location.split('-')
            SupTitle_text = 'Small Earthquakes M>' + str(completeness_mag) + ' near ' + City_Location 

    if region_type == 'Polygon':

        if Polygon_Location == 'None':
            SupTitle_text = 'Earthquake Epicenters for M>' + str(MagLo) + ' in ' + Location + ' at Depth < ' + str(earthquake_depth) +  'km'

        if Polygon_Location != 'None':
            Polygon_Location_actual = Polygon_Location.split('-')
            SupTitle_text = 'Earthquakes M>' + str(completeness_mag) + ' in ' + Location + ' at Depth < ' + str(earthquake_depth) + ' km'

    plt.suptitle(SupTitle_text, fontsize=14)

    if z_limit > 0:
        Title_text = 'Since M'+ str(mag[range_limit]) + ' on ' + date_string[range_limit] + ' ' + \
                ' To: ' + date_string[0] + ' R<'+ str(int(CircleRadius)) + ' km, D<'+ str(int(earthquake_depth)) +' km'
        plt.title(Title_text, fontsize=12)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')

    figure_name = './Data/' + Circle_Location + '_Seismicity_Local.png'

    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

    #plt.show()

    return None

    ######################################################################        

    
def plot_swarm_event_counts(min_daily_events_factor, burst_min_size_factor, completeness_mag, MagLo, Location, \
        sd_factor, burst_print_flag, regional_rate, cluster_percent):
   
   #   -------------------------------------------
   #
   #    Download the data in the circular region for the local correlated seismicity
   #        and centroid plot

    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])


    #.................................................................
    
    #   Open input file 

    input_file = open("USGS_circle.catalog", "r")

    yrs                 =   []
    lng                 =   []
    lat                 =   []
    mag                 =   []
    dep                 =   []
    time_string         =   []
    date_string         =   []

    # Loop over lines and extract variables of interest

    for line in input_file:
        line = line.strip()
        data = line.split()
        data_array = np.asarray(data)

        date_string.append(data_array[0])
        time_string.append(data_array[1])

        yrs.append(float(data_array[2]))
        lng.append(float(data_array[3]))
        lat.append(float(data_array[4]))
        mag.append(float(data_array[5]))
        dep.append(float(data_array[6]))
        
    input_file.close()  
    
   #.................................................................

    #   Now find the value of min_daily_events to assure a Poisson probability = 99.99% that an active day is a swarm

    n_ev = float(len(yrs))
    n_days = 365.0 * (max(yrs)-min(yrs))
#   mean_daily_rate = float(n_ev)/float(n_days)
    mean_daily_rate = regional_rate
    
    Poisson_prob = 0.0
    n_exp = -1

    while Poisson_prob < cluster_percent: 
    
        n_exp += 1
        Poisson_prob = poisson.cdf(n_exp, mean_daily_rate)
        
    min_daily_events = int(n_exp * min_daily_events_factor)
    burst_min_size = burst_min_size_factor * min_daily_events
    if burst_min_size < 2:
        burst_min_size = 2
    
    Poisson_prob = poisson.cdf(min_daily_events, mean_daily_rate)

    number_of_stddev = (min_daily_events - mean_daily_rate) / math.sqrt(mean_daily_rate)
    
    min_daily_events = 2
    burst_min_size = 2
    
    print('')
    print('----------------------------------------------------------')
    print('')
    print('Minimum daily events required: ', min_daily_events)
    print('')
    print('Mean regional rate of events/day: ', str(round(mean_daily_rate,2)) + ' Events M>' + str(completeness_mag) + ' per Day')
    print('')
#     print 'Poisson probability that more than '+ str(min_daily_events)+ \
#             ' events in a day is a swarm = ', str(round(100.*(Poisson_prob),4)) + '%'
#     print ''
#     print 'Corresponding to: ', str(round(number_of_stddev,1)) + ' Standard Deviations from the Mean'
    print('----------------------------------------------------------')

    print('')
    
    number_events, index_events, burst_index = \
            BSTUtilities.build_daily_counts(burst_min_size, min_daily_events, mag, lat, lng, date_string, time_string, yrs, \
                    sd_factor, burst_print_flag)
    
    #.................................................................
    
    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    number_index = list(range(0,len(number_events)))
    
    plt.plot(number_index, number_events, '.b', ms = 2)
    
    #   Plot intensity threshold
    
    SupTitle_text = 'Event Numbers, $M_c=$' + str(completeness_mag) + ' near ' + City_Location + ' Since ' + date_string[0]
    Title_text = 'R<'+ str(int(Circle_Radius)) + ' km, D<'+ str(int(earthquake_depth)) +' km' 
    
    plt.suptitle(SupTitle_text, fontsize=12)
    plt.title(Title_text, fontsize=10)
    
    plt.ylabel('Event Numbers')
    
    plt.xlabel('Days Since ' + date_string[0])
    
    plt.grid(linestyle = '--', lw=0.8)
    
    figure_name = './Data/' + Circle_Location + '_Number_Events.png'
    
    plt.savefig(figure_name,dpi=600)
    
#     plt.show()
    plt.close('all')

    
#    plt.show()
    
    return burst_index, number_events, min_daily_events, burst_min_size, mag, lat, lng, date_string, time_string, yrs
    
    ######################################################################
    
def map_swarm_epicenters_all(Origin_Location, burst_min_size, completeness_mag, MagLo, Location, N_Steps, \
        burst_index, mag, lat, lng, date, time, years):
                
                
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])

#    mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = read_circle_catalog()

    #.................................................................
    #
    #   Draw basic map
    
    lat_0_center = Circle_Lat
    lon_0_center = Circle_Lng
    
    delta_lat = change_in_latitude(Circle_Radius)
    delta_lat = abs(delta_lat)
    delta_lng = change_in_longitude(Circle_Lat, Circle_Radius)
    delta_lng = abs(delta_lng)
    
    llcrnrlat  = SWLat_local = Circle_Lat - delta_lat - 0.5
    llcrnrlon  = SWLng_local = Circle_Lng - delta_lng - 0.5
    urcrnrlat  = NELat_local = Circle_Lat + delta_lat + 0.5
    urcrnrlon  = NELng_local = Circle_Lng + delta_lng + 0.5

    m = Basemap(projection='cyl',llcrnrlat=SWLat_local, urcrnrlat=NELat_local,
            llcrnrlon=SWLng_local, urcrnrlon=NELng_local, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,1.0),labels=[0,0,0,1], color='k', textcolor='k', linewidth=0.4)
    m.drawparallels(np.arange(-90,90,1.0),labels=[1,0,0,0],color='k', textcolor='k', linewidth=0.4)


    m.fillcontinents(color='tan', lake_color='aqua')
    try:
        m.drawcoastlines(linewidth=0.4, linestyle='solid', color='k')
    except:
        pass
    m.drawcountries()
    m.drawstates()
    m.drawrivers()

    #.................................................................

    color_start =   0.0
    color_stop  =   1.0
    
    z_limit = len(burst_index)

    #   z_limit is the number of bursts/swarms, all points in a burst/swarm will have the same color

    cm_subsection = linspace(color_start, color_stop, z_limit) 
    colors = [ cm.rainbow(x) for x in cm_subsection ]

    if z_limit >0:

        for z in range(z_limit):
        
            if len(burst_index[z]) >= burst_min_size:
        
                index_list = burst_index[z]
                lat_list = []
                lng_list = []
                for i in range(len(burst_index[z])):

                    lat_list.append(lat[index_list[i]])
                    lng_list.append(lng[index_list[i]])
                
                mark_size = 2
            
    #   ----------------------
     
#           m.plot(lng_list, lat_list, ".", ms=mark_size, color=colors[z])

                m.plot(lng_list, lat_list, "ro", ms=mark_size, clip_on=False, color=colors[z])
#               m.plot(lng_list, lat_list, "ro", ms=mark_size, clip_on=False, fillstyle='none', color='k', lw=0.01)
                
#           m.plot(lng_list, lat_list, "ro", ms=mark_size, clip_on=False, fillstyle='none', color=colors[z], lw=0.01)
            
    #.................................................................

    #   Draw outer circle defining circular region
    
    x_circle_dg, y_circle_dg = BSTUtilities.createCircleAroundWithRadius(Circle_Lat, Circle_Lng, Circle_Radius)
        
    m.plot(x_circle_dg, y_circle_dg, "b--", lw=0.9)
    
    origin_lat = [Origin_Location[0]]
    origin_lng = [Origin_Location[1]]
    
    m.plot(origin_lng, origin_lat, 'k*', ms = 6)
    m.plot(origin_lng, origin_lat, 'y*', ms = 4)
        
    #   ------------------------------------------------------

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]

    if Circle_Location != 'None':
        Circle_Location_actual = Circle_Location.split('-')
        SupTitle_text = 'Events M>' + str(completeness_mag) + ' and Swarm Size > ' + str (burst_min_size) + ' Events near ' + City_Location 

    plt.suptitle(SupTitle_text, fontsize=12)
    
    start_date = date[0]
    end_date   = date[len(date)-1]
# 
    if z_limit > 0:
        Title_text = 'Within R<'+ str(int(Circle_Radius)) + ' km, D<'+ str(int(earthquake_depth)) +' km, From ' + \
            start_date + ' to ' + end_date
        plt.title(Title_text, fontsize=10)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')

    figure_name = './Data/' + Circle_Location + '_Swarms_Map.png'

    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

    #plt.show()

    return None

#     .................................................................

def map_swarm_epicenters_one_burst(completeness_mag, MagLo, Location, N_Steps, \
        burst_index, burst_number, mag, lat, lng, date, time, years):
                
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])

    mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = read_circle_catalog()
    
    #   ---------------------
    #

    z = burst_number
    
#     cm_subsection = linspace(color_start, color_stop, len(burst_index[z])) 
#     colors = [ cm.rainbow(x) for x in cm_subsection ]
    
    index_list = burst_index[z]
    lat_list = []
    lng_list = []
    mag_list = []
    
    great_circle_distance_2 =   []
    
    for i in range(len(burst_index[z])):

        lat_list.append(lat[index_list[i]])
        lng_list.append(lng[index_list[i]])
        mag_list.append(mag[index_list[i]])
        
    #.................................................................
        
        burst_centroid_lat = mean_val(lat_list)
        burst_centroid_lng = mean_val(lng_list)
        
    #   Compute Radius of Gyration
    
    for i in range(len(lat_list)):

        lat1 = burst_centroid_lat
        lng1 = burst_centroid_lng
        lat2 = lat_list[i]
        lng2 = lng_list[i]
        dist = BSTUtilities.compute_great_circle_distance(lat1, lng1, lat2, lng2)
            
        great_circle_distance_2.append(math.pow(dist,2))
        
    Radius_Gyration = math.pow(mean_val(great_circle_distance_2), 0.5)

    #.................................................................
    #
    #   Draw basic map
    
    fig = plt.figure(figsize=(8, 6))        #   Define large figure and thermometer - 4 axes needed

    gs = gridspec.GridSpec(1,2,width_ratios=[1,1], wspace = 0.2) 
    ax0 = plt.subplot(gs[0])
    
    lat_0_center = Circle_Lat
    lon_0_center = Circle_Lng
    
#     delta_lat = change_in_latitude(Circle_Radius)
#     delta_lat = abs(delta_lat)
#     delta_lng = change_in_longitude(Circle_Lat, Circle_Radius)
#     delta_lng = abs(delta_lng)
    
    delta_lat = change_in_latitude(1.5*Radius_Gyration)
    delta_lat = abs(delta_lat)
    delta_lng = change_in_longitude(burst_centroid_lat, 1.5*Radius_Gyration)
    delta_lng = abs(delta_lng)
    
    llcrnrlat  = SWLat_local = burst_centroid_lat - delta_lat - 0.25
    llcrnrlon  = SWLng_local = burst_centroid_lng - delta_lng - 0.25
    urcrnrlat  = NELat_local = burst_centroid_lat + delta_lat + 0.25
    urcrnrlon  = NELng_local = burst_centroid_lng + delta_lng + 0.25

    ax0 = Basemap(projection='cyl',llcrnrlat=SWLat_local, urcrnrlat=NELat_local,
            llcrnrlon=SWLng_local, urcrnrlon=NELng_local, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    ax0.drawmeridians(np.arange(0,360,0.5),labels=[0,0,0,1], color='k', textcolor='k', linewidth=0.4, fontsize=8)
    ax0.drawparallels(np.arange(-90,90,0.5),labels=[1,0,0,0],color='k', textcolor='k', linewidth=0.4, fontsize=8)


    ax0.fillcontinents(color='tan', lake_color='aqua')
    try:
        ax0.drawcoastlines(linewidth=0.4, linestyle='solid', color='k')
    except:
        pass
    ax0.drawcountries()
    ax0.drawstates()
    ax0.drawrivers()
    
#     ax0.plot(burst_centroid_lng, burst_centroid_lat, 'k*', ms = 8)
#     ax0.plot(burst_centroid_lng, burst_centroid_lat, '*', ms = 6, mfc='y', mec='k', mew=0.5)
    
    color_start =   0.0
    color_stop  =   1.0
    
    cm_subsection = linspace(color_start, color_stop, len(index_list)) 
    colors = [ cm.rainbow(x) for x in cm_subsection ]
    
    mark_size = 3
    
    for i in range(len(lng_list)):
        if mag_list[i] >= completeness_mag:
            mark_size = 3
        if mag_list[i] > completeness_mag+1:
            mark_size = 5
        if mag_list[i] > completeness_mag+2:
            mark_size = 7
        if mag_list[i] > completeness_mag+3:
            mark_size = 9
        if mag_list[i] > completeness_mag+4:
            mark_size = 11
        if mag_list[i] > completeness_mag+5:
            mark_size = 13
            

        ax0.plot(lng_list[i], lat_list[i], "o", color = colors[i], ms=mark_size, clip_on=False)
        ax0.plot(lng_list[i], lat_list[i], "ko", ms=mark_size, clip_on=False, fillstyle='none', lw=0.01)
        
    #   ------------------------------------------------------

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    first_index     = burst_index[z][0]
    second_index    = len(burst_index[z]) - 1
    second_index    = burst_index[z][second_index]
    
    start_date = date[first_index]
    end_date   = date[second_index]
    
    number_swarm_events = len(mag_list)
    min_mag = round(min(mag_list),1)
    max_mag = round(max(mag_list),1)
    
    centroid_lat = str(round(burst_centroid_lat,2))
    centroid_lng = str(round(burst_centroid_lng,2))

    if Circle_Location != 'None':
#         Circle_Location_actual = Circle_Location.split('-')
#         SupTitle_text = str(number_swarm_events) + ' Swarm Events Centered on ' + City_Location + ' with $M_c\geq$'+ str(completeness_mag)
        SupTitle_text = str(number_swarm_events) + ' Events Centered on ' + \
                '('+ centroid_lat+ '$^o$' + ', ' + centroid_lng + '$^o$' + ')' + ' with $M_c\geq$'+ str(completeness_mag)


    plt.suptitle(SupTitle_text, y=0.95, fontsize=12)
    # 

    Title_text =  'From: ' + start_date + ' To: ' + end_date + ', ' + str(min_mag) +'$\leq M \leq $' + str(max_mag)
    plt.title(Title_text, fontsize=8)

    ax0.drawmapboundary(fill_color='#e6e6ff')
    
#     ax0.plot(burst_centroid_lng, burst_centroid_lat, 'k*', ms = 8)
#     ax0.plot(burst_centroid_lng, burst_centroid_lat, 'y*', ms = 6)
#     
    #   ------------------------------------------------------
    
    ax2 = plt.subplot(gs[1])
    frame1 = plt.gca()
    plt.ylim([int(round(completeness_mag,0)) - 1,8])                         #   Show the y-axis labels
    
    mag_index = list(range(0,len(mag_list)))
    
    for i in range(len(mag_list)):
        x = [mag_index[i], mag_index[i]]
        y = [round(completeness_mag,1), mag_list[i]]
        ax2.plot(x, y, 'b-', lw=0.5)
        
    if len(mag_index) > 0:
        mark_size = 3
    elif len(mag_index) > 100:
        mark_size = 2.5
    elif len(mag_list) > 200:
        mark_size = 2
    elif len(mag_list) > 300:
        mark_size = 1.5
    elif len(mag_list) > 400:
        mark_size = 1.0
    
    ax2.plot(mag_index, mag_list, 'ro', ms = mark_size)
    
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    
    zero_line_x =   [xmin,xmax]
    zero_line_y =   [completeness_mag, completeness_mag]
    
    ax2.plot(zero_line_x, zero_line_y, 'k-', lw=0.5)
        
    ax2.grid(True, ls = '--', lw = 0.5)
    
    Title_text = 'Magnitude vs. Event Number'
    plt.title(Title_text, fontsize=8)
    
    plt.ylabel('Magnitude', fontsize = 8)
    plt.xlabel('Event Number', fontsize = 8)
    
    #   ------------------------------------------------------
    #   
    text_x      = (xmax-xmin)*0.165 + xmin
    text_y      = (ymax-ymin) * 0.95 + ymin
    text_string = 'Radius of Gyration: ' + str(round(Radius_Gyration,2)) + ' Km'
    plt.text(text_x,text_y,text_string, fontsize=9)    
    
    #   ------------------------------------------------------

    figure_name = './Data/SwarmData/' + Circle_Location + '_Swarm_Map_Number_' + str(z) + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

#    plt.show()

    return None
    
    ######################################################################
    
    
def log_freq_bursts(burst_index, burst_min_size, completeness_mag, Location):

    #   This method plots the GR relation for the bursts, i.e., the survivor function
    #       for number vs log_10(burst mass)
    
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])
    
    mass_list = []
    mag_list  = []
        
    for i in range(len(burst_index)):
        mass_list.append(int(len(burst_index[i])))
        
    for i in range(len(mass_list)):
        mag_list.append(math.log(1.5*float(mass_list[i]),10))
        
    min_burst_mass =  min(mass_list)
    max_burst_mass =  max(mass_list)
    
    min_mag_burst =  min(mag_list)
    max_mag_burst =  max(mag_list) 
    
    number_mag_bins = int ((max_mag_burst - min_mag_burst + 1) * 10)  #   0.1 magnitude bins
    
    freq_mag_bins_sdf = [0.0 for i in range(number_mag_bins)]
    
    bins_list       = [int(min_mag_burst) + float(i)*0.1 for i in range(number_mag_bins)]
    bins_list       = [float(int(i*10))*0.1 for i in range(len(bins_list))]
    
    bin_diff = 0.1
    min_mag = int(min_mag_burst)
    
    for i in range(len(mag_list)):
        bin_number      = int((mag_list[i]-min_mag)/bin_diff + 1 )
        for j in range(bin_number):
            freq_mag_bins_sdf[j] += 1
            
    City = Location.split('_')
    Number_Spaces = len(City)
    
    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]    #   Is a list of strings, length 1
#     City_Location_actual = City_Location.split('-')
    
    SupTitle_text = 'Cumulative Frequency-Magnitude for Swarms in ' + Location_Actual + ' Region'

    Title_text = 'Events M>' + str(completeness_mag) + ' and Swarm Size > ' + str (burst_min_size)
                
    plt.plot(bins_list, freq_mag_bins_sdf, 'k--', lw=0.8)
    plt.plot(bins_list, freq_mag_bins_sdf, 'bo', ms = 3)
    
    plt.suptitle(SupTitle_text, fontsize=12)
    plt.title(Title_text, fontsize=8)
    
    plt.ylabel('Cumulative Number $\geq$ M', fontsize = 8)
    plt.xlabel('Swarm Magnitude (Mass) M', fontsize = 8)

    figure_name = './Data/' + Circle_Location + '_FreqMag_Bursts' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')
    
#    number_bursts   =   

    return
    
        ######################################################################
        
def plot_swarm_migration_radgyr(burst_min_size, completeness_mag, Location, Origin_Location, \
        circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years):
        
    #   This routine plots the space-time migration of the swarms, symbol size denotes radius of gyration
    
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])
    
    mean_lat        =   []
    mean_lng        =   []
    mean_yrs        =   []
    swarm_distance  =   []
    year_list       =   []
    
    origin_lat = Origin_Location[0]
    origin_lng = Origin_Location[1]

    for i in range(len(burst_index)):        #   Iterating over the bursts/swarms
    
        if len(burst_index[i]) >= burst_min_size:
        
            index_list = burst_index[i]
            lat_list = []
            lng_list = []
            yrs_list = []
                
            for j in range(len(burst_index[i])):
            
                lat_list.append(lat[index_list[j]])
                lng_list.append(lng[index_list[j]])
                yrs_list.append(years[index_list[j]])
                
            mean_lat = np.mean(lat_list)
            mean_lng = np.mean(lng_list)
            mean_yrs = np.mean(yrs_list)
            
            great_circle_distance = BSTUtilities.compute_great_circle_distance(origin_lat, origin_lng, mean_lat, mean_lng)
            
            radius_gyration, burst_centroid_lat, burst_centroid_lng = \
                    BSTUtilities.compute_burst_radius_gyration(burst_index[i], mag, lat, lng, date, time, years)
            
            size_ratio = float(radius_gyration)
            mark_size = 3 * math.log(10.*size_ratio,10)
            mark_size = 3.*math.pow(size_ratio,.3)
            if mark_size <= 3:
                mark_size = 3
            
            plt.plot([mean_yrs],[great_circle_distance],'ro', ms = mark_size)
            
    plt.grid(True, ls = '--', lw = 0.5)
    
    xmin,xmax = plt.xlim()
    xmin,xmax = plt.xlim(2010.0, 2020.0)
    
    xmin,xmax = plt.xlim(circle_catalog_date_start - 1., circle_catalog_date_end + 1.)
    
    
    ymin, ymax = plt.ylim()
    ymin, ymax = plt.ylim(-5,ymax)
    
    x_EMC = [2010.255, 2010.255]
    y_EMC = [0,ymax-5]
    
    plt.plot(x_EMC,y_EMC,'b--', lw=0.8)
    
    SupTitle_text = 'Migration of Swarm Events Relative to Lat, Lng: ' + str(origin_lat) + '$^o$, ' + str(origin_lng) + '$^o$'

    plt.suptitle(SupTitle_text, fontsize=12)
    
    Title_text = 'Distance vs. Time for Minimum Swarm Size of ' + str(burst_min_size) + ' Events'
    plt.title(Title_text, fontsize=8)
    
    plt.ylabel('Distance (Km)', fontsize = 8)
    plt.xlabel('Time (Years)', fontsize = 8)
    
    text_x      = (xmax-xmin)*0.10 + xmin
    text_y      = (ymax-ymin) * 0.93 + ymin
    text_string = 'Circle Size Proportional to Swarm Radius of Gyration'
    plt.text(text_x,text_y,text_string, fontsize=8)    
    
    figure_name = './Data/' + Circle_Location + '_Swarm_Migration_RGyr' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

    return
    
    ######################################################################
    
def plot_swarm_migration_mass(burst_min_size, completeness_mag, Location, Origin_Location, \
        circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years):
        
    #   This routine plots the space-time migration of the swarms
    
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])
    
    mean_lat        =   []
    mean_lng        =   []
    mean_yrs        =   []
    swarm_distance  =   []
    year_list       =   []
    
    origin_lat = Origin_Location[0]
    origin_lng = Origin_Location[1]

    for i in range(len(burst_index)):        #   Iterating over the bursts/swarms
    
        if len(burst_index[i]) >= burst_min_size:
        
            index_list = burst_index[i]
            lat_list = []
            lng_list = []
            yrs_list = []
                
            for j in range(len(burst_index[i])):
            
                lat_list.append(lat[index_list[j]])
                lng_list.append(lng[index_list[j]])
                yrs_list.append(years[index_list[j]])
                
            mean_lat = np.mean(lat_list)
            mean_lng = np.mean(lng_list)
            mean_yrs = np.mean(yrs_list)
            
            great_circle_distance = BSTUtilities.compute_great_circle_distance(origin_lat, origin_lng, mean_lat, mean_lng)
            
            size_ratio = float(len(lat_list))/float(burst_min_size)
            mark_size = 3 * math.log(10.*size_ratio,10)
            mark_size = 3.*math.pow(size_ratio,.5)
            
            plt.plot([mean_yrs],[great_circle_distance],'ro', ms = mark_size)
            
    plt.grid(True, ls = '--', lw = 0.5)
    
    xmin,xmax = plt.xlim()
    xmin,xmax = plt.xlim(circle_catalog_date_start - 1., circle_catalog_date_end + 1.)
    
    ymin, ymax = plt.ylim()
    ymin, ymax = plt.ylim(-5,ymax)
    
    x_EMC = [2010.255, 2010.255]
    y_EMC = [0,ymax-5]
    
    plt.plot(x_EMC,y_EMC,'b--', lw=0.8)
    
    SupTitle_text = 'Migration of Swarm Events Relative to Lat, Lng: ' + str(origin_lat) + '$^o$, ' + str(origin_lng) + '$^o$'

    plt.suptitle(SupTitle_text, fontsize=12)
    
    Title_text = 'Distance vs. Time for Minimum Swarm Size of ' + str(burst_min_size) + ' Events'
    plt.title(Title_text, fontsize=8)
    
    plt.ylabel('Distance (Km)', fontsize = 8)
    plt.xlabel('Time (Years)', fontsize = 8)
    
    text_x      = (xmax-xmin)*0.10 + xmin
    text_y      = (ymax-ymin) * 0.93 + ymin
    text_string = 'Circle Size Proportional to Swarm Mass'
    plt.text(text_x,text_y,text_string, fontsize=8)    
    
    figure_name = './Data/' + Circle_Location + '_Swarm_Migration_Mass' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

    return
    
    ######################################################################
    
def plot_radgyr_time(burst_min_size, completeness_mag, Location, Origin_Location,\
        circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years):
        
    #   This routine plots the space-time migration of the swarms, symbol size denotes radius of gyration
    
    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])
    
    mean_lat        =   []
    mean_lng        =   []
    mean_yrs        =   []
    swarm_distance  =   []
    year_list       =   []
    radgyr_list     =   []
    year_swarm   =   []
    
    origin_lat = Origin_Location[0]
    origin_lng = Origin_Location[1]

    for i in range(len(burst_index)):        #   Iterating over the bursts/swarms
    
        if len(burst_index[i]) >= burst_min_size:
        
            index_list = burst_index[i]
            lat_list = []
            lng_list = []
            yrs_list = []
                
            for j in range(len(burst_index[i])):
            
                lat_list.append(lat[index_list[j]])
                lng_list.append(lng[index_list[j]])
                yrs_list.append(years[index_list[j]])
                
            year_swarm.append(years[index_list[0]])
            radius_gyration, burst_centroid_lat, burst_centroid_lng = \
                    BSTUtilities.compute_burst_radius_gyration(burst_index[i], mag, lat, lng, date, time, years)
            
            radgyr_list.append(radius_gyration)     #   This is the time series we want to use EMA on
            
    radgyr_index = list(range(len(radgyr_list)))

    mark_size = 3

    plt.plot(year_swarm,radgyr_list,'--', lw=0.8)
    plt.plot(year_swarm,radgyr_list,'bo', ms = mark_size)
            
    plt.grid(True, ls = '--', lw = 0.5)
        
    SupTitle_text = 'Radius of Gyration for Swarms vs. Time'

    plt.suptitle(SupTitle_text, fontsize=12)
    
    Title_text = 'Radius (Km) vs. Time for Minimum Swarm Size of ' + str(burst_min_size) + ' Events near ' + City_Location
    plt.title(Title_text, fontsize=8)
    
    plt.ylabel('Radius of Gyration (Km)', fontsize = 8)
    plt.xlabel('Time (Year)', fontsize = 8)
    
#     text_x      = (xmax-xmin)*0.10 + xmin
#     text_y      = (ymax-ymin) * 0.93 + ymin
#     text_string = 'Circle Size Proportional to Swarm Radius of Gyration'
#     plt.text(text_x,text_y,text_string, fontsize=8)    
    
    figure_name = './Data/' + Circle_Location + '_RadGyr_Time' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

    return
    
    ######################################################################
    
def plot_radgyr_EMA_time(burst_min_size, completeness_mag, Location, Origin_Location, N_Steps,\
        circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
        ratio_limit, sd_factor, cutoff_start_year, regional_rate):
        
    #   This routine plots the space-time migration of the swarms, symbol size denotes radius of gyration
    
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])
    
    mean_lat        =   []
    mean_lng        =   []
    mean_yrs        =   []
    swarm_distance  =   []
    year_list       =   []
    radgyr_list     =   []
    radgyr_EMA_list =   []
    year_swarm      =   []
    
    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    for i in range(3,len(burst_index)):        #   Iterating over the bursts/swarms
    
        if len(burst_index[i]) >= burst_min_size:
        
            index_list = burst_index[i]
            lat_list = []
            lng_list = []
            yrs_list = []
                
            for j in range(len(burst_index[i])):
            
                lat_list.append(lat[index_list[j]])
                lng_list.append(lng[index_list[j]])
                yrs_list.append(years[index_list[j]])
                
            radius_gyration, burst_centroid_lat, burst_centroid_lng =\
                     BSTUtilities.compute_burst_radius_gyration(burst_index[i], mag, lat, lng, date, time, years)
            
            radgyr_list.append(radius_gyration)     #   This is the time series we want to use EMA on
            year_swarm.append(years[index_list[0]])

#   Filter the clusters now, we filtered the events in the clusters previously
    
#     radgyr_list, year_swarm, = BSTUtilities.reject_outliers_clusters(radgyr_list, year_swarm)

    radgyr_list, year_swarm, index_swarm, centroid_lat, centroid_lng = \
            BSTUtilities.reject_outliers_clusters_alternate(ratio_limit, burst_index, burst_min_size, mag, lat, lng, date, time, years)
    
    average_rate_start_date =   circle_catalog_date_start
    average_rate_end_date   =   circle_catalog_date_end
    
    if average_rate_end_date < 1990.0:
        average_rate_end_date = 1990.0      #   The beginning of decent digital data
    
    year_interval = average_rate_end_date - average_rate_start_date
    
    number_earthquakes = 0
    number_bursts_in_plot = 0
    
    for i in range(len(burst_index)):
        if years[burst_index[i][0]] >= average_rate_start_date:
            number_earthquakes += 1
            number_bursts_in_plot += 1
        
    for i in range(1,len(radgyr_list)+1):
        RG_list_raw = []
        
        for j in range(i):
            RG_list_raw.append(radgyr_list[j])
            
#         RG_list_raw = [math.log(i,10) for i in RG_list_raw]
        
        radgyr_EMA = BSTUtilities.EMA_weighted_time_series(RG_list_raw, N_Steps)
        
        radgyr_EMA_list.append(radgyr_EMA)
        
#     radgyr_EMA_list = [10**i for i in radgyr_EMA_list]    #   Convert from Log back to km
        
    #   Determine which bursts occurred after 1990.0
    
#     kk = 0
    
    year_swarm_reduced  =   []
    radgyr_EMA_reduced  =   []
    
    for i in range(len(year_swarm)):
        if year_swarm[i] > cutoff_start_year:
            year_swarm_reduced.append(year_swarm[i])
            radgyr_EMA_reduced.append(radgyr_EMA_list[i])
            
    year_swarm = year_swarm_reduced
    radgyr_EMA_list = radgyr_EMA_reduced
    
    earthquake_year     =   []
    earthquake_mags     =   []
    earthquake_date     =   []
    earthquake_time     =   []
    earthquake_lat      =   []
    earthquake_lng      =   []
    
    for i in range(len(years)):
        if float(mag[i]) >= 6.0:
            earthquake_year.append(years[i])
            earthquake_mags.append(mag[i])
            earthquake_date.append(date[i])
            earthquake_time.append(time[i])
            earthquake_lat.append(lat[i])
            earthquake_lng.append(lng[i])
            
    mark_size = 2
    
    EQ_Table = []

    for i in range(0,len(earthquake_date)):
        EQ_Table.append([i+1, earthquake_date[i], \
                earthquake_year[i],earthquake_time[i], earthquake_mags[i], earthquake_lat[i], earthquake_lng[i]])
    
    table_data = tabulate(EQ_Table,headers=['#', 'Date', 'Decimal Year', 'Time(Z)', 'Mag', 'Lat', 'Lng'], numalign="center")
    print(table_data)
    
#     for i in range(len(earthquake_year)):
#         print '#, Date, Decimal Year, Time(Z), Magnitude, Lat, Lng: ', i, earthquake_date[i], \
#                 earthquake_year[i],earthquake_time[i], earthquake_mags[i], earthquake_lat[i], earthquake_lng[i]
    
    print('')
    print('----------------------------------------------------------')
    print('')
    print('ENF, CLF, len(year_swarm), len(radgyr_EMA_list): ', ratio_limit, sd_factor, len(year_swarm), len(radgyr_EMA_list))
    

# 
# ---------------------------------------

                
    plt.plot(year_swarm,radgyr_EMA_list, linestyle='--', lw=0.5, color='b', zorder=3)
    plt.plot(year_swarm,radgyr_EMA_list,'b.', ms = 4, zorder=3)
    
    plt.gca().invert_yaxis()
            
    plt.grid(True, lw = 0.5, linestyle='dotted', zorder=0, axis = 'y')
    
    xmin,xmax = plt.xlim()
#     xmin,xmax = plt.xlim(1990,xmax)
    ymin, ymax = plt.ylim()
    
    min_plot_line = [ymin for i in range(len(year_swarm))]
    plt.fill_between(year_swarm,min_plot_line, radgyr_EMA_list, color='c', alpha=0.1, zorder=0)
    
    for i in range(len(earthquake_year)):
        x_eq = [earthquake_year[i], earthquake_year[i]]
        y_eq = [ymin,ymax]
        
#         if float(earthquake_mags[i]) < 6.5 and float(earthquake_mags[i]) >= 6.0 and float(earthquake_year[i]) >= 1980 :
#             plt.plot(x_eq, y_eq, '--', color='c', lw=0.8, zorder=0)
#                         
        if float(earthquake_mags[i]) >= 6.0 and float(earthquake_mags[i]) < 7.0  and float(earthquake_year[i]) >= cutoff_start_year:
            plt.plot(x_eq, y_eq, linestyle='dotted', color='k', lw=0.7, zorder=1)
            
        if float(earthquake_mags[i]) >= 7.0 and float(earthquake_year[i]) >= cutoff_start_year:
            plt.plot(x_eq, y_eq, 'r', linestyle='--', lw=0.7, zorder=2)
            
    xmin,xmax = plt.xlim()
    ymin, ymax = plt.ylim()
#     plt.ylim(30,ymax)
#     
    plt.minorticks_on()
    
    SupTitle_text = 'EMA Radius of Gyration for Swarms M $\geq$' + str(completeness_mag) + ' vs. Time Within ' + str(Circle_Radius) + ' Km of ' +  City_Location

    plt.suptitle(SupTitle_text, fontsize=10, y = 0.96)
    
    if ratio_limit > -100.0:
        Title_text = 'N $\geq$ ' + str(burst_min_size) + ' Events;  Depth $\leq$ '\
                +   str(earthquake_depth) + ' Km; N_Steps = ' + str(N_Steps) \
                + '; ENF = ' + str(ratio_limit) + '; CLF = ' + str(sd_factor)
    else:
        Title_text = 'Swarms > ' + str(burst_min_size) + ' Events Within ' + str(Circle_Radius) + ' Km, at Depth < '\
                +   str(earthquake_depth) + ' Km, with N_Steps = ' + str(N_Steps) \
                + ' and sd_factor = ' + str(sd_factor)
            
    plt.title(Title_text, fontsize=9)
    
    plt.ylabel('Radius of Gyration (Km), Inverted Scale', fontsize = 12)
    plt.xlabel('Time (Year)', fontsize = 12)
    
#     text_x      = (xmax-xmin)*0.10 + xmin
#     text_y      = (ymax-ymin) * 0.93 + ymin
#     text_string = 'Circle Size Proportional to Swarm Radius of Gyration'
#     plt.text(text_x,text_y,text_string, fontsize=8)    
    
    figure_name = './Data/' + Circle_Location + '_RadGyr_Time_EMA' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

    return N_Steps, radgyr_list, year_swarm, index_swarm, centroid_lat, centroid_lng 
    
    ######################################################################
    
def map_local_seismicity(eq_start, eq_end, Location, Local_Circle_Radius, completeness_mag,\
        mag, lat, lng, date, time, years):

    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius       =   Local_Circle_Radius
    
    #   ---------------------
    #
    
    #   Filter the events
    
    lat_local   =   []
    lng_local   =   []
    
    lng_circle_dg, lat_circle_dg = BSTUtilities.createCircleAroundWithRadius(Circle_Lat, Circle_Lng, Local_Circle_Radius)

    number_polygon_vertices = len(lng_circle_dg)
    
    point_list = []

    for i in range(number_polygon_vertices):
        point_list.append((float(lat_circle_dg[i]),float(lng_circle_dg[i])))
    
    polygon = Polygon(point_list)
    
    delta_days = 0.0030 #   1 day = 0.00274 years
    
    for i in range(len(years)):
        
        point = Point((float(lat[i]),float(lng[i])))
        
        if (polygon.contains(point) == True and years[i] > (eq_start+delta_days) and years[i] < (eq_end - delta_days)\
                 and float(mag[i]) >= completeness_mag):
            lat_local.append(lat[i])
            lng_local.append(lng[i])

#  Earthquake Dates:  [Landers, Hector Mine, El Major Cucupah, Ridgecrest]
    earthquake_dates = [1992.491803, 1994.046448, 1999.789617, 2010.256831, 2019.510929]
    
#     color_start =   0.0
#     color_stop  =   1.0
#     
#     cm_subsection = linspace(color_start, color_stop, len(earthquake_dates)) 
#     colors = [ cm.rainbow(x) for x in cm_subsection]
#         

#   z_limit is the number of bursts/swarms, all points in a burst/swarm will have the same color

    color_start = 0.0
    color_stop  = 1.0

    cm_subsection = linspace(color_start, color_stop, len(lat_local)) 
    colors = [ cm.rainbow(x) for x in cm_subsection ]
    
    #.................................................................
    
    #
    #   Draw map of burst centroids
    
    lat_0_center = Circle_Lat
    lon_0_center = Circle_Lng
    
    delta_lat = change_in_latitude(Local_Circle_Radius)
    delta_lat = abs(delta_lat)
    delta_lng = change_in_longitude(Circle_Lat, Local_Circle_Radius)
    delta_lng = abs(delta_lng)
    
    llcrnrlat  = SWLat_local = Circle_Lat - delta_lat - 0.5
    llcrnrlon  = SWLng_local = Circle_Lng - delta_lng - 0.5
    urcrnrlat  = NELat_local = Circle_Lat + delta_lat + 0.5
    urcrnrlon  = NELng_local = Circle_Lng + delta_lng + 0.5

    m = Basemap(projection='cyl',llcrnrlat=SWLat_local, urcrnrlat=NELat_local,
            llcrnrlon=SWLng_local, urcrnrlon=NELng_local, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,1.0),labels=[0,0,0,1], color='k', textcolor='k', linewidth=0.4)
    m.drawparallels(np.arange(-90,90,1.0),labels=[1,0,0,0],color='k', textcolor='k', linewidth=0.4)

    m.fillcontinents(color='tan', lake_color='aqua')
    m.drawcoastlines(linewidth=0.4, linestyle='solid', color='k')
    m.drawcountries()
    m.drawstates()
    m.drawrivers()
    
    mark_size = 2
    
    for i in range(len(lat_local)):
        x_lng = [lng_local[i]]
        y_lat = [lat_local[i]]
        m.plot(x_lng, y_lat, '.', color=colors[i], ms=mark_size)
        
#    m.plot(lng_local, lat_local, "r.", ms=mark_size)

    x_circle_dg, y_circle_dg = BSTUtilities.createCircleAroundWithRadius(Circle_Lat, Circle_Lng, Local_Circle_Radius)
    m.plot(x_circle_dg, y_circle_dg, "b--", lw=0.9)
            
    #   ------------------------------------------------------

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    if Circle_Location != 'None':
        Circle_Location_actual = Circle_Location.split('-')
        SupTitle_text = 'Local Seismicity within ' + str(Local_Circle_Radius) + ' Km of ' + City_Location 

    plt.suptitle(SupTitle_text, y=0.97, fontsize=12)
    # 
    
    Title_text =  'From: ' + str(eq_start+0.2) + ' To: ' + str(eq_end-0.2)
    plt.title(Title_text, fontsize=10)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')
    

    figure_name = './Data/' + Circle_Location + '_Local_PreSeismic_Activity' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

#    plt.show()

    return None
    
    ######################################################################
    
def plot_cluster_density_rgyr(burst_min_size, completeness_mag, Location, Origin_Location, N_Steps,\
        circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
        ratio_limit, cutoff_start_year):
        
    #   This plots the density against the radius of gyration, to allow identification of families of burst clusters
    
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])
    
    mean_lat        =   []
    mean_lng        =   []
    mean_yrs        =   []
    swarm_distance  =   []
    radgyr_list     =   []
    radgyr_EMA_list =   []
    mass_list       =   []
    
    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    mass_list   = []
    
    for i in range(3,len(burst_index)):        #   Iterating over the bursts/swarms
    
        if len(burst_index[i]) >= burst_min_size:
        
            index_list = burst_index[i]
            lat_list    = []
            lng_list    = []
            mass_list.append(float(len(index_list)))

                
            for j in range(len(burst_index[i]) and years[index_list[0]] >= cutoff_start_year):
            
                lat_list.append(lat[index_list[j]])
                lng_list.append(lng[index_list[j]])

                
            radius_gyration, burst_centroid_lat, burst_centroid_lng = \
                    BSTUtilities.compute_burst_radius_gyration(burst_index[i], mag, lat, lng, date, time, years)
            
            radgyr_list.append(radius_gyration)     #   This is the time series we want to use EMA on
            
    density_list   =   []
    for i in range(len(mass_list)):
        density = mass_list[i]/((radgyr_list[i]))
        density_list.append(density)
        
    mean_density_list = BSTUtilities.mean_val(density_list)
    
    median_density_list = median(density_list)
    
    Log_10_median_density_list = math.log(median_density_list,10)
    
    Log_10_mean_density_list = math.log(mean_density_list,10)
        
    density_list = [math.log(density_list[i],10) for i in range(len(density_list))]
    
    mean_log_10_density = BSTUtilities.mean_val(density_list)
    mean_log_10_density = round(mean_log_10_density,2)
    
    std_dev_log_10_density, variance_dev_log_10_density = BSTUtilities.std_var_val(density_list)
    std_dev_log_10_density = round(std_dev_log_10_density,2)
    
    print('')
    print('----------------------------------------------------------')
    print('')
    print('Mean(Density): ', mean_density_list)
    print('')
    print('Log_10_mean_density_list: ', Log_10_mean_density_list)
    print('')
#    print 'Mean Log_10(Density):', round(sum(density_list)/float(len(density_list)),2)
    print('Mean Log_10(Density):', mean_log_10_density, '+/-', std_dev_log_10_density)
    print('')
    print('Median Log_10(Density):', median(density_list))
    print('')
    print('Log_10_median_density_list: ', Log_10_median_density_list)
    print('')
    print('----------------------------------------------------------')
    print('')
        
#     print 'len(density_list), len(radgyr_list): ', len(density_list), len(radgyr_list)

    plt.plot(radgyr_list, density_list,'b.', ms = 4, zorder=1)
            
    plt.grid(True, lw = 0.5, linestyle='dotted', zorder=0, axis = 'both')
    
    
    xmin,xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    
    x_ratio = [xmin,xmax]
    y_ratio = [ratio_limit,ratio_limit]
    plt.plot(x_ratio,y_ratio, '--',color='brown', lw=1.0)
#     
    plt.minorticks_on()
    
    SupTitle_text = 'Burst Cluster Density vs. Radius of Gyration for M $\geq$' + str(completeness_mag) + ' near ' +  City_Location

    plt.suptitle(SupTitle_text, fontsize=10)
    
    Title_text = 'Swarms Size > ' + str(burst_min_size) + ' Events Within ' + str(Circle_Radius) + ' Km, at Depth < '\
            +str(earthquake_depth) + ' Km'
            
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

    plt.title(Title_text, fontsize=8)
    
    plt.ylabel('$Log_{10}$ (Mass Ratio, N/Km)', fontsize = 8)
    plt.xlabel('Radius of Gyration (Km)', fontsize = 8)
    
#     text_x      = (xmax-xmin)*0.10 + xmin
#     text_y      = (ymax-ymin) * 0.93 + ymin
#     text_string = 'Circle Size Proportional to Swarm Radius of Gyration'
#     plt.text(text_x,text_y,text_string, fontsize=8)    
    
    figure_name = './Data/' + Circle_Location + '_Density_RadGyr' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

    return
    
    ######################################################################
    
    ######################################################################
    
def calc_radgyr_EMA_time_dev(burst_min_size, completeness_mag, Location, Origin_Location, N_Steps,\
        circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
        ratio_limit, sd_factor, cutoff_start_year, regional_rate):
        
    #   This routine calculates a single radius of gyration time series for given EVF and CLF values
    
    mean_lat        =   []
    mean_lng        =   []
    mean_yrs        =   []
    swarm_distance  =   []
    year_list       =   []
    radgyr_list     =   []
    radgyr_EMA_list =   []
    year_swarm      =   []
    index_swarm     =   []
    

    for i in range(3,len(burst_index)):        #   Iterating over the bursts/swarms
    
        if len(burst_index[i]) >= burst_min_size:
        
            index_list = burst_index[i]
            lat_list = []
            lng_list = []
            yrs_list = []
                
            for j in range(len(burst_index[i])):
            
                lat_list.append(lat[index_list[j]])
                lng_list.append(lng[index_list[j]])
                yrs_list.append(years[index_list[j]])
                
            radius_gyration, burst_centroid_lat, burst_centroid_lng = \
                    BSTUtilities.compute_burst_radius_gyration(burst_index[i], mag, lat, lng, date, time, years)
            
            radgyr_list.append(radius_gyration)     #   This is the time series we want to use EMA on
            year_swarm.append(years[index_list[0]])
            index_swarm.append(index_list[0])

#   Filter the clusters now, we filtered the events in the clusters previously
    
#     radgyr_list, year_swarm, = BSTUtilities.reject_outliers_clusters(radgyr_list, year_swarm)

    radgyr_list, year_swarm, index_swarm, centroid_lat, centroid_lng = \
            BSTUtilities.reject_outliers_clusters_alternate(ratio_limit, burst_index, burst_min_size, mag, lat, lng, date, time, years)

    #   Independent of any start date down to here -------------------------
            
    average_rate_start_date =   circle_catalog_date_start
    average_rate_end_date   =   circle_catalog_date_end
    
    if average_rate_end_date < 1990.0:
        average_rate_end_date = 1990.0      #   The beginning of decent digital data
    
    year_interval = average_rate_end_date - average_rate_start_date
    
    number_earthquakes = 0
    number_bursts_in_plot = 0
    
    for i in range(len(burst_index)):
        if years[burst_index[i][0]] >= average_rate_start_date:
            number_earthquakes += 1
            number_bursts_in_plot += 1
        
    for i in range(1,len(radgyr_list)+1):
        RG_list_raw = []
        
        for j in range(i):
            RG_list_raw.append(radgyr_list[j])
            
        radgyr_EMA = BSTUtilities.EMA_weighted_time_series(RG_list_raw, N_Steps)
        
        radgyr_EMA_list.append(radgyr_EMA)
        
    #   Determine which bursts occurred after 1990.0
    
#     kk = 0
    
    year_swarm_reduced  =   []
    radgyr_EMA_reduced  =   []
    index_swarm_reduced  =   []
    
    for i in range(len(year_swarm)):
        if year_swarm[i] > cutoff_start_year:
            year_swarm_reduced.append(year_swarm[i])
            radgyr_EMA_reduced.append(radgyr_EMA_list[i])
            index_swarm_reduced.append(index_swarm[i])
            
    year_swarm = year_swarm_reduced
    radgyr_EMA_list = radgyr_EMA_reduced
    
    print('')
    print('----------------------------------------------------------')
    print('')
    print('ENF, CLF, len(year_swarm), len(radgyr_EMA_list): ', ratio_limit, sd_factor, len(year_swarm), len(radgyr_EMA_list))
    
    
#   -----------------------------------------------------------------
#   Here we build a continuous, constant length, regular interval, time series for radgyr as a function of time
#
#     radgyr_EMA_list = [round(i,2) for i in radgyr_EMA_list]
#     print 'radgyr_EMA_list: ', radgyr_EMA_list
# #     
        
    for i in range(len(years)):
        if years[i] <= cutoff_start_year:       #   Find the index corresponding to cutoff_year
            cutoff_index = i
        else:
            pass
    
    length_lists = len(years) - cutoff_index
    
    year_swarm_regular  =   []
    radgyr_list_regular =   []
    
    for i in range(length_lists):
#         print 'len(years), i+length_lists: ', len(years), i+cutoff_index
        year_swarm_regular.append(years[i+cutoff_index])
        radgyr_list_regular.append(0.0)
        
    for i in range(len(radgyr_EMA_list)):
        radgyr_list_regular[index_swarm_reduced[i]-cutoff_index] = radgyr_EMA_list[i]
    
    radgyr_value = radgyr_EMA_list[0]
    
    for i in range(len(radgyr_list_regular)):

        if radgyr_list_regular[i] == 0.0:
            radgyr_list_regular[i] += round(radgyr_value,2)
        
        if radgyr_list_regular[i] > 0.0:
            radgyr_value = round(radgyr_list_regular[i],2)
            

    return year_swarm_regular,radgyr_list_regular
    
    ######################################################################
    
    ######################################################################
    
def plot_radgyr_EMA_time_dev(burst_min_size, completeness_mag, Location, Origin_Location, N_Steps,\
        circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
        ratio_limit, sd_factor, cutoff_start_year, regional_rate):
        
    #   This routine plots the space-time migration of the swarms, symbol size denotes radius of gyration
    
#   -----------------------------------------------------------------

    

#     ENF_list = [-1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2]
#     
#     CLF_list = [25, 50, 75, 100, 125]
#     
    ENF_list = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2]
#     
#     CLF_list = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

#     ENF_list = [-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0]
#     
    CLF_list = [5, 10, 15, 20, 25]
#     
#     ENF_list = [ -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
#     
#     CLF_list = [50, 75, 100, 125, 150]

#     ENF_list = [ -0.025, 0.025]
#     
#     CLF_list = [24, 26]

    for n in range(len(ENF_list)):
    
        ratio_limit = ENF_list[n]
        
        for k in range(len(CLF_list)):
        
            sd_factor   = CLF_list[k]
        
            print('')
            print(' Now working on ENF value: ', ratio_limit)
            print(' Now working on CLF value: ', sd_factor)
        
            year_swarm_temp,radgyr_swarm_temp = \
                    calc_radgyr_EMA_time_dev(burst_min_size, completeness_mag, Location, Origin_Location, N_Steps,\
                    circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
                    ratio_limit, sd_factor, cutoff_start_year, regional_rate)
                
            if n == 0 and k == 0:
                year_swarm_regular  = [i for i in year_swarm_temp]
                radgyr_swarm_regular = [[] for i in range(len(year_swarm_temp))]
                radgyr_swarm_regular_mean = [0.0 for i in range(len(year_swarm_temp))]
                radgyr_swarm_regular_sdev = [0.0 for i in range(len(year_swarm_temp))]
            
            for j in range(len(radgyr_swarm_temp)):
                radgyr_swarm_regular[j].append(round(radgyr_swarm_temp[j],2))
#               print 'j, radgyr_swarm_regular[j]: ', j, radgyr_swarm_regular[j]
            
#    print 'radgyr_swarm_regular: ', radgyr_swarm_regular
    
    for i in range(len(radgyr_swarm_regular)):
        radgyr_swarm_regular_mean[i] = np.mean(radgyr_swarm_regular[i])
        radgyr_swarm_regular_sdev[i] = np.std(radgyr_swarm_regular[i])
        
    radgyr_lower_sigma  =   []
    radgyr_higher_sigma =   []
    
    for i in range(len(radgyr_swarm_regular_mean)):
        lower_value = radgyr_swarm_regular_mean[i] - radgyr_swarm_regular_sdev[i]
        radgyr_lower_sigma.append(lower_value)
        higher_value = radgyr_swarm_regular_mean[i] + radgyr_swarm_regular_sdev[i]
        radgyr_higher_sigma.append(higher_value) 

    for i in range(len(radgyr_higher_sigma)):
        if radgyr_lower_sigma[i] < 0.0:
            radgyr_lower_sigma[i] = 0.0
            
        if radgyr_higher_sigma[i] < 0.0:
            radgyr_higher_sigma[i] = 0.0
             
#   -----------------------------------------------------------------
    
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])
    
    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]

#   -----------------------------------------------------------------

    earthquake_year     =   []
    earthquake_mags     =   []
    earthquake_date     =   []
    earthquake_time     =   []
    
    for i in range(len(years)):
        if float(mag[i]) >= 6.0:
            earthquake_year.append(years[i])
            earthquake_mags.append(mag[i])
            earthquake_date.append(date[i])
            earthquake_time.append(time[i])
    mark_size = 2

    comma = ','
    large_earthquakes_output_file = open('Large_Earthquakes.csv', 'w')
    print('Date', comma,'Year',comma, 'Time', comma, 'Magnitude', file=large_earthquakes_output_file)
    for i in range(len(earthquake_year)):
        print(earthquake_date[i],comma,earthquake_year[i],comma, \
                earthquake_time[i], comma, earthquake_mags[i], file=large_earthquakes_output_file)
    large_earthquakes_output_file.close()
                
#   -----------------------------------------------------------------
#   Return plot arrays to non-equidistant data and use points

    year_swarm_neq          =   []
    radgyr_swarm_mean_neq   =   []
    radgyr_swarm_lsig_neq   =   []
    radgyr_swarm_hsig_neq   =   []

    year_swarm_neq.append(year_swarm_regular[0])
    radgyr_swarm_mean_neq.append(radgyr_swarm_regular_mean[0])
    radgyr_swarm_lsig_neq.append(radgyr_lower_sigma[0])
    radgyr_swarm_hsig_neq.append(radgyr_higher_sigma[0])
    
    for i in range(1,len(year_swarm_regular)):
        diff_values = abs(radgyr_swarm_regular_mean[i] - radgyr_swarm_regular_mean[i-1])
        if diff_values >= 0.01:
            year_swarm_neq.append(year_swarm_regular[i])
            radgyr_swarm_mean_neq.append(radgyr_swarm_regular_mean[i])
            radgyr_swarm_lsig_neq.append(radgyr_lower_sigma[i])
            radgyr_swarm_hsig_neq.append(radgyr_higher_sigma[i])
            
#   -----------------------------------------------------------------
#
#   Determine values of radgyr at the time the large M>7 earthquakes occur

    radgyr_M7_EQ    = [0.0 for kk in range(4)]
    year_M7_EQ      = [0.0 for kk in range(4)]
    
    kk = -1

    for ii in range(len(earthquake_year)):
    
        if float(earthquake_mags[ii]) >= 7.0:
        
            kk += 1
        
            for jj in range(len(year_swarm_neq)):
                if year_swarm_neq[jj] < earthquake_year[ii]-0.0192:    # Use date of M>7 earthquake minus 1 week since error in
                                                                       #    year-to-day conversion accumulates over time
                    year_M7_EQ[kk]      = year_swarm_neq[jj]
                    radgyr_M7_EQ[kk]    = radgyr_swarm_mean_neq[jj]
                    
    mean_value_radgyr = np.mean(radgyr_M7_EQ)
    standard_error_of_mean = np.std(radgyr_M7_EQ)
    
    print('')
    print('------------------------------------------------')
    print('')
    print('Year of M>7 EQ: ', earthquake_year)
    print('')
    print('Year of Time Series Value 1 Week Prior to M>7 EQ: ', year_M7_EQ)
    print('')
    print('Radius of Gyration Prior to M>7 EQ: ', radgyr_M7_EQ)
    print('')
    print('Mean + Std Error of the Mean Radius: ', str(round(mean_value_radgyr,4)) + ' Km +/- ' + \
            str(round(standard_error_of_mean,4)) + ' Km')
    print('')
    print('------------------------------------------------')
    print('')

#   -----------------------------------------------------------------
        
#     plt.plot(year_swarm_regular,radgyr_swarm_regular_mean, 'b-', lw=1.0, zorder=3)
#     plt.plot(year_swarm_regular,radgyr_lower_sigma, 'c-', lw=0.5, zorder=3)
#     plt.plot(year_swarm_regular,radgyr_higher_sigma, 'c-', lw=0.5, zorder=3)

#     plt.plot(year_swarm_neq,radgyr_swarm_mean_neq, 'b--', lw=0.6, zorder=2)
#     plt.plot(year_swarm_neq,radgyr_swarm_mean_neq, '.', color='b', ms=3, lw=0.5, zorder=4)
    plt.plot(year_swarm_neq,radgyr_swarm_mean_neq, 'b-', lw=1.0, zorder=5)
    plt.plot(year_swarm_neq,radgyr_swarm_lsig_neq, 'm-', lw=0.4, zorder=4)
    plt.plot(year_swarm_neq,radgyr_swarm_hsig_neq, 'm-', lw=0.4, zorder=4)
    
    plt.gca().invert_yaxis()
            
    plt.grid(True, lw = 0.5, linestyle='dotted', zorder=10, axis = 'y')
    
    xmin,xmax = plt.xlim()
#   xmin,xmax = plt.xlim(1990,xmax)
    ymin, ymax = plt.ylim()
    
#     print 'ymin, ymax: ', ymin, ymax
    
    min_plot_line = [ymin for i in range(len(year_swarm_neq))]
    
    basis_line = [-1.1 for i in range(len(year_swarm_neq))]
    
    zero_line = [0.0 for i in range(len(year_swarm_neq))]
    
    half_line = [-0.5 for i in range(len(year_swarm_neq))]
    
    #   cyan fill
    plt.fill_between(year_swarm_neq,min_plot_line, radgyr_swarm_hsig_neq, color='c', alpha=0.1, zorder=0)
    
    #   white opaque mask at the top
#     plt.fill_between(year_swarm_neq,basis_line, radgyr_swarm_hsig_neq, color='w', alpha=1.0, zorder=2)
    plt.fill_between(year_swarm_neq,radgyr_swarm_lsig_neq, radgyr_swarm_hsig_neq, color='w', alpha=1.0, zorder=2)
    
    #   fill between standard deviation curves
    plt.fill_between(year_swarm_neq,radgyr_swarm_lsig_neq, radgyr_swarm_hsig_neq, color='chartreuse', alpha=0.3, zorder=3)

    
    for i in range(len(earthquake_year)):
        x_eq = [earthquake_year[i], earthquake_year[i]]
        y_eq = [ymin,-0.5]
#                         
        if float(earthquake_mags[i]) >= 6.0 and float(earthquake_mags[i]) < 7.0  and float(earthquake_year[i]) >= cutoff_start_year:
            plt.plot(x_eq, y_eq, linestyle='dotted', color='k', lw=0.7, zorder=1)
            
        if float(earthquake_mags[i]) >= 7.0 and float(earthquake_year[i]) >= cutoff_start_year:
            plt.plot(x_eq, y_eq, color='r', linestyle='--', lw=0.7, zorder=1)

    plt.plot(year_swarm_neq, zero_line, linestyle='dotted', lw=0.5, zorder=10)
            
    xmin,xmax = plt.xlim()
    ymin, ymax = plt.ylim()
#     
    plt.minorticks_on()
    
    SupTitle_text = 'EMA Radius of Gyration for Swarms M $\geq$' + str(completeness_mag) + ' vs. Time Within ' + str(Circle_Radius) + ' Km of ' +  City_Location

    plt.suptitle(SupTitle_text, fontsize=10, y = 0.94)
    
    if ratio_limit > -100.0:
        Title_text = 'N $\geq$ ' + str(burst_min_size) + ' Events;  Depth $\leq$ '\
                +   str(earthquake_depth) + ' Km; N_Steps = ' + str(N_Steps) \
                + '; ENF = ' + str(ratio_limit) + '; CLF = ' + str(sd_factor)
    else:
        Title_text = 'Swarms > ' + str(burst_min_size) + ' Events Within ' + str(Circle_Radius) + ' Km, at Depth < '\
                +   str(earthquake_depth) + ' Km, with N_Steps = ' + str(N_Steps) \
                + ' and sd_factor = ' + str(sd_factor)
            
#     plt.title(Title_text, fontsize=9)
    
    plt.ylabel('Radius of Gyration (Km), Inverted Scale', fontsize = 9)
    plt.xlabel('Time (Year)', fontsize = 9)
    
#     text_x      = (xmax-xmin)*0.10 + xmin
#     text_y      = (ymax-ymin) * 0.93 + ymin
#     text_string = 'Circle Size Proportional to Swarm Radius of Gyration'
#     plt.text(text_x,text_y,text_string, fontsize=8)    
    
    figure_name = './Data/' + Circle_Location + '_RadGyr_Time_EMA_Dev' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

    return 
    
    ######################################################################
    
def map_swarm_centroids(Location, burst_index, mag, lat, lng, date, time, years, cutoff_start_year, cutoff_start_date):
                
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])

#     mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = read_circle_catalog()
    
    #   ---------------------
    
    #   -------------------------------------------------------------
    #
    
    radgyr_list     =   []
    centroid_lat    =   []
    centroid_lng    =   []
    year_swarm      =   []
    
    for i in range(len(burst_index)):
        index_list = burst_index[i]
        radius_gyration, burst_centroid_lat, burst_centroid_lng = \
            BSTUtilities.compute_burst_radius_gyration(index_list, mag, lat, lng, date, time, years)
            
        radgyr_list.append(radius_gyration)
        centroid_lat.append(burst_centroid_lat)
        centroid_lng.append(burst_centroid_lng)
        year_swarm.append(years[index_list[0]])
        
#     print radgyr_list
#     print year_swarm

    #   -------------------------------------------------------------
    #
    
#   Earthquake Dates:  [Landers, Hector Mine, El Major Cucupah, Ridgecrest]
    earthquake_dates    =   []
    

    earthquake_lat  =   []
    earthquake_lng  =   []
    earthquake_mag  =   []
    earthquake_year =   []
    
    for i in range(len(years)):
        if float(mag[i]) >= 6.0:
            earthquake_lat.append(lat[i])
            earthquake_lng.append(lng[i])
            earthquake_mag.append(mag[i])
            earthquake_year.append(years[i])
            
            if float(mag[i]) >= 7.0:
                earthquake_dates.append(years[i])
    
#     color_start =   0.0
#     color_stop  =   1.0
#     
#     cm_subsection = linspace(color_start, color_stop, len(earthquake_dates)) 
#     colors = [ cm.rainbow(x) for x in cm_subsection ]
    
    centroid_lat_list   =   []
    centroid_lng_list   =   []
    centroid_yrs_list   =   []
    
    color_index         =   ['c' for i in range(len(burst_index))]
    
    
    for i in range(len(burst_index)):
    
        lat_list = []
        lng_list = []
        mag_list = []
        
        for j in range(len(burst_index[i])):

            lat_list.append(lat[burst_index[i][j]] )
            lng_list.append(lng[burst_index[i][j]] )
            mag_list.append(mag[burst_index[i][j]] )
        
    #.................................................................
    
        cent_lat = mean_val(lat_list)
        cent_lng = mean_val(lng_list)
        cent_year= years[burst_index[i][0]]
        
        centroid_lat_list.append(cent_lat)
        centroid_lng_list.append(cent_lng)
        centroid_yrs_list.append(cent_year)
        
        
        if centroid_yrs_list[i] < earthquake_dates[0] and centroid_yrs_list[i] >= float(cutoff_start_year):
            color_index[i] = 'c'
            
        elif centroid_yrs_list[i] >= earthquake_dates[0] and centroid_yrs_list[i] < earthquake_dates[1]:
            color_index[i] = 'g'
            
        elif centroid_yrs_list[i] >= earthquake_dates[1] and centroid_yrs_list[i] < earthquake_dates[2]:
            color_index[i] = 'b'

        elif centroid_yrs_list[i] >= earthquake_dates[2] and centroid_yrs_list[i] < earthquake_dates[3]:
            color_index[i] = 'r'

        elif centroid_yrs_list[i] >= earthquake_dates[3]:
            color_index[i] = 'k'
            
#     print centroid_yrs_list
        
#         print 'len(burst_index), i, cent_lat, cent_lng: ', i, cent_lat, cent_lng
        
    #.................................................................
    
    #
    #   Draw map of burst centroids
    
    lat_0_center = Circle_Lat
    lon_0_center = Circle_Lng
    
    delta_lat = change_in_latitude(Circle_Radius)
    delta_lat = abs(delta_lat)
    delta_lng = change_in_longitude(Circle_Lat, Circle_Radius)
    delta_lng = abs(delta_lng)
    
    llcrnrlat  = SWLat_local = Circle_Lat - delta_lat - 0.5
    llcrnrlon  = SWLng_local = Circle_Lng - delta_lng - 0.5
    urcrnrlat  = NELat_local = Circle_Lat + delta_lat + 0.5
    urcrnrlon  = NELng_local = Circle_Lng + delta_lng + 0.5

    m = Basemap(projection='cyl',llcrnrlat=SWLat_local, urcrnrlat=NELat_local,
            llcrnrlon=SWLng_local, urcrnrlon=NELng_local, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,2.0),labels=[0,0,0,1], color='k', textcolor='k', linewidth=0.4, fontsize=8)
    m.drawparallels(np.arange(-90,90,2.0),labels=[1,0,0,0],color='k', textcolor='k', linewidth=0.4, fontsize=8)

    m.fillcontinents(color='tan', lake_color='aqua')
    try:
        m.drawcoastlines(linewidth=0.4, linestyle='solid', color='k')
    except:
        pass
    m.drawcountries()
    m.drawstates()
    m.drawrivers()
    
    mark_size = 1.5
    
    #   ------------------------------------------------------
    
    for i in range(len(color_index)):
    
        if centroid_yrs_list[i] >= cutoff_start_year:
    
            cent_lat = [centroid_lat_list[i]]
            cent_lng = [centroid_lng_list[i]]
        
            m.plot(cent_lng, cent_lat, "o", color=color_index[i], ms=mark_size)
#       m.plot(centroid_lng_list, centroid_lat_list, "ro", ms=mark_size)
#     ax0.plot(centroid_lng_list, centroid_lat_list, "o", ms=mark_size, fillstyle='none', color='k', lw=0.1)

    x_circle_dg, y_circle_dg = BSTUtilities.createCircleAroundWithRadius(Circle_Lat, Circle_Lng, Circle_Radius)
    m.plot(x_circle_dg, y_circle_dg, "b--", lw=0.9)
    
    #   ------------------------------------------------------

    for i in range(len(earthquake_lat)):
        if float(earthquake_mag[i])  >= 6.0 and float(earthquake_mag[i]) < 7.0 and earthquake_year[i] >= float(cutoff_start_year):
#                 m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=4, color = 'k')
                m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=5, mfc='y', mec='k', mew=0.5)
        elif float(earthquake_mag[i]) >= 7.0 and earthquake_year[i] >= float(cutoff_start_year):
#                 m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=8, color = 'k')
                m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=8, mfc='y', mec='k', mew=0.5)
            
    #   ------------------------------------------------------

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    first_index         = burst_index[0][0]
    second_index        = len(burst_index) - 1
    last_event_index    = len(burst_index[second_index]) -1
    second_index    = burst_index[second_index][last_event_index]
    
    start_date = date[first_index]
    end_date   = date[second_index]
    
    number_bursts_since_start = 0
    for i in range(len(burst_index)):
        if centroid_yrs_list[i] >= cutoff_start_year:
            number_bursts_since_start += 1
    
    if Circle_Location != 'None':
        Circle_Location_actual = Circle_Location.split('-')
        SupTitle_text = 'Centroids for ' + str(number_bursts_since_start) + ' Burst Events near ' + City_Location 

    plt.suptitle(SupTitle_text, y=0.97, fontsize=12)
    # 

    Title_text =  'From: ' + cutoff_start_date + '   To: ' + end_date 
    plt.title(Title_text, fontsize=10)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')
    

    figure_name = './Data/' + Circle_Location + '_Swarm_Centroids' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

#    plt.show()

    return None
    
    ######################################################################
    
def map_swarm_centroids_alternate(Location, burst_index, mag, lat, lng, date, time, years, \
        cutoff_start_year, cutoff_start_date, cutoff_radgyr,sd_factor, ratio_limit):
                
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])

#     mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = read_circle_catalog()
    

    #   -------------------------------------------------------------
    #
    
    radgyr_list     =   []
    centroid_lat    =   []
    centroid_lng    =   []
    year_swarm      =   []
    
    for i in range(len(burst_index)):
        index_list = burst_index[i]
        radius_gyration, burst_centroid_lat, burst_centroid_lng = \
            BSTUtilities.compute_burst_radius_gyration(index_list, mag, lat, lng, date, time, years)
            
        radgyr_list.append(radius_gyration)
        centroid_lat.append(burst_centroid_lat)
        centroid_lng.append(burst_centroid_lng)
        year_swarm.append(years[index_list[0]])
        
#     print radgyr_list
#     print year_swarm

    #   -------------------------------------------------------------
    
    
    earthquake_dates    =   []
    

    earthquake_lat  =   []
    earthquake_lng  =   []
    earthquake_mag  =   []
    earthquake_year =   []
    
    for i in range(len(years)):
        if float(mag[i]) >= 6.0:
            earthquake_lat.append(lat[i])
            earthquake_lng.append(lng[i])
            earthquake_mag.append(mag[i])
            earthquake_year.append(years[i])
            
            if float(mag[i]) >= 7.0:
                earthquake_dates.append(years[i])
                
    number_red_circles = 0
    
    color_index = ['None' for i in range(len(year_swarm))]
                
    for i in range(len(year_swarm)):
    
#     color_start =   0.0
#     color_stop  =   1.0
#     
#     cm_subsection = linspace(color_start, color_stop, len(earthquake_dates)) 
#     colors = [ cm.rainbow(x) for x in cm_subsection ]

        if year_swarm[i] < earthquake_dates[0] and year_swarm[i] >= float(cutoff_start_year) \
                and float(radgyr_list[i]) <= cutoff_radgyr:
            color_index[i] = 'c'
            
        elif year_swarm[i] >= earthquake_dates[0] and year_swarm[i] < earthquake_dates[1] \
                and float(radgyr_list[i]) <= cutoff_radgyr:
            color_index[i] = 'g'
            
        elif year_swarm[i] >= earthquake_dates[1] and year_swarm[i] < earthquake_dates[2] \
                and float(radgyr_list[i]) <= cutoff_radgyr:
            color_index[i] = 'b'

        elif year_swarm[i] >= earthquake_dates[2] and year_swarm[i] < earthquake_dates[3] \
                and float(radgyr_list[i]) <= cutoff_radgyr:
            color_index[i] = 'r'
            
            number_red_circles += 1
            
        elif year_swarm[i] >= earthquake_dates[3] and float(radgyr_list[i]) <= cutoff_radgyr:
            color_index[i] = 'k'

    #.................................................................
    
    #
    #   Draw map of burst centroids
    
    lat_0_center = Circle_Lat
    lon_0_center = Circle_Lng
    
    delta_lat = change_in_latitude(Circle_Radius)
    delta_lat = abs(delta_lat)
    delta_lng = change_in_longitude(Circle_Lat, Circle_Radius)
    delta_lng = abs(delta_lng)
    
    llcrnrlat  = SWLat_local = Circle_Lat - delta_lat - 0.5
    llcrnrlon  = SWLng_local = Circle_Lng - delta_lng - 0.5
    urcrnrlat  = NELat_local = Circle_Lat + delta_lat + 0.5
    urcrnrlon  = NELng_local = Circle_Lng + delta_lng + 0.5

    m = Basemap(projection='cyl',llcrnrlat=SWLat_local, urcrnrlat=NELat_local,
            llcrnrlon=SWLng_local, urcrnrlon=NELng_local, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,2.0),labels=[0,0,0,1], color='k', textcolor='k', linewidth=0.4, fontsize=8)
    m.drawparallels(np.arange(-90,90,2.0),labels=[1,0,0,0],color='k', textcolor='k', linewidth=0.4, fontsize=8)

    m.fillcontinents(color='tan', lake_color='aqua')
    try:
        m.drawcoastlines(linewidth=0.4, linestyle='solid', color='k')
    except:
        pass
    m.drawcountries()
    m.drawstates()
    m.drawrivers()
    
    mark_size = 2
    
    #   ------------------------------------------------------
    
    for i in range(len(year_swarm)):
    
        if year_swarm[i] >= cutoff_start_year:
    
            cent_lat = [centroid_lat[i]]
            cent_lng = [centroid_lng[i]]
            
            m.plot(cent_lng, cent_lat, "o", color=color_index[i], ms=mark_size)
#       m.plot(centroid_lng_list, centroid_lat_list, "ro", ms=mark_size)
#     ax0.plot(centroid_lng_list, centroid_lat_list, "o", ms=mark_size, fillstyle='none', color='k', lw=0.1)

    x_circle_dg, y_circle_dg = BSTUtilities.createCircleAroundWithRadius(Circle_Lat, Circle_Lng, Circle_Radius)
    m.plot(x_circle_dg, y_circle_dg, "b--", lw=0.9)
    
    #   ------------------------------------------------------

    for i in range(len(earthquake_lat)):
        if float(earthquake_mag[i])  >= 6.0 and float(earthquake_mag[i]) < 7.0 and earthquake_year[i] >= float(cutoff_start_year):
#                 m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=4, color = 'k')
                m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=5, mfc='y', mec='k', mew=0.5)
        elif float(earthquake_mag[i]) >= 7.0 and earthquake_year[i] >= float(cutoff_start_year):
#                 m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=8, color = 'k')
                m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=8, mfc='y', mec='k', mew=0.5)
            
    #   ------------------------------------------------------

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    first_index         = burst_index[0][0]
    second_index        = len(burst_index) - 1
    last_event_index    = len(burst_index[second_index]) -1
    second_index    = burst_index[second_index][last_event_index]
    
    start_date = date[first_index]
    end_date   = date[second_index]
    
    number_bursts_since_start = 0
    for i in range(len(year_swarm)):
        if year_swarm[i] >= cutoff_start_year and float(radgyr_list[i]) <= cutoff_radgyr:
            number_bursts_since_start += 1
    
    if Circle_Location != 'None':
        Circle_Location_actual = Circle_Location.split('-')
        SupTitle_text = 'Centroids for ' + str(number_bursts_since_start) + ' Burst Events near ' + City_Location 

    plt.suptitle(SupTitle_text, y=0.97, fontsize=12)
    # 

    Title_text =  'From: ' + cutoff_start_date + '   To: ' + end_date + ' for Cutoff Radius: ' + str(cutoff_radgyr) + ' Km; ' +\
            'ENF: ' + str(ratio_limit) + ' Km$^{-1}$; ' + 'CNF: '+ str(sd_factor) + ' Km'
    plt.title(Title_text, fontsize=8)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')
    

    figure_name = './Data/' + Circle_Location + '_Swarm_Centroids_Alt' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

#    plt.show()

    return None
    
    ######################################################################
    
def map_swarm_centroids_dev(Location, burst_index, mag, lat, lng, date, time, years, \
        cutoff_start_year, cutoff_start_date, cutoff_radgyr,sd_factor, ratio_limit):
                
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])

#     mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = read_circle_catalog()
    
    #   -------------------------------------------------------------
    #
    
    radgyr_list     =   []
    centroid_lat    =   []
    centroid_lng    =   []
    year_swarm      =   []
    
    for i in range(len(burst_index)):
        index_list = burst_index[i]
        radius_gyration, burst_centroid_lat, burst_centroid_lng = \
            BSTUtilities.compute_burst_radius_gyration(index_list, mag, lat, lng, date, time, years)
            
        radgyr_list.append(radius_gyration)
        centroid_lat.append(burst_centroid_lat)
        centroid_lng.append(burst_centroid_lng)
        year_swarm.append(years[index_list[0]])
        
#     print radgyr_list
#     print year_swarm

    #   -------------------------------------------------------------
    
    earthquake_dates    =   []
    

    earthquake_lat  =   []
    earthquake_lng  =   []
    earthquake_mag  =   []
    earthquake_year =   []
    color_index     =   []
    
    for i in range(len(years)):
        if float(mag[i]) >= 6.0:
            earthquake_lat.append(lat[i])
            earthquake_lng.append(lng[i])
            earthquake_mag.append(mag[i])
            earthquake_year.append(years[i])
            
            if float(mag[i]) >= 7.0:
                earthquake_dates.append(years[i])
                
    color_index = ['None' for i in range(len(year_swarm))]
    
    number_centroids = 0

    for i in range(len(year_swarm)):
    
#     color_start =   0.0
#     color_stop  =   1.0
#     
#     cm_subsection = linspace(color_start, color_stop, len(earthquake_dates)) 
#     colors = [ cm.rainbow(x) for x in cm_subsection ]
    
#         if year_swarm[i] < earthquake_dates[0] and year_swarm[i] >= float(cutoff_start_year) \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:
#             color_index[i] = 'c'
            
#         if year_swarm[i] >= earthquake_dates[0] and year_swarm[i] < earthquake_dates[1] \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:
#             color_index[i] = 'g'
            
#         if year_swarm[i] >= earthquake_dates[1] and year_swarm[i] < earthquake_dates[2] \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:
#             color_index[i] = 'b'

#         if year_swarm[i] >= earthquake_dates[2] and year_swarm[i] < earthquake_dates[3] \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:

        if year_swarm[i] >= earthquake_dates[3]-5.0 and year_swarm[i] < earthquake_dates[3] \
                and float(radgyr_list[i]) <= cutoff_radgyr:

#             color_index[i] = 'r'
            number_centroids += 1

#         if year_swarm[i] >= earthquake_dates[3] and float(radgyr_list[i]) <= cutoff_radgyr:
#             color_index[i] = 'k'

    color_start =   0.0
    color_stop  =   1.0
    
    cm_subsection = linspace(color_start, color_stop, number_centroids) 
    colors = [ cm.rainbow(x) for x in cm_subsection ]
        
    #.................................................................
    
    #
    #   Draw map of burst centroids
    
    lat_0_center = Circle_Lat
    lon_0_center = Circle_Lng
    
    delta_lat = change_in_latitude(Circle_Radius)
    delta_lat = abs(delta_lat)
    delta_lng = change_in_longitude(Circle_Lat, Circle_Radius)
    delta_lng = abs(delta_lng)
    
    llcrnrlat  = SWLat_local = Circle_Lat - delta_lat - 0.5
    llcrnrlon  = SWLng_local = Circle_Lng - delta_lng - 0.5
    urcrnrlat  = NELat_local = Circle_Lat + delta_lat + 0.5
    urcrnrlon  = NELng_local = Circle_Lng + delta_lng + 0.5

    m = Basemap(projection='cyl',llcrnrlat=SWLat_local, urcrnrlat=NELat_local,
            llcrnrlon=SWLng_local, urcrnrlon=NELng_local, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,2.0),labels=[0,0,0,1], color='k', textcolor='k', linewidth=0.4, fontsize=8)
    m.drawparallels(np.arange(-90,90,2.0),labels=[1,0,0,0],color='k', textcolor='k', linewidth=0.4, fontsize=8)

    m.fillcontinents(color='tan', lake_color='aqua')
    try:
        m.drawcoastlines(linewidth=0.4, linestyle='solid', color='k')
    except:
        pass
    m.drawcountries()
    m.drawstates()
    m.drawrivers()
    
    mark_size = 4
    
    #   ------------------------------------------------------
    
    kk = 0
    
    print('')
    for i in range(len(year_swarm)):
    
#         if year_swarm[i] >= earthquake_dates[2] and year_swarm[i] < earthquake_dates[3] \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:
#     
#         if year_swarm[i] >= cutoff_start_year:

        if year_swarm[i] >= earthquake_dates[3]-5.0 and year_swarm[i] < earthquake_dates[3] \
                and float(radgyr_list[i]) <= cutoff_radgyr:
    
            cent_lat = [centroid_lat[i]]
            cent_lng = [centroid_lng[i]]
            
            centroid_color = colors[kk]
#             centroid_color = 'r'
            
            kk += 1
            
            print('Lat, Lng: ', cent_lat, cent_lng)
        
#             m.plot(cent_lng, cent_lat, "o", color=color_index[i], ms=mark_size)
#             m.plot(cent_lng, cent_lat, "o", color=centroid_color, ms=mark_size)
#             m.plot(cent_lng, cent_lat, "o", mfc='none', mec='r', mew = 0.5, ms=mark_size)
            m.plot(cent_lng, cent_lat, "o", mfc=centroid_color, mec='k', mew = 0.5, ms=mark_size)
#             m.plot(cent_lng, cent_lat, "o", color='r', ms=mark_size)

    print('')
    print('---------------------------------------')
    print('')
    print('Todal Number of Pre-EQ Events: ', kk)
    print('')
    print('---------------------------------------')
    print('')

            
#       m.plot(centroid_lng_list, centroid_lat_list, "ro", ms=mark_size)
#     ax0.plot(centroid_lng_list, centroid_lat_list, "o", ms=mark_size, fillstyle='none', color='k', lw=0.1)

    x_circle_dg, y_circle_dg = BSTUtilities.createCircleAroundWithRadius(Circle_Lat, Circle_Lng, Circle_Radius)
    m.plot(x_circle_dg, y_circle_dg, "b--", lw=0.9)
    
    #   ------------------------------------------------------

    mark_size = 2

    for i in range(len(earthquake_lat)):
        if float(earthquake_mag[i])  >= 6.0 and float(earthquake_mag[i]) < 7.0 and earthquake_year[i] >= float(cutoff_start_year):
#                 m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=4, color = 'k')
                m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=5, mfc='y', mec='k', mew=0.5)
        elif float(earthquake_mag[i]) >= 7.0 and earthquake_year[i] >= float(cutoff_start_year):
#                 m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=8, color = 'k')
                m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=8, mfc='y', mec='k', mew=0.5)
            
    #   ------------------------------------------------------

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    first_index         = burst_index[0][0]
    second_index        = len(burst_index) - 1
    last_event_index    = len(burst_index[second_index]) -1
    second_index    = burst_index[second_index][last_event_index]
    
    start_date = date[first_index]
    end_date   = date[second_index]
    
    number_bursts_since_start = 0
    for i in range(len(year_swarm)):
        if year_swarm[i] >= cutoff_start_year and float(radgyr_list[i]) <= cutoff_radgyr:
            number_bursts_since_start += 1
    
    if Circle_Location != 'None':
        Circle_Location_actual = Circle_Location.split('-')
        SupTitle_text = 'Centroids for ' + str(number_bursts_since_start) + ' Burst Events near ' + City_Location 

    plt.suptitle(SupTitle_text, y=0.97, fontsize=12)
    # 

    Title_text =  'From: ' + cutoff_start_date + '   To: ' + end_date + ' for Cutoff Radius: ' + str(cutoff_radgyr) + ' Km; ' +\
            'ENF: ' + str(ratio_limit) + ' Km$^{-1}$; ' + 'CNF: '+ str(sd_factor) + ' Km'
    plt.title(Title_text, fontsize=8)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')
    

    figure_name = './Data/' + Circle_Location + '_Swarm_Centroids_Dev' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

#    plt.show()

    return None
    
    ######################################################################
    
def map_swarm_centroids_dev1(Location, burst_index, mag, lat, lng, date, time, years, \
        cutoff_start_year, cutoff_start_date, cutoff_radgyr,sd_factor, ratio_limit, offset_year, Large_EQ):
                
    settings_params = get_settings()

    region_type = settings_params[0]
    earthquake_depth    =   float(settings_params[2])
    Circle_Location     =   settings_params[3]
    Circle_Lat          =   float(settings_params[4])
    Circle_Lng          =   float(settings_params[5])
    Circle_Radius        =   float(settings_params[6])

#     mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = read_circle_catalog()
    
    #   -------------------------------------------------------------
    #
    
    radgyr_list     =   []
    centroid_lat    =   []
    centroid_lng    =   []
    year_swarm      =   []
    
    for i in range(len(burst_index)):
        index_list = burst_index[i]
        radius_gyration, burst_centroid_lat, burst_centroid_lng = \
            BSTUtilities.compute_burst_radius_gyration(index_list, mag, lat, lng, date, time, years)
            
        radgyr_list.append(radius_gyration)
        centroid_lat.append(burst_centroid_lat)
        centroid_lng.append(burst_centroid_lng)
        year_swarm.append(years[index_list[0]])
        
#     print radgyr_list
#     print year_swarm

    #   -------------------------------------------------------------
    
    earthquake_dates    =   []
    

    earthquake_lat  =   []
    earthquake_lng  =   []
    earthquake_mag  =   []
    earthquake_year =   []
    color_index     =   []
    LEQ_date        =   []
    LEQ_Burst_Index =   []
    
    for i in range(len(years)):
        if float(mag[i]) >= 6.0:
            earthquake_lat.append(lat[i])
            earthquake_lng.append(lng[i])
            earthquake_mag.append(mag[i])
            earthquake_year.append(years[i])
            
            if float(mag[i]) >= 7.0:
                earthquake_dates.append(years[i])
                LEQ_date.append(date)
                LEQ_Burst_Index.append(i)

    if Large_EQ == 'Landers':
        LEQ_Index = 0
    elif Large_EQ == 'Hector Mine':
        LEQ_Index = 1
    elif Large_EQ == 'El Major Cucupah':
        LEQ_Index = 2
    elif Large_EQ == 'Ridgecrest':
        LEQ_Index = 3
                
    color_index = ['None' for i in range(len(year_swarm))]
    
    number_centroids = 0

    for i in range(len(year_swarm)):
    
#     color_start =   0.0
#     color_stop  =   1.0
#     
#     cm_subsection = linspace(color_start, color_stop, len(earthquake_dates)) 
#     colors = [ cm.rainbow(x) for x in cm_subsection ]
    
#         if year_swarm[i] < earthquake_dates[0] and year_swarm[i] >= float(cutoff_start_year) \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:
#             color_index[i] = 'c'
            
#         if year_swarm[i] >= earthquake_dates[0] and year_swarm[i] < earthquake_dates[1] \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:
#             color_index[i] = 'g'
            
#         if year_swarm[i] >= earthquake_dates[1] and year_swarm[i] < earthquake_dates[2] \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:
#             color_index[i] = 'b'

#         if year_swarm[i] >= earthquake_dates[2] and year_swarm[i] < earthquake_dates[3] \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:

        if year_swarm[i] >= earthquake_dates[LEQ_Index]-offset_year and year_swarm[i] < earthquake_dates[LEQ_Index] \
                and float(radgyr_list[i]) <= cutoff_radgyr:

#             color_index[i] = 'r'
            number_centroids += 1

#         if year_swarm[i] >= earthquake_dates[3] and float(radgyr_list[i]) <= cutoff_radgyr:
#             color_index[i] = 'k'

    color_start =   0.0
    color_stop  =   1.0
    
    cm_subsection = linspace(color_start, color_stop, number_centroids) 
    colors = [ cm.rainbow(x) for x in cm_subsection ]
        
    #.................................................................
    
    #
    #   Draw map of burst centroids
    
    lat_0_center = Circle_Lat
    lon_0_center = Circle_Lng
    
    delta_lat = change_in_latitude(Circle_Radius)
    delta_lat = abs(delta_lat)
    delta_lng = change_in_longitude(Circle_Lat, Circle_Radius)
    delta_lng = abs(delta_lng)
    
    llcrnrlat  = SWLat_local = Circle_Lat - delta_lat - 0.5
    llcrnrlon  = SWLng_local = Circle_Lng - delta_lng - 0.5
    urcrnrlat  = NELat_local = Circle_Lat + delta_lat + 0.5
    urcrnrlon  = NELng_local = Circle_Lng + delta_lng + 0.5

    m = Basemap(projection='cyl',llcrnrlat=SWLat_local, urcrnrlat=NELat_local,
            llcrnrlon=SWLng_local, urcrnrlon=NELng_local, lat_0=lat_0_center, lon_0=lat_0_center, lat_ts=20, resolution='h')

    m.drawmeridians(np.arange(0,360,2.0),labels=[0,0,0,1], color='k', textcolor='k', linewidth=0.4, fontsize=8)
    m.drawparallels(np.arange(-90,90,2.0),labels=[1,0,0,0],color='k', textcolor='k', linewidth=0.4, fontsize=8)

    m.fillcontinents(color='tan', lake_color='aqua')
    try:
        m.drawcoastlines(linewidth=0.4, linestyle='solid', color='k')
    except:
        pass
    m.drawcountries()
    m.drawstates()
    m.drawrivers()
    
    mark_size = 4
    
    #   ------------------------------------------------------
    
    kk = 0
    
    cent_lat_reduced    =   []
    cent_lng_reduced    =   []
    radgyr_reduced      =   []
    year_reduced        =   []
    
    print('')
    for i in range(len(year_swarm)):
    
#         if year_swarm[i] >= earthquake_dates[2] and year_swarm[i] < earthquake_dates[3] \
#                 and float(radgyr_list[i]) <= cutoff_radgyr:
#     
#         if year_swarm[i] >= cutoff_start_year:

        if year_swarm[i] >= earthquake_dates[LEQ_Index]-offset_year and year_swarm[i] < earthquake_dates[LEQ_Index] \
                and float(radgyr_list[i]) <= cutoff_radgyr:
    
            cent_lat = [centroid_lat[i]]
            cent_lng = [centroid_lng[i]]
            
            cent_lat_reduced.append(float(centroid_lat[i]))
            cent_lng_reduced.append(float(centroid_lng[i]))
            radgyr_reduced.append(float(radgyr_list[i]))
            year_reduced.append(float(year_swarm[i]))
            
            centroid_color = colors[kk]
#             centroid_color = 'r'
            
            kk += 1
            
#             m.plot(cent_lng, cent_lat, "o", color=color_index[i], ms=mark_size)
#             m.plot(cent_lng, cent_lat, "o", color=centroid_color, ms=mark_size)
#             m.plot(cent_lng, cent_lat, "o", mfc='none', mec='r', mew = 0.5, ms=mark_size)
            m.plot(cent_lng, cent_lat, "o", mfc=centroid_color, mec='k', mew = 0.5, ms=mark_size)
#             m.plot(cent_lng, cent_lat, "o", color='r', ms=mark_size)

    year_lat = list(zip(year_reduced,cent_lat_reduced))
    year_lng = list(zip(year_reduced,cent_lng_reduced))
    year_rgyr= list(zip(year_reduced,radgyr_reduced))
    
    ylat = sorted(year_lat)
    ylng = sorted(year_lng)
    ygyr = sorted(year_rgyr)
    year_sorted = sorted(year_reduced)
    
    RG_Table = []
    
    print('')
    for i in range(len(radgyr_reduced)):
        RG_Table.append([i+1, round(year_sorted[i],4), round(ygyr[i][1],4), round(ylat[i][1],4), round(ylng[i][1],4)])
        
    table_data = tabulate(RG_Table,headers=['Burst', 'Year', 'Radius of Gyration (Km)',\
            'Centroid Latitude (Deg)', 'Centroid Longitude (Deg)'], numalign="center")
    
    print(table_data)
        
    print('')

    x_circle_dg, y_circle_dg = BSTUtilities.createCircleAroundWithRadius(Circle_Lat, Circle_Lng, Circle_Radius)
    m.plot(x_circle_dg, y_circle_dg, "b--", lw=0.9)
    
    #   ------------------------------------------------------

    mark_size = 2

    for i in range(len(earthquake_lat)):
        if float(earthquake_mag[i])  >= 6.0 and float(earthquake_mag[i]) < 7.0 and earthquake_year[i] >= float(cutoff_start_year):
#                 m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=4, color = 'k')
                m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=5, mfc='y', mec='k', mew=0.5)
        elif float(earthquake_mag[i]) >= 7.0 and earthquake_year[i] >= float(cutoff_start_year):
#                 m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=8, color = 'k')
                m.plot(earthquake_lng[i], earthquake_lat[i], "*", ms=8, mfc='y', mec='k', mew=0.5)
            
    #   ------------------------------------------------------

    City = Location.split('_')
    Number_Spaces = len(City)

    Location_Actual = City[0] + ' '
    for i in range(1,Number_Spaces):
        Location_Actual += City[i] + ' '

    City_Location = Location_Actual[:-1]
    
    first_index         = burst_index[0][0]
    second_index        = len(burst_index) - 1
    last_event_index    = len(burst_index[second_index]) -1
    second_index    = burst_index[second_index][last_event_index]
    
    start_date = date[first_index]
    end_date   = date[second_index]
    
    number_bursts_since_start = 0
    for i in range(len(year_swarm)):
        if year_swarm[i] >= cutoff_start_year and float(radgyr_list[i]) <= cutoff_radgyr:
            number_bursts_since_start += 1
    
    if Circle_Location != 'None':
        Circle_Location_actual = Circle_Location.split('-')
        SupTitle_text = 'Centroids for ' + str(kk) + ' Burst Events near ' + City_Location 

    plt.suptitle(SupTitle_text, y=0.97, fontsize=12)
    # 
    
    Title_text =  'For the ' + str(offset_year) + ' Years Prior to ' + \
            date[LEQ_Burst_Index[LEQ_Index]] + ' for Cutoff Radius: ' + str(cutoff_radgyr) + ' Km; ' +\
            'ENF: ' + str(ratio_limit) + ' Km$^{-1}$; ' + 'CNF: '+ str(sd_factor) + ' Km'
    plt.title(Title_text, fontsize=8)

    #m.fillcontinents(color='#ffe8a0',lake_color='#e6e6ff')
    m.drawmapboundary(fill_color='#e6e6ff')
    

    figure_name = './Data/' + Circle_Location + '_Swarm_Centroids_Dev' + '.png'
    matplotlib.pyplot.savefig(figure_name,dpi=600)
    matplotlib.pyplot.close('all')

#    print ' '
#    print '     Close plot window to continue...'

#    plt.show()

    return None
    
