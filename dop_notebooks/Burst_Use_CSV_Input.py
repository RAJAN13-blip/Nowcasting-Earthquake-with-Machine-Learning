#!/opt/local/bin python

    #   Code Burst.py
    #
    #   Python code to compute a nowcast from bursts, swarms and aftershocks
    #
    #   This code downloads data from the USGS web site 
    #
    #   This code was written on a Mac using Macports python.  A list of the ports needed to run the code are available at:
    #       https://www.dropbox.com/s/8wr5su8d7l7a30z/myports-wailea.txt?dl=0
    
    #   ---------------------------------------------------------------------------------------
    
    # Copyright 2018 by John B Rundle, University of California, Davis, CA USA
    # 
    # Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
    # documentation files (the     "Software"), to deal in the Software without restriction, including without 
    # limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
    # and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    # 
    # The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    # 
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
    # WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
    # COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
    # ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    #   ---------------------------------------------------------------------------------------

import sys
import os
import numpy as np
from array import array

import BSTCalcMethods
import BSTFileMethods
import BSTPlotMethods
import BSTUtilities

import datetime
import dateutil.parser

import time
from time import sleep  #   Added a pause of 30 seconds between downloads

import math


    ######################################
    #
    #   Note that the input circle centers around cities are given in the file "circle_locations.txt"

    ######################################

#datetime.date.today().strftime("%F")
date_info = str(datetime.date.today())
date_data = date_info.split('-')
current_date = str(date_data[0]) + '/' + str(date_data[1]) + '/' + str(date_data[2])
current_time_info = str(datetime.datetime.now())
current_time = current_time_info.split(' ')[1]
current_time = str(current_time.split('.')[0])
#
#   ------------------------------------------------------------

#   SoCal
maxlng = -110.0
minlng = -126.
maxlat =  40.0
minlat = 28.0

depth = 100.0                       #   Defines maximum earthquake depth (variable)
km = 2000.0                         #   Radius Defines regional length scale for statistics in km (default=2000)

start_date = "1975/01/01"       #   Events downloaded occurred after this date
# start_date = "1990/01/01"       #   Events downloaded occurred after this date
cutoff_date = 2100.0

circle_catalog_date_start = 1975.0
circle_catalog_date_end   = 2020.999

city_circle = 400.0     #   <-- Was 400.0
completeness_mag_uniform = 2.99     #   This is for the Worldwide catalog
completeness_mag = 3.29     

# cluster_percent = 0.985  #  ~ 3*sigma -> 0.985 for SoCal with 400 km circle M>2.99 centered on SoCal/USA,33.5,-118.5
                                #   0.992 for Kyoto, 350 km centered on Kyoto, M>3.99
sd_factor  = 25        #   This is CLF,it filters cluster outliersResults seem to be insensitive to this value in the range 50-150
                        #       at a value of median ratio_limit for 600 km radius
                    
cutoff_start_year = 1984.0
cutoff_start_date = "1984/1/1"

cutoff_radgyr = 1.0         #   Cutoff radius for plotting burst centroids on the map (Km)
offset_year = 3.0           #   Value of time in years prior to M7 earthquakes we are using for centroid plot

Large_EQ = 'Landers'
Large_EQ = 'Hector Mine'
Large_EQ = 'El Major Cucupah'
Large_EQ = 'Ridgecrest'

#   .......................................................

N_Steps = 23

ratio_limit = 0.2   #   0.2 seems to give interesting results for the centroid_dev method

                    #   (We also find that this value produces the same radius of gyration for all M>7 earthquakes.
                    #       Physical reason this is good is that if this failure radius is to have some physical meaning,
                    #       it should be the same for all large earthquakes This is our fundamental hypothesis, that
                    #       the minimum radius of gyration has physical meaning.

                    
plot_burst_size = 100000000000  #   large number means no burst plots
#plot_burst_size = 10

burst_print_flag = 'OFF'

#   ------------------------------------------------------------

WW_Download = 'ON'  #   Turn this off if you want to speed up code by not using the data download

if WW_Download == 'ON':
    print('')
    print('Downloading the World Wide catalog...')
    print('')
    BSTFileMethods.get_worldwide_catalog(maxlat, minlat, maxlng, minlng, completeness_mag_uniform, start_date)
#
#   -----------------------------------
#
   
    #.....................................
    #
    # Read pre-defined locations file
    #
    #.....................................

input_file = open("city_locations.txt", "r")
i=0
for line in input_file:
    i +=  1
input_file.close()  # Put the file pointer back at top

number_cities = i

print('')
print('Number of Cities: ', number_cities)

    # Create arrays of length i filled with zeros

Location_file   = ["" for x in range(i)]
Country_file    = ["" for x in range(i)]
Center_Lat_file = np.zeros(i)
Center_Lng_file = np.zeros(i)

MinLat_file     = np.zeros(i)
MaxLat_file     = np.zeros(i)
MinLng_file     = np.zeros(i)
MaxLng_file     = np.zeros(i)

MinLat_local_file     = np.zeros(i)
MaxLat_local_file     = np.zeros(i)
MinLng_local_file     = np.zeros(i)
MaxLng_local_file     = np.zeros(i)

    #   Compute the size of the region and store for later use
    #   Start by assuming a circle of radius "km" around
    #       city (lat,long)


input_file = open("city_locations.txt", "r")

i=-1
for line in input_file:
    i+=1
    line    = line.strip()
    items   = line.split(',')

    #   Read in the city, its lat and long

    items_array = np.asarray(items)

    city_country            = items_array[0]

    city                    = city_country.split('/')[0]
    country                 = city_country.split('/')[1]

    Location_file[i]   = city
    Center_Lat_file[i] = items_array[1]
    Center_Lng_file[i] = items_array[2]
    Country_file[i]         = city_country.split('/')[1]
 
    delta_Lng = BSTUtilities.change_in_longitude(Center_Lat_file[i], km)
    delta_Lat = BSTUtilities.change_in_latitude(km)


    delta_local_Lng = BSTUtilities.change_in_longitude(Center_Lat_file[i], city_circle)
    delta_local_Lat = BSTUtilities.change_in_latitude(city_circle)

    MinLat_file[i]     = Center_Lat_file[i] - delta_Lat
    MaxLat_file[i]     = Center_Lat_file[i] + delta_Lat
    MinLng_file[i]     = Center_Lng_file[i] - delta_Lng
    MaxLng_file[i]     = Center_Lng_file[i] + delta_Lng

    MinLat_local_file[i]     = Center_Lat_file[i] - delta_local_Lat
    MaxLat_local_file[i]     = Center_Lat_file[i] + delta_local_Lat
    MinLng_local_file[i]     = Center_Lng_file[i] - delta_local_Lng
    MaxLng_local_file[i]     = Center_Lng_file[i] + delta_local_Lng

input_file.close()  # Put the file pointer back at top

    #.....................................
    
for i in range(0,number_cities):

    #   For US west coast cities, completeness mag = 2.99
    #       for other global cities, completeness mag = 3.99
    #
    #   List US cities first, and execute following if statement

#     completeness_mag    =   completeness_mag_uniform

#    if (i<number_cities and Country_file[i] == 'USA'):
#        completeness_mag = 2.99

    Circle_Location     =   Location_file[i]
    Circle_Lat          =   Center_Lat_file[i]
    Circle_Lng          =   Center_Lng_file[i]
    Radius_float        =   city_circle
    earthquake_depth    =   depth
    region_type         =   'Circle'

    small_radius = str(Radius_float)
    length_scale   = str(2.*km)
    depth        = str(earthquake_depth)

    settings_params     =   []
    settings_params.append(region_type)
    settings_params.append(completeness_mag)
    settings_params.append(earthquake_depth)
    settings_params.append(Circle_Location)
    settings_params.append(Circle_Lat)
    settings_params.append(Circle_Lng)
    settings_params.append(Radius_float)

    BSTUtilities.save_settings(settings_params)

    settings_params = BSTUtilities.get_settings()

    Location    = Location_file[i]
    SWLat       = MinLat_file[i]
    NELat       = MaxLat_file[i]
    SWLng       = MinLng_file[i]
    NELng       = MaxLng_file[i]

    SWLat_local       = MinLat_local_file[i] - 0.25
    NELat_local       = MaxLat_local_file[i] + 0.25
    SWLng_local       = MinLng_local_file[i] - 0.25
    NELng_local       = MaxLng_local_file[i] + 0.25
    
    # Build the regional catalog from the World Wide catalog
#     BSTFileMethods.get_catalog(NELat, NELng, SWLat, SWLng)
    
    # Build the circle catalog from the World Wide catalog
    Rebuild_flag = 'OFF'
    
    BSTFileMethods.get_circle_catalog(Circle_Lat,Circle_Lng,Radius_float,Rebuild_flag, \
        circle_catalog_date_start, circle_catalog_date_end)
        
#     Plot EQ Epicenters on a Map
#     BSTMethods.map_epicenters(NELat, NELng, SWLat, SWLng, Location)
#      
    
    # Plot Local EQ Epicenters on a Map in US
    if (i<number_cities and Country_file[i] == 'USA' or Country_file[i] == 'NL' or Country_file[i] == 'Japan'):
#         completeness_mag = 0.99  #   Plot events down to M>3 for USA
        settings_params = BSTUtilities.get_settings()
        settings_params[1] = completeness_mag
        BSTUtilities.save_settings(settings_params)
        
        # Build the circle catalog again from the World Wide catalog but with smaller completeness mag
        Rebuild_flag = 'ON'
        
        regional_rate = BSTCalcMethods.get_regional_rate(completeness_mag, depth, cutoff_start_year)
        
        BSTFileMethods.get_circle_catalog(Circle_Lat,Circle_Lng,Radius_float,Rebuild_flag, \
                circle_catalog_date_start, circle_catalog_date_end) 

    # Plot Local EQ Epicenters on a Map in US
#         BSTMethods.map_epicenters_local(NELat_local, NELng_local, SWLat_local, SWLng_local, completeness_mag, Location)

        burst_index, number_events, min_daily_events, burst_min_size, mag, lat, lng, date, time, years = \
                    BSTPlotMethods.plot_swarm_event_counts\
                            (completeness_mag, Location, sd_factor, regional_rate, burst_print_flag)
#                     
#         BSTMethods.map_swarm_centroids(Location, burst_index, mag, lat, lng, date, time, years, cutoff_start_year, cutoff_start_date)
        
#       Earthquake Dates:  [Landers, Hector Mine, El Major Cucupah, Ridgecrest]
        
#       Remember to add Landers if not there: 1992/06/28	11:57:35.000	1992.486339	-116.433	34.217	7.3	1.09
        
#         BSTMethods.log_freq_bursts(burst_index, burst_min_size, completeness_mag, Location)
        
#         BSTMethods.map_swarm_epicenters_all(Origin_Location, burst_min_size, completeness_mag, Location, N_Steps, \
#                 burst_index, mag, lat, lng, date, time, years)
                
 
        BSTPlotMethods.plot_radgyr_time(burst_min_size, completeness_mag, Location, \
                circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years)
                
        N_Steps, radgyr_list, year_swarm, index_swarm, centroid_lat, centroid_lng  = \
                BSTPlotMethods.plot_radgyr_EMA_time(burst_min_size, completeness_mag, Location,  N_Steps,\
                circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
                ratio_limit, sd_factor, cutoff_start_year, regional_rate)
                
        BSTPlotMethods.plot_cluster_density_rgyr(burst_min_size, completeness_mag, Location, N_Steps,\
                circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
                ratio_limit, cutoff_start_year)
                
        BSTPlotMethods.map_swarm_centroids(Location, burst_index, mag, lat, lng, date, time, years, cutoff_start_year, cutoff_start_date)
        
        BSTPlotMethods.map_swarm_centroids_alternate(Location, burst_index, mag, lat, lng, date, time, years, \
                cutoff_start_year, cutoff_start_date, cutoff_radgyr, sd_factor, ratio_limit)
                
        BSTPlotMethods.map_swarm_centroids_dev1(Location, burst_index, mag, lat, lng, date, time, years, \
                cutoff_start_year, cutoff_start_date, cutoff_radgyr, sd_factor, ratio_limit, offset_year, Large_EQ)
                
                
        BSTPlotMethods.plot_radgyr_EMA_time_dev(burst_min_size, completeness_mag, Location, N_Steps,\
                circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
                ratio_limit, sd_factor, cutoff_start_year, regional_rate)
                
#         plot_burst_size = burst_min_size * burst_min_size_plot_factor

#   -----------------------------------------------------------------

        if i == 0:
            for j in range(len(burst_index)):
                burst_number = j
                if len(burst_index[j]) >=  plot_burst_size and len(burst_index[j]) <= 15.*plot_burst_size:
                    BSTPlotMethods.map_swarm_epicenters_one_burst(completeness_mag, Location, N_Steps, \
                            burst_index, burst_number, mag, lat, lng, date, time, years)

#   -----------------------------------------------------------------

#   Years of M7 earthquakes = [1992.491803, 1999.789617, 2010.256831, 2019.510929]
#   Plot the bursts that occur 1 year prior to the M>7 earthquakes

#         if i == 0:
#             for j in range(len(burst_index)):
#                 burst_number = j
#                 burst_year = years[burst_index[j][0]]
#                 swarm_flag = 0
#                 
#                 if burst_year >= 1991.5 and burst_year < 1992.4:
# #                     print('Burst Year: ', burst_year)
#                     swarm_flag = 1
#                     
#                 if burst_year >= 1998.78 and burst_year < 1999.7:
# #                     print('Burst Year: ', burst_year)
#                     swarm_flag = 1
#                     
#                 if burst_year >= 2009.25 and burst_year < 2010.2:
# #                     print('Burst Year: ', burst_year)
#                     swarm_flag = 1
#                     
#                 if burst_year >= 2018.5 and burst_year < 2019.4:
# #                     print('Burst Year: ', burst_year)
#                     swarm_flag = 1
#                 
#                 if swarm_flag == 1:
#                     BSTPlotMethods.map_swarm_epicenters_one_burst(completeness_mag, Location, N_Steps, \
#                             burst_index, burst_number, mag, lat, lng, date, time, years)
# 
# #   -----------------------------------------------------------------