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
import numpy as np
from numpy import *
from array import array

import os

import math
from scipy import special
#import scipy.special

import random

import datetime
from osgeo import ogr, osr

import BSTUtilities

    ######################################################################

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
            
        if float(mag) >= float(completeness_mag) and float(dep) <= float(earthquake_depth) and float(year) >= float(cutoff_start_year):
            number_earthquakes += 1
                
    input_file.close()

    time_interval = (float(year) - cutoff_start_year)*365.0
    
    regional_rate = float(number_earthquakes)/time_interval
    
    return regional_rate
    
def EMA_weights(N_events, N_Steps):

    #   This method computes the weights for the Exponential Weighted Average (EMA)

    alpha = 2./float((N_Steps+1))

    #   time_series_list is the time series of floating point values
    #       arranged in order of first element in list being earliest

    assert 0 < alpha <= 1
    
    weights = []
    
    #   Define the weights
    
    for i in range(0,N_events):
        weight_i = (1.0-alpha)**i
        weights.append(weight_i)
        
    sum_weights = sum(weights)
    weights =  [i/sum_weights for i in weights]
     
    return weights
    
    #.................................................................
    
def EMA_weighted_time_series(time_series, N_Steps):

    #   This method computes the Exponential Weighted Average of a list.  Last
    #       in the list elements are exponentially weighted the most

    N_events = len(time_series)
    
    weights = EMA_weights(N_events, N_Steps)
    
    weights_reversed = list(reversed(weights))

    
#     print 'N_events: ', N_events
#     print ''
#     print 'weights_reversed: ', weights_reversed
#     print ''
#     print 'time_series: ', time_series
#     print ''
    
#     print ''
#     print time_series
#     print ''
#     print weights_reversed
#     print ''

    EMA_weighted_ts = []
    partial_weight_sum = 0.
    
    for i in range(N_events):
        partial_weight_sum += weights[i]
        weighted_ts = round(float(time_series[i])*weights_reversed[i],4)
        
        EMA_weighted_ts.append(weighted_ts)
        
    partial_weight_sum = round(partial_weight_sum,4)
    sum_value = sum(EMA_weighted_ts)
    
    if (float(partial_weight_sum)) <= 0.0:
        sum_value = 0.0001
        partial_weight_sum = 1.
    
    weighted_sum = float(sum_value)/float(partial_weight_sum)
    
    return weighted_sum
    
    #.................................................................
    
def build_daily_counts(burst_min_size, min_daily_events, mag, lat, lng, date, time, yrs, sd_factor, burst_print_flag):

    #   First find lat and lng of burst centroid.  Then find standard deviation.
    #
    
    max_duration = 365.0*float((max(yrs) - min(yrs) ))    # In days
    
    number_events = [0. for i in range(int(max_duration))]
    index_events  = [[] for i in range(int(max_duration))]
    
    for i in range(len(yrs)):
        
        day_of_event = 365.0*float(yrs[i] - min(yrs))

        k = int(day_of_event)-1     #   Index of the day
        
        number_events[k] += 1      #   Tabulates number of events on each day
        
        index_events[k].append(i)
        
    burst_index = combine_daily_count_bursts(min_daily_events, number_events, index_events)
    
    #   Clean the data by removing outliers
    
    burst_index =   clean_data_outliers(burst_index, lat, lng, burst_min_size, sd_factor)
    
    #   ---------------------------------------------------
    
    basic_event_rate = float(len(yrs))/(float(max(yrs)-min(yrs))*365.0)
    
    burst_number = 0
    
    total_number_events_in_bursts = 0
    number_events_in_large_bursts = 0
    
#     for i in range(len(burst_index)):
#         print 'i, burst_index: ', i, burst_index[i]

    save_bursts = 'ON'

    if save_bursts == 'ON':
        output_file_name = 'burst_file.txt'
        burst_file = open('burst_file.txt', 'w')
    
    for i in range(len(burst_index)):
    
        total_number_events_in_bursts += int(len(burst_index[i]))
    
        lat_list = []
        lng_list = []
            
        for j in range(len(burst_index[i])):
        
            k = burst_index[i][j]
            lat_list.append(lat[k])
            lng_list.append(lng[k])
        
        burst_centroid_lat = BSTUtilities.mean_val(lat_list)
        burst_centroid_lng = BSTUtilities.mean_val(lng_list)
        
        std_dev_lat, var_lat = BSTUtilities.std_var_val(lat_list)
        std_dev_lng, var_lng = BSTUtilities.std_var_val(lng_list)
        
        first_index = burst_index[i][0]
        k = len(burst_index[i])-1
        second_index = burst_index[i][k]
        
        start_date = date[first_index]
        end_date   = date[second_index]
        
        start_time = time[first_index]
        end_time   = time[second_index]
            
        number_events_burst = int(len(burst_index[i]))
        
        burst_duration =  int( 365.0*float(yrs[second_index] - yrs[first_index]) ) + 1
        
        
        if start_date != end_date:
            burst_duration += 1
   
        duration_rate = float(number_events_burst)/float(burst_duration)
        
        if len(burst_index[i]) >= burst_min_size and burst_print_flag == 'ON':
        
            number_events_in_large_bursts += int(len(burst_index[i]))
        
            print('')
            print('Burst Number:', i, 'Number of Events in Swarm:', len(burst_index[i]))
            print('')
            print('Start Date & Time: ', start_date,' @ ',start_time)
            print('')
            print('  End Date & Time: ', end_date, ' @ ', end_time)
            print('')
            print('Swarm centered at (degrees in lat,lng): ', round(burst_centroid_lat,4), round(burst_centroid_lng,4))
            print('')
            print('Swarm size (standard deviation in degrees lat,lng): ', round(std_dev_lat,4), round(std_dev_lng,4))
            print('')
            print('Swarm duration (days): ', int(burst_duration))
            print('')
            print('Basic rate, Duration rate (per day): ', round(basic_event_rate,4), round(duration_rate,4))
            print('')
            print('Minimum Burst Size in Plots: ', burst_min_size)
            print('')
            print('----------------------------------------')
            
            if save_bursts == 'ON':
                print('', file=burst_file)
                print('Burst Number:', i, 'Number of Events in Swarm:', len(burst_index[i]), file=burst_file)
                print('', file=burst_file)
                print('Start Date & Time: ', start_date,' @ ',start_time, file=burst_file)
                print('', file=burst_file)
                print('  End Date & Time: ', end_date, ' @ ', end_time, file=burst_file)
                print('', file=burst_file)
                print('Swarm centered at (degrees in lat,lng): ', round(burst_centroid_lat,4), round(burst_centroid_lng,4), file=burst_file)
                print('', file=burst_file)
                print('Swarm size (standard deviation in degrees lat,lng): ', round(std_dev_lat,4), round(std_dev_lng,4), file=burst_file)
                print('', file=burst_file)
                print('Swarm duration (days): ', int(burst_duration), file=burst_file)
                print('', file=burst_file)
                print('Basic rate, Duration rate (per day): ', round(basic_event_rate,4), round(duration_rate,4), file=burst_file)
                print('', file=burst_file)
                print('Minimum Burst Size in Plots: ', burst_min_size, file=burst_file)
                print('', file=burst_file)
                print('----------------------------------------', file=burst_file)
        
        burst_number += 1
        
    print('----------------------------------------')
    print('')
    print('Number of Bursts for the Entire USGS Circle Catalog: ', len(burst_index))
    print('')
    print('Percent of Events Contained in Bursts of Any Size: ', str(round(100.0 * \
            float(total_number_events_in_bursts)/float(len(yrs)),2))+'%')
#     print ''
#     print 'Percent of Events Contained in Large Bursts: ', str(round(100.0 * float(number_events_in_large_bursts)/float(len(yrs)),2))+'%'
    print('')
    print('----------------------------------------')
    print('----------------------------------------')

    print('')
    
    if save_bursts == 'OFF':
        burst_file.close()
    
    #   ---------------------------------------------------
    
    return number_events, index_events, burst_index
        
    #   ---------------------------------------------------
    
def combine_daily_count_bursts(min_daily_events, number_events, index_events):

    #   This method combines the daily counts into coherent burst event swarms
    
    temp_burst_index    =   []
    burst_index_prelim  =   []
    burst_index         =   []
    
    for i in range(len(number_events)):
         #   This is a burst
        if i == 0 and number_events[0] >= min_daily_events:    
            for j in range(len(index_events[0])):
                temp_burst_index.append(index_events[0][j])
            
        # Allows possibility of including foreshocks
        elif i > 0 and i< len(number_events)-1 and number_events[i] < min_daily_events and number_events[i+1] >= min_daily_events:
            for j in range(len(index_events[i])):
                temp_burst_index.append(index_events[i][j])
            
        #   Start a new burst
        elif i > 0 and number_events[i] >= min_daily_events and number_events[i-1] < min_daily_events: #   Add this to burst
            for j in range(len(index_events[i])):
                temp_burst_index.append(index_events[i][j])
        
        #   Add this to the current burst
        elif i > 0 and number_events[i] >= min_daily_events and number_events[i-1] >= min_daily_events: #   Add this to burst
            for j in range(len(index_events[i])):
                temp_burst_index.append(index_events[i][j])

        #   Terminate the current burst and reset the temp_burst_index
        elif i > 0 and number_events[i] < min_daily_events and number_events[i-1] >= min_daily_events and len(temp_burst_index) > 0:
            burst_index_prelim.append(temp_burst_index)
            temp_burst_index    =   []
            
        else:
            pass
    
    # Ensure that the number of events in the burst is at least the required minimum
    
    for i in range(len(burst_index_prelim)):
        if len(burst_index_prelim[i]) >= min_daily_events:
            burst_index.append(burst_index_prelim[i])
            
    return  burst_index
    
    #.................................................................
    
def reject_outliers_earthquakes(index_list, data_lat, data_lng, sd_factor):

    #   This method removes lat-lng points that are far from the median of the burst
    
    median_lat = np.median(data_lat)
    median_lng = np.median(data_lng)
    
    great_circle_distance   =   []
    data_lat_removed        =   []
    data_lng_removed        =   []
    index_list_removed      =   []
    
    for i in range(len(data_lat)):
        lat1 = median_lat
        lng1 = median_lng
        lat2 = data_lat[i]
        lng2 = data_lng[i]
        dist = BSTUtilities.compute_great_circle_distance(lat1, lng1, lat2, lng2)
        great_circle_distance.append(dist)
        
    median_gcdist = np.median(great_circle_distance)
    
    for i in range(len(data_lat)):
        if great_circle_distance[i] < sd_factor * median_gcdist:
            data_lat_removed.append(data_lat[i])
            data_lng_removed.append(data_lng[i])
            index_list_removed.append(index_list[i])
            
#     index_list_removed = index_list
            
    return index_list_removed
    
    #.................................................................
    
def clean_data_outliers(burst_index_test, lat, lng, burst_min_size, sd_factor):

    burst_index    =   []
        
    for i in range(len(burst_index_test)):
        index_list = burst_index_test[i]
        lat_list = []
        lng_list = []
            
        for j in range(len(index_list)):
            lat_list.append(lat[index_list[j]])
            lng_list.append(lng[index_list[j]])
    
        index_list_removed = reject_outliers_earthquakes(index_list, lat_list, lng_list, sd_factor)
        
        if len(index_list_removed) >= burst_min_size:
            burst_index.append(index_list_removed)
        
    #   Make sure all bursts are non-empty
            
    return burst_index
    
    #.................................................................
    
def compute_burst_radius_gyration(index_list, mag, lat, lng, date, time, years):
    
    lat_list = []
    lng_list = []
    mag_list = []
    
    great_circle_distance_2 =   []
    
    for i in range(len(index_list)):

        lat_list.append(lat[index_list[i]])
        lng_list.append(lng[index_list[i]])
        mag_list.append(mag[index_list[i]])
        
        burst_centroid_lat = BSTUtilities.mean_val(lat_list)
        burst_centroid_lng = BSTUtilities.mean_val(lng_list)
        
    #   Compute Radius of Gyration
    
    for i in range(len(lat_list)):

        lat1 = burst_centroid_lat
        lng1 = burst_centroid_lng
        lat2 = lat_list[i]
        lng2 = lng_list[i]
        dist = BSTUtilities.compute_great_circle_distance(lat1, lng1, lat2, lng2)
            
        great_circle_distance_2.append(math.pow(dist,2))
        
    radius_gyration = math.pow(BSTUtilities.mean_val(great_circle_distance_2), 0.5)
    
    return radius_gyration, burst_centroid_lat, burst_centroid_lng

    #.................................................................
    
def reject_outliers_clusters(rgyr_cluster, year_cluster, sd_cluster):

    #   This method removes lat-lng points that are far from the median of the burst
    
    rgyr_cluster_filtered   =   []
    year_cluster_filtered   =   []
    
    median_rgyr = np.median(rgyr_cluster)
    
    print('')
    print('----------------------------------------')
    print('')
    print('Median Radius of Gyration: ', str(round(median_rgyr,2)) + ' (km)')
    print('')


    
    number_all_clusters = len(rgyr_cluster)

    for i in range(len(rgyr_cluster)):
    
        if rgyr_cluster[i] < sd_cluster * median_rgyr:
            rgyr_cluster_filtered.append(rgyr_cluster[i])
            year_cluster_filtered.append(year_cluster[i])
            
    number_reduced_clusters = len(rgyr_cluster_filtered)
    
    print('Fraction of Clusters Removed: ', str(100.0*(1 - float(number_reduced_clusters)/float(number_all_clusters))) + '%')
    print('')
    print('----------------------------------------')
    print('')

    return rgyr_cluster_filtered, year_cluster_filtered
    
    #.................................................................
    
def reject_outliers_clusters_alternate(ratio_limit, burst_index, burst_min_size, mag, lat, lng, date, time, years):

    #   This method removes clusters that are not compact
    
    mass_list           =   []
    year_list           =   []
    radgyr_list         =   []
    index_cluster       =   []
    centroid_lat_list   =   []
    centroid_lng_list   =   []
    
#     for i in range(3,len(burst_index)):        #   Iterating over the bursts/swarms
    for i in range(len(burst_index)):        #   Iterating over the bursts/swarms
    
        if len(burst_index[i]) >= burst_min_size:
        
            index_list = burst_index[i]
            lat_list    = []
            lng_list    = []

            
            mass_list.append(float(len(burst_index[i])))    #   Number of events in the burst

            year_list.append(years[burst_index[i][0]])      #   Year of the first event in the burst
            
            for j in range(len(burst_index[i])):
            
                lat_list.append(lat[index_list[j]])
                lng_list.append(lng[index_list[j]])

            radius_gyration, burst_centroid_lat, burst_centroid_lng = \
                    compute_burst_radius_gyration(burst_index[i], mag, lat, lng, date, time, years)
            
            radgyr_list.append(radius_gyration)     #   This is the time series we want to use EMA on
            
            index_cluster.append(index_list[0])
            
            centroid_lat_list.append(burst_centroid_lat)
            centroid_lng_list.append(burst_centroid_lng)
            
    density_list   =   []
    
    for i in range(len(mass_list)):
        density = mass_list[i]/((radgyr_list[i]))
        density_list.append(density)
        
    density_list = [math.log(density_list[i],10) for i in range(len(density_list))]
    
    rgyr_cluster_filtered   =   []
    year_cluster_filtered   =   []
    index_cluster_filtered  =   []
    centroid_lat_filtered   =   []
    centroid_lng_filtered   =   []
    
    for i in range(len(mass_list)):
        if density_list[i] >= ratio_limit:
            rgyr_cluster_filtered.append(radgyr_list[i])
            year_cluster_filtered.append(year_list[i])
            index_cluster_filtered.append(index_cluster[i])
            centroid_lat_filtered.append(centroid_lat_list[i])
            centroid_lng_filtered.append(centroid_lng_list[i])
        
    return rgyr_cluster_filtered, year_cluster_filtered, index_cluster_filtered, centroid_lat_filtered, centroid_lng_filtered
    
    #.................................................................
    
def calc_radgyr_EMA_time_dev(burst_min_size, completeness_mag, Location, N_Steps,\
        circle_catalog_date_start, circle_catalog_date_end, burst_index, mag, lat, lng, date, time, years, \
        ratio_limit, sd_factor, cutoff_start_year, regional_rate):
        
    #   This routine calculates a single radius of gyration time series for given ENF and CLF values
    
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
                    compute_burst_radius_gyration(burst_index[i], mag, lat, lng, date, time, years)
            
            radgyr_list.append(radius_gyration)     #   This is the time series we want to use EMA on
            year_swarm.append(years[index_list[0]])
            index_swarm.append(index_list[0])

#   Filter the clusters now, we filtered the events in the clusters previously
    
    radgyr_list, year_swarm, index_swarm, centroid_lat, centroid_lng = \
            reject_outliers_clusters_alternate(ratio_limit, burst_index, burst_min_size, mag, lat, lng, date, time, years)

    #   Independent of any start date down to here -------------------------
            
    average_rate_start_date =   circle_catalog_date_start
    average_rate_end_date   =   circle_catalog_date_end
    
    if average_rate_end_date < 1985.0:
        average_rate_end_date = 1985.0      #   The beginning of decent digital data
    
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
            
        radgyr_EMA = EMA_weighted_time_series(RG_list_raw, N_Steps)
        
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
#   Here we build a continuous, constant length, regular interval, 
#       time series for radgyr as a function of time
#
#   We are basically interpolating between points here
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
            

#     print('Length of year_swarm_regular time series', len(year_swarm_regular))
#   Here we write the values of the interpolated radii of gyration into a file.
#       Values are 1 week prior to the major earthquakes.  We will use this
#       data to optimize the weights for the time series

    write_precursor_radii_to_file(years, mag, year_swarm_regular,radgyr_list_regular, ratio_limit, sd_factor)

    return year_swarm_regular,radgyr_list_regular
    
    #.................................................................

def write_precursor_radii_to_file(years, mag, year_swarm_regular,radgyr_list_regular, ratio_limit, sd_factor):

#   ratio_limit is ENF.  sd_factor is CLF

    interpolation_file_name = 'interpolation_file.txt'
    interpolation_file = open(interpolation_file_name, 'a')
    
#   Find dates of large earthquakes

    earthquake_year_M7  =   []
    
    for i in range(len(years)):
        if float(mag[i]) >= 7.0:
            earthquake_year_M7.append(years[i])

    radgyr_M7_EQ    = [0.0 for kk in range(4)]
    year_M7_EQ      = [0.0 for kk in range(4)]
    
    kk = -1

    for ii in range(4):
        
            for jj in range(len(year_swarm_regular)):
                if year_swarm_regular[jj] < earthquake_year_M7[ii]-0.0192:     #  Use date of M>7 earthquake minus 1 week since error in
                                                                        #    year-to-day conversion accumulates over time
                    year_M7_EQ[ii]      = year_swarm_regular[jj]        #  Overwrite the data until the large EQ date is passed
                    radgyr_M7_EQ[ii]    = radgyr_list_regular[jj]

    data_list = ratio_limit, sd_factor, year_M7_EQ[0], year_M7_EQ[1],year_M7_EQ[2],year_M7_EQ[3], \
            radgyr_M7_EQ[0], radgyr_M7_EQ[1],radgyr_M7_EQ[2],radgyr_M7_EQ[3]
    print(ratio_limit, sd_factor, year_M7_EQ[0], year_M7_EQ[1],year_M7_EQ[2],year_M7_EQ[3], \
            radgyr_M7_EQ[0], radgyr_M7_EQ[1],radgyr_M7_EQ[2],radgyr_M7_EQ[3], file=interpolation_file)    

    interpolation_file.close()
    
    return

    #.................................................................
    
def weight_optimizer():

#   This code optimizes the combination of the radii of gyration curves using 
#       a brute force optimization

    interpolation_file_name = 'interpolation_file.txt'
    interpolation_file = open(interpolation_file_name, 'r')
    
    count = len(open(interpolation_file_name).readlines(  ))
    print('')
    print('length of interpolation file: ',count)
    print('')
    
    ENF_list    =   []
    CLF_list    =   []

    radii_list  =   [[] for i in range(count)]  #   A list of lists
    
    i=0
    for line in interpolation_file:

        items = line.strip().split()
#         print(items)
#         print(items[6],items[7],items[8],items[9])
        
        ENF_list.append(items[0])
        CLF_list.append(items[1])
        
        radii_list[i].extend((items[6],items[7],items[8],items[9]))
        
        i += 1
        
#     for i in range(len(radii_list)):
#         print('radii_list: ', radii_list[i])
        
    #   print initial non-normalized weights

    initial_weights =   [1.0/float(count) for i in range(count)]     #   We start with equal weighting
    
    adjusted_weights=   [1.0/float(count) for i in range(count)]
#     
# #     print('')
# #     print('Initial Non-Normalized Weights: ', initial_weights)
# #     print('')
#     
# #     l = [random.randint(0,1) for i in range(100)]
# #     print('l: ',l)
# #     
#     
#     for kk in range(len(ENF_list)):
#     
#         random_weights = random.sample(range(0, len(initial_weights)), len(initial_weights))
#         sum_of_weights = sum(random_weights)
#         random_weights = [float(random_weights[i])/sum_of_weights for i in range(len(random_weights))]
    
    mean_radius_4 = [0. for i in range(4)]
    
    for i in range(4):
        for j in range(count):
#             mean_radius[i] += random_weights[j]*float(radii_list[j][i])
            mean_radius_4[i] += initial_weights[j]*float(radii_list[j][i])
    
    mean_radius_init = np.mean(mean_radius_4)
    mean_radius_init = round(mean_radius_init,2)
    
    std_radius_init  = np.std(mean_radius_4)
    std_radius_init = round(std_radius_init,2)
    
    mean_radius_4_init = [round(mean_radius_4[i],2) for i in range(4)]
    
    trials = 20000
    
    min_mean_radius = mean_radius_init
    min_std_radius  = std_radius_init
    
    for k in range(trials):       
    
        random_weights = random.sample(range(0, count), count)
        sum_of_weights = sum(random_weights)
        random_weights = [float(random_weights[i])/sum_of_weights for i in range(len(random_weights))]
        
        mean_radius_4 = [0. for i in range(4)]
        
        for i in range(4):

            for j in range(count):
                mean_radius_4[i] += random_weights[j]*float(radii_list[j][i])
    
        mean_radius = np.mean(mean_radius_4)
        mean_radius = round(mean_radius,2)
    
        std_radius  = np.std(mean_radius_4)
        std_radius = round(std_radius,2)
        
        mean_radius_4 = [round(mean_radius_4[i],2) for i in range(4)]
#         print('Adjusted Mean_Radii (Km): ', mean_radius_4)
#         print('Mean Adjusted Radius Value: ', mean_radius, '+/-', std_radius, 'Km)')
#         print('')
        
        if (std_radius < min_std_radius):
            min_std_radius = std_radius
            min_mean_radius = mean_radius
            min_random_weights = random_weights
            
#     print(min_random_weights)


        
    print('')
    print('------------------------------------------------')
    print('')
    print('Non-adjusted Mean_Radii (Km): ', mean_radius_4_init)
    print('')
    print('Mean Non-Adjusted Radius Value: ', mean_radius_init, '+/-', std_radius_init, 'Km')
    print('')
    print('Minimum Mean Adjusted Radius Value: ', min_mean_radius, '+/-', min_std_radius, 'Km')
    print('')
    
    for i in range(4):
        values_list = []
        
        for j in range(count):
            values_list.append(float(radii_list[j][i]))
        
        mean_ts, std_ts = time_series_mean_std(values_list, min_random_weights)
        
        print('EQ: ',i+1,' Minimum Mean Adjusted Radius Value: ', round(mean_ts,2), '+/-', round(std_ts,2), 'Km')
        print('')
    
    print('------------------------------------------------')
    print('')

    
    interpolation_file.close()

#     return initial_weights and adjusted weights

    return initial_weights, adjusted_weights
    
    #.................................................................
    
def time_series_mean_std(values_list, weight_list):

#   This method computes the mean and standard deviation of the time series
#       at every time step using the adjusted time series weights

    mean_ts = 0.
    std_ts = 0.
    
    mean_ts_list = []
    
    for i in range(len(values_list)):
    
        mean_term = float(values_list[i])*float(weight_list[i])
        mean_ts_list.append( mean_term )
        
    mean_ts = sum(mean_ts_list)
    
    variance_ts_list = []
    
    for i in range(len(values_list)):
    
        variance_term = float(weight_list[i]) * math.pow(float(values_list[i])-float(mean_ts),2)
        variance_ts_list.append( variance_term )
            
    std_ts = math.pow(sum(variance_ts_list),0.5)

    return mean_ts, std_ts
    
    #.................................................................
    
