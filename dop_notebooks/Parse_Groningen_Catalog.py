#   Parse Groningen catalog and write in standard input catalog format

import sys
import os
import numpy as np
from array import array

import datetime
import dateutil.parser

import time
from time import sleep  #   Added a pause of 30 seconds between downloads

import math

# data_file_in    = 'Groningen_Data_2000-2019.csv'
data_file_out   = 'Groningen_Data.catalog'

#print ''
name_of_catalog = raw_input('Enter name of catalog: ')
data_file_in = str(name_of_catalog)

dateSTD     =   []
timeUTC     =   []
decimalY    =   []
lat         =   []
lng         =   []
mag         =   []
dep         =   []

input_file  = open(data_file_in,'r')

number_events = 0

for line in input_file:
    items = line.strip().split(',')
     
    date = dateutil.parser.parse(items[0].split('T')[0])     # gives year-month-day (space) hours:minutes:seconds
#   date = dateutil.parser.parse(items[0])    <-- gives year-month-day only
    ts = date.year + float(date.strftime("%j"))/366
     
    print(items)
     
    raw_date_time = items[0].split('T')
    YY = raw_date_time[0].split('-')[0]
    MM = raw_date_time[0].split('-')[1]
    DD = raw_date_time[0].split('-')[2]
    date_string = YY + '/' + MM + '/' + DD
    dateSTD.append(date_string)
     
    timeUTC .append(raw_date_time[1])
    decimalY.append(ts)
    lat.append(items[2])
    lng.append(items[3])
    mag.append(items[1])
    dep.append(items[4])
    number_events += 1
     
input_file.close()
     
output_file = open(data_file_out, 'w')

for j in range(number_events):
    i = number_events - j - 1
    print >> output_file, dateSTD[i],timeUTC[i],decimalY[i],lng[i],lat[i],mag[i],dep[i]

output_file.close()
     