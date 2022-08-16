#!/usr/bin/env python3

# Authour: Adam Wright
# Purpose: Calculate tides at L4 so can time dives to high and low water

# Import packages
import utide
import matplotlib.dates as mdates
from scipy.signal import argrelextrema

import numpy as np
import datetime
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import astropy.utils.iers
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_moon, Angle
from astroplan import moon
import timeit
import calendar

import subprocess
from dateutil import parser
import SpectralSplit
import scipy.integrate
import os
import sys
import pdb
import warnings
warnings.simplefilter('ignore', np.RankWarning)

# Empirical function for deriving L4 offsets from Devonport
# Derived by Reg Uncles
def tidal_offset(tide_heights, tide_time):

   slack_times = np.array([], dtype=object)
   counter = 0
   for HW in tide_heights:
      # All Tides with No Wind
      # Slack before HW Devonport (hours) = 3.569 - 0.67909*(HW) + 0.0876*(HW)^2
      slack_time_offset_before = 3.569 - 0.67909*HW + 0.0876*HW**2
      slack_time = str(datetime.datetime.strptime(str(tide_time[counter]), '%Y-%m-%dT%H:%M:%S.%f000') - datetime.timedelta(hours=-slack_time_offset_before))
      slack_times = np.append(slack_times,slack_time[:-7])
      #slack_offsets = np.append(slack_offsets, slack_time_offset_before)
      # Slack after HW Devonport (hours) = 6.862 - 1.1*(HW) + 0.08264*(HW)^2
      slack_time_offset_after = 6.862 - 1.1*HW + 0.08264*HW**2
      slack_time = str(datetime.datetime.strptime(str(tide_time[counter]), '%Y-%m-%dT%H:%M:%S.%f000') - datetime.timedelta(hours=+slack_time_offset_after))
      slack_times = np.append(slack_times,slack_time[:-7])
      #slack_offsets = np.append(slack_offsets, slack_time_offset_after)
      counter += 1

   return slack_times

def format_time(time_array):
   formatted_time_array = np.array([])
   for time in time_array:
      formatted_time_array = time.strftime('%Y-%m-%d %H:%M:%S')
   return formatted_time_array

def main():

##  Create Parser instructions
    aparser = argparse.ArgumentParser()
    # TIDE
    aparser.add_argument("-t", "--tidal", action="store_true", help="Tidal model - REQUIRED: TEXT FILE CONTAINING TIDAL DATA FROM COASTAL OBSERVATORY")
    aparser.add_argument("-d", "--datum", action="store", type=float, help="Depth in metres to datum below Mean Spring tide")
    aparser.add_argument("-dp", "--datum_percentage", action="store", type=float, help="Depth percentage of tidal range for Mean Spring tide")

    aparser.add_argument("-o", "--output", action="store_true", help="generate output file(s)")

    # MODEL PARAMETERS
    aparser.add_argument("-T", "--time", action="store", type=float, help="Time increment in decimal hours i.e. 0.25 = 15 minutes")

    aparser.add_argument("-start", "--start", action="store", type=str, help="Start Date (yyyy-mm-dd)")
    aparser.add_argument("-end", "--end", action="store", type=str, help="End Date (yyyy-mm-dd)")

    args = aparser.parse_args()

    # Assign empty variables
    TL = []; t = []; waterdepth = []; tide_h = []; ftide_h = []; Zc = []
    #all_tide_times = []

    if args.start:
       start_date = args.start; print("Start date: ", start_date)
    else:
       start_date = "2001-01-01"; print("Start date: ", start_date)
    if args.end:
       end_date = args.end; print("End date: ", end_date)
    else:  
       end_date = "2001-01-14"; print("End date: ", end_date)
           
    data_start_date = subprocess.getoutput(f"date --date '{start_date}' +'%F %H:%M:%S'")
    print("Start date", data_start_date)
    data_end_date = subprocess.getoutput(f"date --date '{end_date}' +'%F %H:%M:%S'")
    print("End date", data_end_date)

    date_range_start = parser.parse(data_start_date)
    date_range_end = parser.parse(data_end_date)
    date_range_start_string_without_delta = date_range_start.strftime("%F %T")
    date_start = datetime.datetime.strptime(date_range_start_string_without_delta, '%Y-%m-%d %H:%M:%S')
    date_range_end_string_without_delta = date_range_end.strftime("%F %T")
    date_end = datetime.datetime.strptime(date_range_end_string_without_delta, '%Y-%m-%d %H:%M:%S')
    tdiff_month = date_end.month-date_start.month
    if tdiff_month == 0:
        tdiff_month = 1
    tt_e = date_end.timetuple()
    tt_e = tt_e.tm_yday
    tt_s = date_start.timetuple()
    tt_s = tt_s.tm_yday
    delta = date_end-date_start
    tdiff_day = delta.days
    year = date_start.year
    start_month = date_start.month

    geo_location = 'Plymouth_L4'
    latitude_deg = 50.27; longitude_deg = -4.13
    Tide_fname = "TidePlymouth_L1.csv"
    
    t_incr = args.time

    if (args.tidal):
        
        print('Running tidal model...')

        if args.datum_percentage:
            datum_percentage = args.datum_percentage
        else:
            datum_percentage = 0.1
        
##      Read in tidal data 
##      Tidal data from a variety of sources 
##      1. PSMSL
##         https://www.psmsl.org/data/
##      2. University of Hawaii Sea Level Centre
##         http://uhslc.soest.hawaii.edu/data/
##      3. Generate data from the TPXO data files 
##         https://tpxows.azurewebsites.net/  
##         ==> generate weekly tidal data, multiple times and then create a ?_L1.csv file from about a month's worth of data
##      Use routines found in the Tide_PrePro directory and create a suitable ?_L1.csv file
##      Put data in Model/Required/TideGauge directory

##      https://www.tide-forecast.com/ - good site for eyeballing how accurately the tide is recreated.

        print('     using tidegauge data for tidal coefficients')
        tidepath = os.getcwd() + "/Required/TideGauge/" + Tide_fname # path of data file output
        tides_ = pd.read_csv(tidepath, delimiter=',', engine='python') #usecols=np.arange(16,48), engine='python')
        df = pd.DataFrame(tides_)
        for i in range(len(tides_)):

##       Convert/alter DataFrame Variables to type datetime and float
           tide_date = df.iloc[i,0]
           tide_level =df.iloc[i,1]
           tide_level = float(tide_level)

##       Convert date to datetime format
           T = datetime.datetime.strptime(tide_date, '%Y/%m/%d %H:%M:%S') # All stations to a standard timestamp
##       append variables to lists
           TL.append(tide_level)
           t.append(T)

##      Convert time back to DataFrame
        t_df = pd.to_datetime(t)
##      Convert time to format of UTide ##'date2num' function only seems to work with pandas DF ##
        time = mdates.date2num(t_df.to_pydatetime())

        # Calculate the tidal coefficients from tide gauge data 
        tide = np.array(TL, dtype=float)
        c = utide.solve(time, u=tide, v=None, lat=latitude_deg,
                        nodal=False,
                        trend=False,
                        method='ols',
                        conf_int='linear',
                        Rayleigh_min=0.95)


        # 1. Generate list of date & time for time period of interest at 1 hour resolution
        dts = (pd.DataFrame(columns=['NULL'],index=pd.date_range(start_date, end_date,freq='1T')).between_time('00:00','23:59').index.strftime('%Y-%m-%d %H:%M:%S').tolist())
        tidetime = pd.to_datetime(dts)
        time = mdates.date2num(tidetime.to_pydatetime())
        # 2. Reconstruct the tide over that period so can determine datum more satisfactorily
        TCreconst = utide.reconstruct(time, c)
        TL = TCreconst.h

        if args.datum:
            datum = args.datum
        else:
            datum = round((min(TL) + (max(TL) - min(TL))*datum_percentage),2)
            
        print('     finished tidal model...')
        
        # now determine times of high and low water
        low_tide = argrelextrema(TL, np.less)
        high_tide = argrelextrema(TL, np.greater)
        
        low_tide_times = tidetime[low_tide]; low_tide_values = TL[low_tide]
        high_tide_times = tidetime[high_tide]; high_tide_values = TL[high_tide]
        
        # join the tide times together into an unsorted array
        all_tide_times = np.concatenate((low_tide_times, high_tide_times), axis=None)
        # join the tide level together
        all_tide_heights = np.concatenate((TL[low_tide], TL[high_tide]), axis=None)

        # find the sorted indices
        sorted_tide_time_order = np.argsort(all_tide_times)

        sorted_tide_time = np.sort(all_tide_times)
        sorted_tide_heights = all_tide_heights[sorted_tide_time_order]
        
        low_tide_label = ['Low'] *len(TL[low_tide])
        high_tide_label = ['High'] *len(TL[high_tide])
        all_tide_labels = np.concatenate((low_tide_label, high_tide_label))
        
        sorted_tide_labels = all_tide_labels[sorted_tide_time_order]
        
        # Output dataframe
        df_out = pd.DataFrame()
        df_out['Time(UTC)'] = sorted_tide_time
        df_out['Height(m)'] = sorted_tide_heights
        df_out['Hi/Lo'] = sorted_tide_labels
        
        # Calculate the tidal offset at L4
        # Only want to pass through the high water times here
        slack_times = tidal_offset(sorted_tide_heights[sorted_tide_labels == 'High'], 
                                     sorted_tide_time[sorted_tide_labels == 'High'])
        
        # issue if there are more slack waters in the given period than Hi / Lo tides
        sorted_slack_times = np.sort(slack_times)
        if (len(sorted_slack_times) > len(sorted_tide_time)):
           sorted_slack_times = sorted_slack_times[:-1]
         
        df_out['Slack(UTC)'] = pd.to_datetime(pd.Series(sorted_slack_times))
        
        if args.output:
           df_out.to_csv('Output/L4buoy_tides.txt', sep=' ', index=False, float_format='%.2f')
        
        sys.exit(0)
  
# Run script if called from command line.   
if __name__=='__main__':
    main()
    #plt.show()
  
        
       
                     





