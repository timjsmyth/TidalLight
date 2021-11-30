#!/usr/bin/env python
# Description:
# The script is designed to plot the data output from the model: IntensityModel_Bb-DateTicks2.py
# The output data in file DATA_1-366 contains the modelled lunar irradiance at the equator and  at a latitude of 50.25 N (degrees)
# The radiative transfer model utilises a number of python packages such as UTide, pysolar and astropy for ephemeral calculations. 
# Resulting values are attenuated according to the location via Beers law, assuming a point source and bull diffusion coefficient.  

# Author: Adam Wright
# Date: 21/10/2020
# Institution: Plymouth Marine Laboratory

import sys
import pandas as pd
import numpy as np
import datetime
import matplotlib.dates as mdates

import matplotlib.pyplot as plt
import pandas as pd
import os

import timeit
import calendar

import subprocess
from dateutil import parser

def time_me():
    timeline = input("Do you wish to enter a date range [Y/n]:\n"
                    "Y: Enter date range for model\n"
                    "n: Default to synodic month: 6th June- 7th July 2020 \n")

    if timeline == "Y":
        # Ask the user to input date range
        start_date_input = input("Please enter start date: ")
        end_date_input = input("Please enter end date: ")
    else:
        # Otherwise just run synodic month
        dt_start = "2020-06-06" # 6th of June 2020
        dt_end = "2020-07-07" # 7th of July 2020
        start_date_input = datetime.datetime.strptime(dt_start, "%Y-%m-%d") # format string  as a datetime object
        end_date_input = datetime.datetime.strptime(dt_end, "%Y-%m-%d") # format string  as a datetime object

    data_start_date = subprocess.getoutput(f"date --date '{start_date_input}' +'%F %H:%M:%S'")
    data_end_date = subprocess.getoutput(f"date --date '{end_date_input}' +'%F %H:%M:%S'")

    date_range_start = parser.parse(data_start_date)
    date_range_end = parser.parse(data_end_date)
    date_range_start_string_without_delta = date_range_start.strftime("%F %T")
    date_start = datetime.datetime.strptime(date_range_start_string_without_delta, '%Y-%m-%d %H:%M:%S') # Output start date as datetime object
    date_range_end_string_without_delta = date_range_end.strftime("%F %T")
    date_end = datetime.datetime.strptime(date_range_end_string_without_delta, '%Y-%m-%d %H:%M:%S') # Output end date as datetime object
    tdiff_month = date_end.month-date_start.month # Find difference between dates in number of months
    if tdiff_month == 0:
        tdiff_month = 1
    tt_e = date_end.timetuple() # day number of End
    tt_e = tt_e.tm_yday
    tt_s = date_start.timetuple() # day number of Start
    tt_s = tt_s.tm_yday
    tdiff_day = tt_e-tt_s
    year = date_start.year # Year of interest
    return tt_s, tt_e, year, date_end, date_start

    
def Figure2_plots(): 
    figure_name1 = ('Fig2A-B.png')
    figure_name2 = ('Fig2C-F.png')

    # Set date interval
    tt_s, tt_e, year, enddate, startdate = time_me()
    TT_S = startdate.timetuple()
    TT_E = enddate.timetuple()
    
    # Set datum:
    datum_set_choice = input("Would you like to set your own datum? \n Y: Input datum depth in metres below sea level \n n: Default to 1.7m \n")
    if datum_set_choice=="Y": 
        datum_set = float(input("Input datum: \n"))
    else:
        datum_set = 1.7
    df = pd.read_csv(sys.argv[1]) # input 50N dataframe (same data frame for 2C-F)
    EQ_df = pd.read_csv(sys.argv[2]) # input EQ dataframe
    df['date'] = pd.to_datetime(df['date'])
    df['date'].tolist()
    df['Lunar_Int_SeaLevel_EQ(PPFD)'] = EQ_df['Lunar_Int_SeaLevel_EQ(PPFD)'] # assign varibale explicit values from secondary dataset.

    #########################################################################
    #######                          2A-B                             #######
    #########################################################################

    # This plot is figure 2A-B in Marine light pollution research guide, 2020
    # Plot the full year of lunar data for lat = 50 and lat = 0.1 equator
    fig, ax = plt.subplots(2, figsize=(13,7.34))
    months = mdates.MonthLocator()
    months_fmt = mdates.DateFormatter('%B')
    
    ax[0].plot(df['decimal_day'], df['Lunar_Int_SeaLevel(PPFD)'], color='black')
    ax[0].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax[0].set_title('Lunar Irradiance at latitude 50$^\circ$N')
    ax[0].xaxis.set_major_locator(months)
    ax[0].xaxis.set_major_formatter(months_fmt)
    ax[0].set_xticklabels([])

    ax[1].plot(df['decimal_day'], df['Lunar_Int_SeaLevel_EQ(PPFD)'], color='black', alpha=0.7)
    ax[1].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax[1].set_title('Lunar Irradiance at the equator')
    ax[1].xaxis.set_major_locator(months)
    ax[1].xaxis.set_major_formatter(months_fmt)

    ax[1].set_xlabel('Month')
    ax[1].xaxis.set_label_position('bottom')
    ax[1].xaxis.set_ticks_position('bottom')

    plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.3)
    plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=25, fontsize=8)
    fig.savefig(figure_name1, dpi=450)

    #########################################################################
    #######                          2C-F                             #######
    #########################################################################

    # This plot is figure 2C-F in Marine light pollution research guide, 2020
    # Plot the synodic month of Jun-July 2020 for lunar irradiance at sea level, datum, tide and lunar phase 
    # Plymouth: Latitude 50.25 N, Longitude 4.13 W (degrees)

    fig, ax = plt.subplots(4, figsize=(13,8))
    loc = mdates.DayLocator(bymonthday=None, interval=3, tz=None)
    loc_minor = mdates.DayLocator(bymonthday=None, interval=1, tz=None)
    loc_fmt = mdates.DateFormatter('%d-%B')
    AA = np.ones(len(df['decimal_day']), dtype=int)
    aa = 2300*AA # extend day-night colouring arbitrarily above highest value
    bb = -300*AA

    # LUNAR INT SEALEVEL
    ax0i = ax[0].twiny()
    ax0 = ax[0].twinx()

    ax0.fill_between(df['decimal_day'], aa, bb, where= df['binary_day'] < AA, facecolor='grey', alpha=0.2)
    ax0.axes.get_yaxis().set_visible(False)
    ax0.set_xticklabels([]); ax0.set_xticks([])

    ax0i.plot(df['date'], df['Lunar_Int_SeaLevel(PPFD)'], color='none')
    ax0i.set_xlim([startdate, enddate])
    ax0i.xaxis.set_major_locator(loc)
    ax0i.xaxis.set_major_formatter(loc_fmt)
    ax0i.xaxis.set_label_position('bottom')
    ax0i.xaxis.set_ticks_position('bottom')
    ax0i.xaxis.set_minor_locator(loc_minor)
    ax0i.set_xticklabels([])

    ax[0].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax[0].set_title('Lunar irradiance at sea level')
    ax[0].plot(df['decimal_day'], df['Lunar_Int_SeaLevel(PPFD)'], label='Lunar irradiance at sea level', color='sienna')
    ax[0].set_ylim([(-max(df['Lunar_Int_SeaLevel(PPFD)'])/10), (max(df['Lunar_Int_SeaLevel(PPFD)'])+max(df['Lunar_Int_SeaLevel(PPFD)'])/10)])
    ax[0].set_xlim([TT_S.tm_yday, TT_E.tm_yday])

    # LUNAR INT SEABED 
    ax1 = ax[1].twinx()
    ax1i = ax[1].twiny()

    ax1.fill_between(df['decimal_day'], aa, bb, where= df['binary_day'] < AA, facecolor='grey', alpha=0.2)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_xticklabels([]); ax1.set_xticks([])

    ax[1].plot(df['decimal_day'], df['Lunar_Int_BelowTide(PPFD)'], label='Lunar irradiance at seabed', color='lightsalmon')
    ax[1].set_title('Lunar irradiance at datum ')
    ax[1].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
    ax[1].set_xlim([TT_S.tm_yday, TT_E.tm_yday])
    ax[1].set_ylim([(-max(df['Lunar_Int_BelowTide(PPFD)'])/10), (max(df['Lunar_Int_BelowTide(PPFD)'])+max(df['Lunar_Int_BelowTide(PPFD)'])/10)])

    ax1i.plot(df['date'], df['Lunar_Int_BelowTide(PPFD)'], color='none', alpha=0.1)
    ax1i.set_xlim([startdate, enddate])
    ax1i.xaxis.set_major_locator(loc)
    ax1i.xaxis.set_major_formatter(loc_fmt)
    ax1i.xaxis.set_label_position('bottom')
    ax1i.xaxis.set_ticks_position('bottom')
    ax1i.xaxis.set_minor_locator(loc_minor)
    ax1i.set_xticklabels([])

    # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
    ax2i = ax[2].twinx()
    ax2 = ax[2].twinx()
    ax2j = ax[2].twiny()

    ax2i.fill_between(df['decimal_day'], aa, bb, where= df['binary_day'] < AA, facecolor='grey', alpha=0.2)
    ax2i.axes.get_yaxis().set_visible(False)

    ax[2].plot(df['decimal_day'], df['water_height_above_datum(m)'], label='water depth', color='royalblue')
    ax[2].set_ylabel(f'Water column\n above {datum_set}m datum')
    ax[2].set_xlim([TT_S.tm_yday, TT_E.tm_yday])
    ax[2].set_ylim([(-max(df['tidal_range'])/10), (max(df['tidal_range'])+max(df['tidal_range'])/10)])
    ax[2].set_title('Height of water column') 

    ax2.axhline(y=datum_set, xmin=0, xmax=366, color='cadetblue', alpha = 0.5, linestyle='dashed', label='datum')
    ax2.legend(prop={"size":8}, loc='upper right')
    ax2.plot(df['decimal_day'], df['tidal_range'], label='tide', color='cadetblue', alpha=0.5)
    ax2.set_ylabel('Tidal range (m)', color='cadetblue')
    ax2.set_xticklabels([]); ax2.set_xticks([])
    
    ax2j.plot(df['date'], df['binary_day'], color='none')
    ax2j.set_xlim([startdate, enddate])
    ax2j.xaxis.set_major_locator(loc)
    ax2j.xaxis.set_major_formatter(loc_fmt)
    ax2j.xaxis.set_label_position('bottom')
    ax2j.xaxis.set_ticks_position('bottom')
    ax2j.xaxis.set_minor_locator(loc_minor)
    ax2j.set_xticklabels([])

    # LUNAR PHASE
    ax[3].plot(df['decimal_day'], df['phase(normalised)'], label='Lunar phase', color='darkgrey')
    ax[3].set_title('Lunar phase cycle')
    ax[3].set_ylabel('Normalised phase')
    ax[3].set_xlim([int(TT_S.tm_yday), int(TT_E.tm_yday)])
    ax[3].axes.get_xaxis().set_visible(False)
    
    ax3i = ax[3].twiny()
    ax3 = ax[3].twinx()
    ax3.plot(df['date'], df['phase_label'], color='none', alpha=0.1)
    ax3.axes.get_xaxis().set_visible(False)

    ax3i.plot(df['date'], df['phase(normalised)'], color='none', alpha=0.1)
    ax3i.xaxis.set_major_locator(loc)
    ax3i.xaxis.set_major_formatter(loc_fmt)
    ax3i.xaxis.set_minor_locator(loc_minor)
    ax3i.set_xlabel("Date")
    ax3i.xaxis.set_label_position('bottom')
    ax3i.xaxis.set_ticks_position('bottom')
    ax3i.set_xlim([startdate, enddate])

    plt.tight_layout() 
    plt.setp(ax3i.xaxis.get_majorticklabels(), rotation=25, fontsize=8)
    fig.savefig(figure_name2, dpi=450)
    plt.show()
    
if __name__ == '__main__':
    if len(sys.argv)< 1: 
        print('Enter data filename of dataset / ensure the datatset is in the current directory')
        
    else:
        print("Working on it...")
        Figure2_plots()

    