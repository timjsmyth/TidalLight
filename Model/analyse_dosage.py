#!/bin/python3
import os
import glob
import numpy as np
import pandas as pd
import pdb
import calendar
import argparse
import re
import matplotlib.pyplot as plt
from braceexpand import braceexpand
import matplotlib.ticker as ticker
import sys

uW_W = 1.0e-6 # scaling from uW to W
indir = "/users/rsg/tjsm/ALAN/tidal/TidalLight/Model/Output/annual"
outdir = "/users/rsg/tjsm/ALAN/tidal/TidalLight/Model/Output/annual/analysis"

# Function: dosage_opf
# Description: Function to read in the ALAN, Solar and Lunar files and create
#              monthly statistics of the different dosages in separate csv files
def dosage_opf(location):
   print("Location: ", location)

   for ALAN_fname in glob.glob(indir+'/ALAN_'+location+'_clear_'+'*.csv'):
      print(ALAN_fname)
      df_ALAN_irradiances = pd.read_csv(ALAN_fname)
      # find the Lunar equivalent
      Lunar_fname = ALAN_fname.replace("ALAN_","Lunar_")
      Lunar_fname = Lunar_fname.replace("_clear_","_")
      print(Lunar_fname)
      df_lunar_irradiances = pd.read_csv(Lunar_fname)
      # find the Solar equivalent
      Solar_fname = ALAN_fname.replace("ALAN_","Solar_")
      Solar_fname = Solar_fname.replace("_clear_","_")
      print(Solar_fname)
      df_solar_irradiances = pd.read_csv(Solar_fname)
      
      # Extract the year
      yyyy = pd.DatetimeIndex(df_solar_irradiances['date']).year[0]
      
      # Integrate over different time periods
      solar_t_sec = df_solar_irradiances['time_increment(hr)'][0]*60*60
      lunar_t_sec = df_lunar_irradiances['time_increment(hr)'][0]*60*60
      ALAN_t_sec = df_ALAN_irradiances['time_increment(hr)'][0]*60*60
      
      general_t_sec = solar_t_sec

      # Broadband surface
      solar_I_BB = df_solar_irradiances['I_BB(W/m2)']*solar_t_sec
      lunar_I_BB = df_lunar_irradiances['I_BB(uW/m2)']*lunar_t_sec*uW_W
      ALAN_I_BB = df_ALAN_irradiances['ALAN_BB(uW/m^2)']*ALAN_t_sec*uW_W
      # Spectral surface
      solar_I_Red = df_solar_irradiances['I_Red(W/m2)']*solar_t_sec
      solar_I_Green = df_solar_irradiances['I_Green(W/m2)']*solar_t_sec
      solar_I_Blue = df_solar_irradiances['I_Blue(W/m2)']*solar_t_sec
      lunar_I_Red = df_lunar_irradiances['I_Red(uW/m2)']*lunar_t_sec*uW_W
      lunar_I_Green = df_lunar_irradiances['I_Green(uW/m2)']*lunar_t_sec*uW_W
      lunar_I_Blue = df_lunar_irradiances['I_Blue(uW/m2)']*lunar_t_sec*uW_W
      ALAN_I_Red = df_ALAN_irradiances['I_Red(uW/m2)']*ALAN_t_sec*uW_W
      ALAN_I_Green = df_ALAN_irradiances['I_Green(uW/m2)']*ALAN_t_sec*uW_W
      ALAN_I_Blue = df_ALAN_irradiances['I_Blue(uW/m2)']*ALAN_t_sec*uW_W
      
      # Broadband at datum
      solar_I_BB_datum = df_solar_irradiances['I_BB_Datum(W/m2)']*solar_t_sec
      lunar_I_BB_datum = df_lunar_irradiances['I_BB_datum(uW/m2)']*lunar_t_sec*uW_W
      ALAN_I_BB_datum = df_ALAN_irradiances['ALAN_BB_datum(uW/m^2)']*ALAN_t_sec*uW_W
      # Spectral at datum
      solar_I_Red_datum = df_solar_irradiances['I_Red_datum(W/m2)']*solar_t_sec
      solar_I_Green_datum = df_solar_irradiances['I_Green_datum(W/m2)']*solar_t_sec
      solar_I_Blue_datum = df_solar_irradiances['I_Blue_datum(W/m2)']*solar_t_sec
      lunar_I_Red_datum = df_lunar_irradiances['I_Red_datum(uW/m2)']*lunar_t_sec*uW_W
      lunar_I_Green_datum = df_lunar_irradiances['I_Green_datum(uW/m2)']*lunar_t_sec*uW_W
      lunar_I_Blue_datum = df_lunar_irradiances['I_Blue_datum(uW/m2)']*lunar_t_sec*uW_W
      ALAN_I_Red_datum = df_ALAN_irradiances['I_Red_datum(uW/m2)']*ALAN_t_sec*uW_W
      ALAN_I_Green_datum = df_ALAN_irradiances['I_Green_datum(uW/m2)']*ALAN_t_sec*uW_W
      ALAN_I_Blue_datum = df_ALAN_irradiances['I_Blue_datum(uW/m2)']*ALAN_t_sec*uW_W

      # Broadband surface
      monthly_solar_I_BB = [] ; monthly_lunar_I_BB = [] ; monthly_ALAN_I_BB = []
      monthly_twilight_I_BB = [];
      # Spectral surface
      monthly_solar_I_Red = [] ; monthly_solar_I_Green = [] ; monthly_solar_I_Blue = [] ;
      monthly_lunar_I_Red = [] ; monthly_lunar_I_Green = [] ; monthly_lunar_I_Blue = [] ;
      monthly_ALAN_I_Red = [] ; monthly_ALAN_I_Green = [] ; monthly_ALAN_I_Blue = [] ; 
      monthly_twilight_I_Red = [] ; monthly_twilight_I_Green = [] ; monthly_twilight_I_Blue = [] ;
      
      # Broadband at datum
      monthly_solar_I_BB_datum = [] ; monthly_lunar_I_BB_datum = [] ; monthly_ALAN_I_BB_datum = [] ;
      monthly_twilight_I_BB_datum = []
      # Spectral at datum
      monthly_solar_I_Red_datum = [] ; monthly_solar_I_Green_datum = [] ; monthly_solar_I_Blue_datum = [] ;
      monthly_lunar_I_Red_datum = [] ; monthly_lunar_I_Green_datum = [] ; monthly_lunar_I_Blue_datum = [] ;
      monthly_ALAN_I_Red_datum = [] ;  monthly_ALAN_I_Green_datum = [] ; monthly_ALAN_I_Blue_datum = [] ;
      monthly_twilight_I_Red_datum = []; monthly_twilight_I_Green_datum = [] ; monthly_twilight_I_Blue_datum = [] ;

      # Broadband surface
      monthly_max_solar_I_BB = [] ; monthly_max_lunar_I_BB = [] ; monthly_max_ALAN_I_BB = []
      monthly_max_twilight_I_BB = [];
      # Spectral surface
      monthly_max_solar_I_Red = [] ; monthly_max_solar_I_Green = [] ; monthly_max_solar_I_Blue = [] ;
      monthly_max_lunar_I_Red = [] ; monthly_max_lunar_I_Green = [] ; monthly_max_lunar_I_Blue = [] ;
      monthly_max_ALAN_I_Red = [] ; monthly_max_ALAN_I_Green = [] ; monthly_max_ALAN_I_Blue = [] ; 
      monthly_max_twilight_I_Red = [] ; monthly_max_twilight_I_Green = [] ; monthly_max_twilight_I_Blue = [] ;
      
      # Broadband at datum
      monthly_max_solar_I_BB_datum = [] ; monthly_max_lunar_I_BB_datum = [] ; monthly_max_ALAN_I_BB_datum = [] ;
      monthly_max_twilight_I_BB_datum = []
      # Spectral at datum
      monthly_max_solar_I_Red_datum = [] ; monthly_max_solar_I_Green_datum = [] ; monthly_max_solar_I_Blue_datum = [] ;
      monthly_max_lunar_I_Red_datum = [] ; monthly_max_lunar_I_Green_datum = [] ; monthly_max_lunar_I_Blue_datum = [] ;
      monthly_max_ALAN_I_Red_datum = [] ;  monthly_max_ALAN_I_Green_datum = [] ; monthly_max_ALAN_I_Blue_datum = [] ;
      monthly_max_twilight_I_Red_datum = []; monthly_max_twilight_I_Green_datum = [] ; monthly_max_twilight_I_Blue_datum = [] ;

      ## Monthly
      mm = [1]
      mm_str = []
      for m in range(1,13):
         mm.append(calendar.monthrange(yyyy, m)[1])
         mm_str.append(calendar.month_abbr[m])
      mm_str.append('Tot')

      #last_days_of month 
      mm = np.cumsum(mm)
      for m in range(0,13):
         # Final run through is total of year
         if (m == 12):
            print("Summing between days: ", mm[0],mm[12])
            mm_index = np.where((df_solar_irradiances['Jday(decimal)'] >= mm[0]) & (df_solar_irradiances['Jday(decimal)'] <= mm[12]))
            mm_twilight_index = np.where((df_solar_irradiances['Jday(decimal)'] >= mm[0]) & (df_solar_irradiances['Jday(decimal)'] <= mm[12]) & (df_solar_irradiances['night(0)_day(1)'] > 0.) & (df_solar_irradiances['night(0)_day(1)'] < 1.0))
         else:
            print("Summing between days: ", mm[m],mm[m+1]-1)
            mm_index = np.where((df_solar_irradiances['Jday(decimal)'] >= mm[m]) & (df_solar_irradiances['Jday(decimal)'] <= mm[m+1])) 
            mm_twilight_index = np.where((df_solar_irradiances['Jday(decimal)'] >= mm[m]) & (df_solar_irradiances['Jday(decimal)'] <= mm[m+1]) & (df_solar_irradiances['night(0)_day(1)'] > 0.) & (df_solar_irradiances['night(0)_day(1)'] < 1.0))

         ### Monthly totals ###
         ## Broadband ##
         # Surface #
         monthly_solar_I_BB.append(np.sum(solar_I_BB[mm_index[0]]))
         monthly_lunar_I_BB.append(np.sum(lunar_I_BB[mm_index[0]]))
         monthly_ALAN_I_BB.append(np.sum(ALAN_I_BB[mm_index[0]]))
         monthly_twilight_I_BB.append(np.sum(solar_I_BB[mm_twilight_index[0]]))
         
         # Datum #
         monthly_solar_I_BB_datum.append(np.sum(solar_I_BB_datum[mm_index[0]]))
         monthly_lunar_I_BB_datum.append(np.sum(lunar_I_BB_datum[mm_index[0]]))
         monthly_ALAN_I_BB_datum.append(np.sum(ALAN_I_BB_datum[mm_index[0]]))
         monthly_twilight_I_BB_datum.append(np.sum(solar_I_BB_datum[mm_twilight_index[0]]))
         
         ## Spectral ##
         # Surface #
         monthly_solar_I_Red.append(np.sum(solar_I_Red[mm_index[0]]))
         monthly_solar_I_Green.append(np.sum(solar_I_Green[mm_index[0]])) 
         monthly_solar_I_Blue.append(np.sum(solar_I_Blue[mm_index[0]]))
         monthly_lunar_I_Red.append(np.sum(lunar_I_Red[mm_index[0]]))
         monthly_lunar_I_Green.append(np.sum(lunar_I_Green[mm_index[0]])) 
         monthly_lunar_I_Blue.append(np.sum(lunar_I_Blue[mm_index[0]])) 
         monthly_ALAN_I_Red.append(np.sum(ALAN_I_Red[mm_index[0]])) 
         monthly_ALAN_I_Green.append(np.sum(ALAN_I_Green[mm_index[0]])) 
         monthly_ALAN_I_Blue.append(np.sum(ALAN_I_Blue[mm_index[0]])) 
         
         monthly_twilight_I_Red.append(np.sum(solar_I_Red[mm_twilight_index[0]]))
         monthly_twilight_I_Green.append(np.sum(solar_I_Green[mm_twilight_index[0]])) 
         monthly_twilight_I_Blue.append(np.sum(solar_I_Blue[mm_twilight_index[0]]))
         
         # Datum #
         monthly_solar_I_Red_datum.append(np.sum(solar_I_Red_datum[mm_index[0]]))
         monthly_solar_I_Green_datum.append(np.sum(solar_I_Green_datum[mm_index[0]])) 
         monthly_solar_I_Blue_datum.append(np.sum(solar_I_Blue_datum[mm_index[0]]))
         monthly_lunar_I_Red_datum.append(np.sum(lunar_I_Red_datum[mm_index[0]]))
         monthly_lunar_I_Green_datum.append(np.sum(lunar_I_Green_datum[mm_index[0]])) 
         monthly_lunar_I_Blue_datum.append(np.sum(lunar_I_Blue_datum[mm_index[0]])) 
         monthly_ALAN_I_Red_datum.append(np.sum(ALAN_I_Red_datum[mm_index[0]])) 
         monthly_ALAN_I_Green_datum.append(np.sum(ALAN_I_Green_datum[mm_index[0]])) 
         monthly_ALAN_I_Blue_datum.append(np.sum(ALAN_I_Blue_datum[mm_index[0]])) 
         monthly_twilight_I_Red_datum.append(np.sum(solar_I_Red_datum[mm_twilight_index[0]]))
         monthly_twilight_I_Green_datum.append(np.sum(solar_I_Green_datum[mm_twilight_index[0]])) 
         monthly_twilight_I_Blue_datum.append(np.sum(solar_I_Blue_datum[mm_twilight_index[0]]))

         ### Monthly maximums ###
         ## Broadband ##
         # Surface #
         monthly_max_solar_I_BB.append(np.max(solar_I_BB[mm_index[0]])/general_t_sec)
         monthly_max_lunar_I_BB.append(np.max(lunar_I_BB[mm_index[0]])/general_t_sec)
         monthly_max_ALAN_I_BB.append(np.max(ALAN_I_BB[mm_index[0]])/general_t_sec)
         monthly_max_twilight_I_BB.append(np.max(solar_I_BB[mm_twilight_index[0]])/general_t_sec)
         
         # Datum #
         monthly_max_solar_I_BB_datum.append(np.max(solar_I_BB_datum[mm_index[0]])/general_t_sec)
         monthly_max_lunar_I_BB_datum.append(np.max(lunar_I_BB_datum[mm_index[0]])/general_t_sec)
         monthly_max_ALAN_I_BB_datum.append(np.max(ALAN_I_BB_datum[mm_index[0]])/general_t_sec)
         monthly_max_twilight_I_BB_datum.append(np.max(solar_I_BB_datum[mm_twilight_index[0]])/general_t_sec)
         
         ## Spectral ##
         # Surface #
         monthly_max_solar_I_Red.append(np.max(solar_I_Red[mm_index[0]])/general_t_sec)
         monthly_max_solar_I_Green.append(np.max(solar_I_Green[mm_index[0]])/general_t_sec) 
         monthly_max_solar_I_Blue.append(np.max(solar_I_Blue[mm_index[0]])/general_t_sec)
         monthly_max_lunar_I_Red.append(np.max(lunar_I_Red[mm_index[0]])/general_t_sec)
         monthly_max_lunar_I_Green.append(np.max(lunar_I_Green[mm_index[0]])/general_t_sec) 
         monthly_max_lunar_I_Blue.append(np.max(lunar_I_Blue[mm_index[0]])/general_t_sec) 
         monthly_max_ALAN_I_Red.append(np.max(ALAN_I_Red[mm_index[0]])/general_t_sec) 
         monthly_max_ALAN_I_Green.append(np.max(ALAN_I_Green[mm_index[0]])/general_t_sec) 
         monthly_max_ALAN_I_Blue.append(np.max(ALAN_I_Blue[mm_index[0]])/general_t_sec) 
         
         monthly_max_twilight_I_Red.append(np.max(solar_I_Red[mm_twilight_index[0]])/general_t_sec)
         monthly_max_twilight_I_Green.append(np.max(solar_I_Green[mm_twilight_index[0]])/general_t_sec) 
         monthly_max_twilight_I_Blue.append(np.max(solar_I_Blue[mm_twilight_index[0]])/general_t_sec)
         
         # Datum #
         monthly_max_solar_I_Red_datum.append(np.max(solar_I_Red_datum[mm_index[0]])/general_t_sec)
         monthly_max_solar_I_Green_datum.append(np.max(solar_I_Green_datum[mm_index[0]])/general_t_sec) 
         monthly_max_solar_I_Blue_datum.append(np.max(solar_I_Blue_datum[mm_index[0]])/general_t_sec)
         monthly_max_lunar_I_Red_datum.append(np.max(lunar_I_Red_datum[mm_index[0]])/general_t_sec)
         monthly_max_lunar_I_Green_datum.append(np.max(lunar_I_Green_datum[mm_index[0]])/general_t_sec) 
         monthly_max_lunar_I_Blue_datum.append(np.max(lunar_I_Blue_datum[mm_index[0]])/general_t_sec) 
         monthly_max_ALAN_I_Red_datum.append(np.max(ALAN_I_Red_datum[mm_index[0]])/general_t_sec) 
         monthly_max_ALAN_I_Green_datum.append(np.max(ALAN_I_Green_datum[mm_index[0]])/general_t_sec) 
         monthly_max_ALAN_I_Blue_datum.append(np.max(ALAN_I_Blue_datum[mm_index[0]])/general_t_sec) 
         monthly_max_twilight_I_Red_datum.append(np.max(solar_I_Red_datum[mm_twilight_index[0]])/general_t_sec)
         monthly_max_twilight_I_Green_datum.append(np.max(solar_I_Green_datum[mm_twilight_index[0]])/general_t_sec) 
         monthly_max_twilight_I_Blue_datum.append(np.max(solar_I_Blue_datum[mm_twilight_index[0]])/general_t_sec)

      # Output 
      outfname = outdir+'/Dosage_'+location+'_'+str(yyyy)+'.csv'
      print("** Output file:",outfname)
      outdf = pd.DataFrame()

      outdf['Month'] = mm_str
      outdf['Year'] = pd.DatetimeIndex(df_solar_irradiances['date']).year[0]
      outdf['Location'] = location
      # Broadband surface
      outdf['solar_I_BB(J/m2)'] = monthly_solar_I_BB 
      outdf['monthly_lunar_I_BB(J/m2)'] = monthly_lunar_I_BB
      outdf['monthly_ALAN_I_BB(J/m2)'] = monthly_ALAN_I_BB
      outdf['monthly_twilight_I_BB(J/m2)'] = monthly_twilight_I_BB
      
      # Spectral surface
      outdf['monthly_solar_I_Red(J/m2)'] = monthly_solar_I_Red 
      outdf['monthly_solar_I_Green(J/m2)'] = monthly_solar_I_Green
      outdf['monthly_solar_I_Blue(J/m2)'] = monthly_solar_I_Blue 
      outdf['monthly_lunar_I_Red(J/m2)'] = monthly_lunar_I_Red 
      outdf['monthly_lunar_I_Green(J/m2)'] = monthly_lunar_I_Green 
      outdf['monthly_lunar_I_Blue(J/m2)'] = monthly_lunar_I_Blue 
      outdf['monthly_ALAN_I_Red(J/m2)'] = monthly_ALAN_I_Red 
      outdf['monthly_ALAN_I_Green(J/m2)'] = monthly_ALAN_I_Green 
      outdf['monthly_ALAN_I_Blue(J/m2)'] = monthly_ALAN_I_Blue  
      outdf['monthly_twilight_I_Red(J/m2)'] = monthly_twilight_I_Red
      outdf['monthly_twilight_I_Green(J/m2)'] = monthly_twilight_I_Green 
      outdf['monthly_twilight_I_Blue(J/m2)'] = monthly_twilight_I_Blue 
      
      # Broadband at datum
      outdf['monthly_solar_I_BB_datum(J/m2)'] = monthly_solar_I_BB_datum 
      outdf['monthly_lunar_I_BB_datum(J/m2)'] = monthly_lunar_I_BB_datum 
      outdf['monthly_ALAN_I_BB_datum(J/m2)'] = monthly_ALAN_I_BB_datum
      outdf['monthly_twilight_I_BB_datum(J/m2)'] = monthly_twilight_I_BB_datum
      # Spectral at datum
      outdf['monthly_solar_I_Red_datum(J/m2)'] = monthly_solar_I_Red_datum 
      outdf['monthly_solar_I_Green_datum(J/m2)'] = monthly_solar_I_Green_datum 
      outdf['monthly_solar_I_Blue_datum(J/m2)'] = monthly_solar_I_Blue_datum
      outdf['monthly_lunar_I_Red_datum(J/m2)'] = monthly_lunar_I_Red_datum 
      outdf['monthly_lunar_I_Green_datum(J/m2)'] = monthly_lunar_I_Green_datum 
      outdf['monthly_lunar_I_Blue_datum(J/m2)'] = monthly_lunar_I_Blue_datum 
      outdf['monthly_ALAN_I_Red_datum(J/m2)'] = monthly_ALAN_I_Red_datum 
      outdf['monthly_ALAN_I_Green_datum(J/m2)'] = monthly_ALAN_I_Green_datum
      outdf['monthly_ALAN_I_Blue_datum(J/m2)'] = monthly_ALAN_I_Blue_datum 
      outdf['monthly_twilight_I_Red_datum(J/m2)'] = monthly_twilight_I_Red_datum
      outdf['monthly_twilight_I_Green_datum(J/m2)'] = monthly_twilight_I_Green_datum 
      outdf['monthly_twilight_I_Blue_datum(J/m2)'] = monthly_twilight_I_Blue_datum 

      # Broadband surface
      outdf['solar_I_BB(W/m2)'] = monthly_max_solar_I_BB 
      outdf['monthly_max_lunar_I_BB(W/m2)'] = monthly_max_lunar_I_BB
      outdf['monthly_max_ALAN_I_BB(W/m2)'] = monthly_max_ALAN_I_BB
      outdf['monthly_max_twilight_I_BB(W/m2)'] = monthly_max_twilight_I_BB
      
      # Spectral surface
      outdf['monthly_max_solar_I_Red(W/m2)'] = monthly_max_solar_I_Red 
      outdf['monthly_max_solar_I_Green(W/m2)'] = monthly_max_solar_I_Green
      outdf['monthly_max_solar_I_Blue(W/m2)'] = monthly_max_solar_I_Blue 
      outdf['monthly_max_lunar_I_Red(W/m2)'] = monthly_max_lunar_I_Red 
      outdf['monthly_max_lunar_I_Green(W/m2)'] = monthly_max_lunar_I_Green 
      outdf['monthly_max_lunar_I_Blue(W/m2)'] = monthly_max_lunar_I_Blue 
      outdf['monthly_max_ALAN_I_Red(W/m2)'] = monthly_max_ALAN_I_Red 
      outdf['monthly_max_ALAN_I_Green(W/m2)'] = monthly_max_ALAN_I_Green 
      outdf['monthly_max_ALAN_I_Blue(W/m2)'] = monthly_max_ALAN_I_Blue  
      outdf['monthly_max_twilight_I_Red(W/m2)'] = monthly_max_twilight_I_Red
      outdf['monthly_max_twilight_I_Green(W/m2)'] = monthly_max_twilight_I_Green 
      outdf['monthly_max_twilight_I_Blue(W/m2)'] = monthly_max_twilight_I_Blue 
      
      # Broadband at datum
      outdf['monthly_max_solar_I_BB_datum(W/m2)'] = monthly_max_solar_I_BB_datum 
      outdf['monthly_max_lunar_I_BB_datum(W/m2)'] = monthly_max_lunar_I_BB_datum 
      outdf['monthly_max_ALAN_I_BB_datum(W/m2)'] = monthly_max_ALAN_I_BB_datum
      outdf['monthly_max_twilight_I_BB_datum(W/m2)'] = monthly_max_twilight_I_BB_datum
      # Spectral at datum
      outdf['monthly_max_solar_I_Red_datum(W/m2)'] = monthly_max_solar_I_Red_datum 
      outdf['monthly_max_solar_I_Green_datum(W/m2)'] = monthly_max_solar_I_Green_datum 
      outdf['monthly_max_solar_I_Blue_datum(W/m2)'] = monthly_max_solar_I_Blue_datum
      outdf['monthly_max_lunar_I_Red_datum(W/m2)'] = monthly_max_lunar_I_Red_datum 
      outdf['monthly_max_lunar_I_Green_datum(W/m2)'] = monthly_max_lunar_I_Green_datum 
      outdf['monthly_max_lunar_I_Blue_datum(W/m2)'] = monthly_max_lunar_I_Blue_datum 
      outdf['monthly_max_ALAN_I_Red_datum(W/m2)'] = monthly_max_ALAN_I_Red_datum 
      outdf['monthly_max_ALAN_I_Green_datum(W/m2)'] = monthly_max_ALAN_I_Green_datum
      outdf['monthly_max_ALAN_I_Blue_datum(W/m2)'] = monthly_max_ALAN_I_Blue_datum 
      outdf['monthly_max_twilight_I_Red_datum(W/m2)'] = monthly_max_twilight_I_Red_datum
      outdf['monthly_max_twilight_I_Green_datum(W/m2)'] = monthly_max_twilight_I_Green_datum 
      outdf['monthly_max_twilight_I_Blue_datum(W/m2)'] = monthly_max_twilight_I_Blue_datum 

      outdf.to_csv(outfname, index=False)
      print("============================================")
   return

def dosage_seasons(position):

   # create empty matrices for DJF, MAM, JJA, SON for Twilight, Lunar, ALAN
   DJF_Twilight = np.array([])
   DJF_Lunar = np.array([])
   DJF_ALAN = np.array([])
   
   MAM_Twilight = np.array([])
   MAM_Lunar = np.array([])
   MAM_ALAN = np.array([])
   
   JJA_Twilight = np.array([])
   JJA_Lunar = np.array([])
   JJA_ALAN = np.array([])

   SON_Twilight = np.array([])
   SON_Lunar = np.array([])
   SON_ALAN = np.array([])

   loc_str = []
   short_loc_str = []

   print("Reading in files:")
   for dosage_fname in braceexpand(outdir+'/Dosage_{Tokyo,Mumbai,Shanghai,Lagos,NewYork,LosAngeles,BuenosAires,Plymouth_Dockyard}_2020.csv'):
      print(dosage_fname)
      df_dosage = pd.read_csv(dosage_fname)
      smonth = df_dosage['Month'].to_numpy()
      syear = df_dosage['Year'].to_numpy()
      monthly_solar_BB = df_dosage['solar_I_BB(J/m2)'].to_numpy()
      monthly_lunar_BB = df_dosage['monthly_lunar_I_BB(J/m2)'].to_numpy()
      monthly_ALAN_BB = df_dosage['monthly_ALAN_I_BB(J/m2)'].to_numpy()
      monthly_twilight_BB = df_dosage['monthly_twilight_I_BB(J/m2)'].to_numpy()

      monthly_solar_I_Red = df_dosage['monthly_solar_I_Red(J/m2)'].to_numpy()
      monthly_solar_I_Green = df_dosage['monthly_solar_I_Green(J/m2)'].to_numpy()
      monthly_solar_I_Blue = df_dosage['monthly_solar_I_Blue(J/m2)'].to_numpy()
      monthly_lunar_I_Red = df_dosage['monthly_lunar_I_Red(J/m2)'].to_numpy()
      monthly_lunar_I_Green = df_dosage['monthly_lunar_I_Green(J/m2)'].to_numpy()
      monthly_lunar_I_Blue = df_dosage['monthly_lunar_I_Blue(J/m2)'].to_numpy()
      monthly_ALAN_I_Red = df_dosage['monthly_ALAN_I_Red(J/m2)'].to_numpy()
      monthly_ALAN_I_Green = df_dosage['monthly_ALAN_I_Green(J/m2)'].to_numpy()
      monthly_ALAN_I_Blue = df_dosage['monthly_ALAN_I_Blue(J/m2)'].to_numpy()
      monthly_twilight_I_Red = df_dosage['monthly_twilight_I_Red(J/m2)'].to_numpy()
      monthly_twilight_I_Green = df_dosage['monthly_twilight_I_Green(J/m2)'].to_numpy()
      monthly_twilight_I_Blue = df_dosage['monthly_twilight_I_Blue(J/m2)'].to_numpy()

      monthly_solar_I_BB_datum = df_dosage['monthly_solar_I_BB_datum(J/m2)'].to_numpy()
      monthly_lunar_I_BB_datum = df_dosage['monthly_lunar_I_BB_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_BB_datum = df_dosage['monthly_ALAN_I_BB_datum(J/m2)'].to_numpy()
      monthly_twilight_I_BB_datum = df_dosage['monthly_twilight_I_BB_datum(J/m2)'].to_numpy()
      monthly_solar_I_Red_datum = df_dosage['monthly_solar_I_Red_datum(J/m2)'].to_numpy()
      monthly_solar_I_Green_datum = df_dosage['monthly_solar_I_Green_datum(J/m2)'].to_numpy()
      monthly_solar_I_Blue_datum = df_dosage['monthly_solar_I_Blue_datum(J/m2)'].to_numpy()
      monthly_lunar_I_Red_datum = df_dosage['monthly_lunar_I_Red_datum(J/m2)'].to_numpy()
      monthly_lunar_I_Green_datum = df_dosage['monthly_lunar_I_Green_datum(J/m2)'].to_numpy()
      monthly_lunar_I_Blue_datum = df_dosage['monthly_lunar_I_Blue_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_Red_datum = df_dosage['monthly_ALAN_I_Red_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_Green_datum = df_dosage['monthly_ALAN_I_Green_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_Blue_datum = df_dosage['monthly_ALAN_I_Blue_datum(J/m2)'].to_numpy()
      monthly_twilight_I_Red_datum = df_dosage['monthly_twilight_I_Red_datum(J/m2)'].to_numpy()
      monthly_twilight_I_Green_datum = df_dosage['monthly_twilight_I_Green_datum(J/m2)'].to_numpy()
      monthly_twilight_I_Blue_datum = df_dosage['monthly_twilight_I_Blue_datum(J/m2)'].to_numpy()
       
      smonth = smonth[0:12]
      syear = str(syear[0])
      monthly_solar_BB = monthly_solar_BB[0:12]
      monthly_lunar_BB = monthly_lunar_BB[0:12]
      monthly_ALAN_BB = monthly_ALAN_BB[0:12]
      monthly_twilight_BB = monthly_twilight_BB[0:12]

      monthly_solar_I_Red = monthly_solar_I_Red[0:12]
      monthly_solar_I_Green = monthly_solar_I_Green[0:12]
      monthly_solar_I_Blue = monthly_solar_I_Blue[0:12]
      monthly_lunar_I_Red = monthly_lunar_I_Red[0:12]
      monthly_lunar_I_Green = monthly_lunar_I_Green[0:12]
      monthly_lunar_I_Blue = monthly_lunar_I_Blue[0:12]
      monthly_ALAN_I_Red = monthly_ALAN_I_Red[0:12]
      monthly_ALAN_I_Green = monthly_ALAN_I_Green[0:12]
      monthly_ALAN_I_Blue = monthly_ALAN_I_Blue[0:12]
      monthly_twilight_I_Red = monthly_twilight_I_Red[0:12]
      monthly_twilight_I_Green = monthly_twilight_I_Green[0:12]
      monthly_twilight_I_Blue = monthly_twilight_I_Blue[0:12]

      monthly_solar_I_BB_datum = monthly_solar_I_BB_datum[0:12]
      monthly_lunar_I_BB_datum = monthly_lunar_I_BB_datum[0:12]
      monthly_ALAN_I_BB_datum = monthly_ALAN_I_BB_datum[0:12]
      monthly_twilight_I_BB_datum = monthly_twilight_I_BB_datum[0:12]

      monthly_solar_I_Red_datum = monthly_solar_I_Red_datum[0:12]
      monthly_solar_I_Green_datum = monthly_solar_I_Green_datum[0:12]
      monthly_solar_I_Blue_datum = monthly_solar_I_Blue_datum[0:12]
      monthly_lunar_I_Red_datum = monthly_lunar_I_Red_datum[0:12]
      monthly_lunar_I_Green_datum = monthly_lunar_I_Green_datum[0:12]
      monthly_lunar_I_Blue_datum = monthly_lunar_I_Blue_datum[0:12]
      monthly_ALAN_I_Red_datum = monthly_ALAN_I_Red_datum[0:12]
      monthly_ALAN_I_Green_datum = monthly_ALAN_I_Green_datum[0:12]
      monthly_ALAN_I_Blue_datum = monthly_ALAN_I_Blue_datum[0:12]
      monthly_twilight_I_Red_datum = monthly_twilight_I_Red_datum[0:12]
      monthly_twilight_I_Green_datum = monthly_twilight_I_Green_datum[0:12]
      monthly_twilight_I_Blue_datum = monthly_twilight_I_Blue_datum[0:12]
   
      DJF_index = np.where((smonth == 'Dec') | (smonth == 'Jan') | (smonth == 'Feb'))
      MAM_index = np.where((smonth == 'Mar') | (smonth == 'Apr') | (smonth == 'May'))
      JJA_index = np.where((smonth == 'Jun') | (smonth == 'Jul') | (smonth == 'Aug'))
      SON_index = np.where((smonth == 'Sep') | (smonth == 'Oct') | (smonth == 'Nov'))

      if position == "intertidal":
         DJF_Twilight = np.append(DJF_Twilight,np.mean(monthly_twilight_I_BB_datum[DJF_index]))
         DJF_Lunar = np.append(DJF_Lunar,np.mean(monthly_lunar_I_BB_datum[DJF_index]))
         DJF_ALAN = np.append(DJF_ALAN,np.mean(monthly_ALAN_I_BB_datum[DJF_index]))
   
         MAM_Twilight = np.append(MAM_Twilight,np.mean(monthly_twilight_I_BB_datum[MAM_index]))
         MAM_Lunar = np.append(MAM_Lunar,np.mean(monthly_lunar_I_BB_datum[MAM_index]))
         MAM_ALAN = np.append(MAM_ALAN,np.mean(monthly_ALAN_I_BB_datum[MAM_index]))
   
         JJA_Twilight = np.append(JJA_Twilight,np.mean(monthly_twilight_I_BB_datum[JJA_index]))
         JJA_Lunar = np.append(JJA_Lunar,np.mean(monthly_lunar_I_BB_datum[JJA_index]))
         JJA_ALAN = np.append(JJA_ALAN,np.mean(monthly_ALAN_I_BB_datum[JJA_index]))

         SON_Twilight = np.append(SON_Twilight,np.mean(monthly_twilight_I_BB_datum[SON_index]))
         SON_Lunar = np.append(SON_Lunar,np.mean(monthly_lunar_I_BB_datum[SON_index]))
         SON_ALAN = np.append(SON_ALAN,np.mean(monthly_ALAN_I_BB_datum[SON_index]))

      if position == "surface":
         DJF_Twilight = np.append(DJF_Twilight,np.mean(monthly_twilight_BB[DJF_index]))
         DJF_Lunar = np.append(DJF_Lunar,np.mean(monthly_lunar_BB[DJF_index]))
         DJF_ALAN = np.append(DJF_ALAN,np.mean(monthly_ALAN_BB[DJF_index]))
   
         MAM_Twilight = np.append(MAM_Twilight,np.mean(monthly_twilight_BB[MAM_index]))
         MAM_Lunar = np.append(MAM_Lunar,np.mean(monthly_lunar_BB[MAM_index]))
         MAM_ALAN = np.append(MAM_ALAN,np.mean(monthly_ALAN_BB[MAM_index]))
   
         JJA_Twilight = np.append(JJA_Twilight,np.mean(monthly_twilight_BB[JJA_index]))
         JJA_Lunar = np.append(JJA_Lunar,np.mean(monthly_lunar_BB[JJA_index]))
         JJA_ALAN = np.append(JJA_ALAN,np.mean(monthly_ALAN_BB[JJA_index]))

         SON_Twilight = np.append(SON_Twilight,np.mean(monthly_twilight_BB[SON_index]))
         SON_Lunar = np.append(SON_Lunar,np.mean(monthly_lunar_BB[SON_index]))
         SON_ALAN = np.append(SON_ALAN,np.mean(monthly_ALAN_BB[SON_index]))
         
      loc_str.append(df_dosage['Location'][0])
      short_loc_str.append(df_dosage['Location'][0][0:2])

   barWidth = 0.25
   
   fig, axs = plt.subplots(2,2, figsize=(10,7))
   fig.suptitle("Seasonal - "+position, fontsize=16)

   # DJF - sort in (descending) order of ALAN dosage 
   short_loc_array = np.array(short_loc_str)
   alan_sort = np.argsort(DJF_ALAN)[::-1]

   bars1 = DJF_ALAN[alan_sort]
   bars2 = DJF_Lunar[alan_sort]
   bars3 = DJF_Twilight[alan_sort]/1e+3
   short_loc_array = short_loc_array[alan_sort]

   r1 = np.arange(len(DJF_Lunar))
   r2 = [x + barWidth for x in r1]
   r3 = [x + barWidth for x in r2]

   pcm = axs[0,0].bar(r1, bars1, color='#ff7f0e', width=barWidth, edgecolor='k', label='ALAN')
   pcm = axs[0,0].bar(r2, bars2, color='#C5C9C7', width=barWidth, edgecolor='k', label='Lunar')
   pcm = axs[0,0].bar(r3, bars3, color='#00008B', width=barWidth, edgecolor='k', label='Twilight')
  
   axs[0,0].xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(short_loc_array))))
   axs[0,0].set(ylabel='Dosage [J] or [KJ]'); axs[0,0].set(xticklabels=short_loc_array)
   axs[0,0].set(title='A) DJF')
   axs[0,0].set_ylim(top = 550, bottom = 0);   
   axs[0,0].legend()
   
   # MAM - sort in (descending) order of ALAN dosage 
   short_loc_array = np.array(short_loc_str)
   alan_sort = np.argsort(MAM_ALAN)[::-1]

   bars1 = MAM_ALAN[alan_sort]
   bars2 = MAM_Lunar[alan_sort]
   bars3 = MAM_Twilight[alan_sort]/1e+3
   short_loc_array = short_loc_array[alan_sort]

   r1 = np.arange(len(MAM_Lunar))
   r2 = [x + barWidth for x in r1]
   r3 = [x + barWidth for x in r2]

   pcm = axs[0,1].bar(r1, bars1, color='#ff7f0e', width=barWidth, edgecolor='k', label='ALAN')
   pcm = axs[0,1].bar(r2, bars2, color='#C5C9C7', width=barWidth, edgecolor='k', label='Lunar')
   pcm = axs[0,1].bar(r3, bars3, color='#00008B', width=barWidth, edgecolor='k', label='Twilight')
  
   axs[0,1].xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(short_loc_array))))
   axs[0,1].set(ylabel='Dosage [J] or [KJ]'); axs[0,1].set(xticklabels=short_loc_array)
   axs[0,1].set(title='B) MAM')
   axs[0,1].set_ylim(top = 550, bottom = 0);   
   #axs[0,1].legend()
   
   # JJA - sort in (descending) order of ALAN dosage 
   short_loc_array = np.array(short_loc_str)
   alan_sort = np.argsort(JJA_ALAN)[::-1]

   bars1 = JJA_ALAN[alan_sort]
   bars2 = JJA_Lunar[alan_sort]
   bars3 = JJA_Twilight[alan_sort]/1e+3
   short_loc_array = short_loc_array[alan_sort]

   r1 = np.arange(len(JJA_Lunar))
   r2 = [x + barWidth for x in r1]
   r3 = [x + barWidth for x in r2]

   pcm = axs[1,0].bar(r1, bars1, color='#ff7f0e', width=barWidth, edgecolor='k', label='ALAN')
   pcm = axs[1,0].bar(r2, bars2, color='#C5C9C7', width=barWidth, edgecolor='k', label='Lunar')
   pcm = axs[1,0].bar(r3, bars3, color='#00008B', width=barWidth, edgecolor='k', label='Twilight')
  
   axs[1,0].xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(short_loc_array))))
   axs[1,0].set(ylabel='Dosage [J] or [KJ]'); axs[1,0].set(xlabel='City', xticklabels=short_loc_array)
   axs[1,0].set(title='C) JJA')
   axs[1,0].set_ylim(top = 550, bottom = 0);   
   #axs[1,0].legend()

   # SON - sort in (descending) order of ALAN dosage 
   short_loc_array = np.array(short_loc_str)
   alan_sort = np.argsort(SON_ALAN)[::-1]

   bars1 = SON_ALAN[alan_sort]
   bars2 = SON_Lunar[alan_sort]
   bars3 = SON_Twilight[alan_sort]/1e+3
   short_loc_array = short_loc_array[alan_sort]

   r1 = np.arange(len(SON_Lunar))
   r2 = [x + barWidth for x in r1]
   r3 = [x + barWidth for x in r2]

   pcm = axs[1,1].bar(r1, bars1, color='#ff7f0e', width=barWidth, edgecolor='k', label='ALAN')
   pcm = axs[1,1].bar(r2, bars2, color='#C5C9C7', width=barWidth, edgecolor='k', label='Lunar')
   pcm = axs[1,1].bar(r3, bars3, color='#00008B', width=barWidth, edgecolor='k', label='Twilight')
  
   axs[1,1].xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(short_loc_array))))
   axs[1,1].set(ylabel='Dosage [J] or [KJ]'); axs[1,1].set(xlabel='City', xticklabels=short_loc_array)
   axs[1,1].set(title='D) SON')
   axs[1,1].set_ylim(top = 550, bottom = 0);   
   #axs[1,1].legend()
   
   pdb.set_trace()

   fig.savefig(outdir+'/Dosage_seasonal_'+position+'_'+syear+'.png')
   plt.close(fig)
   print('Produced image: '+outdir+'/Dosage_seasonal_'+position+'_'+syear+'.png')
   
   return   

def max_seasons(position):

   # create empty matrices for DJF, MAM, JJA, SON for Twilight, Lunar, ALAN
   DJF_Twilight = np.array([])
   DJF_Lunar = np.array([])
   DJF_ALAN = np.array([])
   
   MAM_Twilight = np.array([])
   MAM_Lunar = np.array([])
   MAM_ALAN = np.array([])
   
   JJA_Twilight = np.array([])
   JJA_Lunar = np.array([])
   JJA_ALAN = np.array([])

   SON_Twilight = np.array([])
   SON_Lunar = np.array([])
   SON_ALAN = np.array([])

   loc_str = []
   short_loc_str = []

   print("Reading in files:")
   for dosage_fname in braceexpand(outdir+'/Dosage_{Tokyo,Mumbai,Shanghai,Lagos,NewYork,LosAngeles,BuenosAires,Plymouth_Dockyard}_2020.csv'):
      print(dosage_fname)
      df_dosage = pd.read_csv(dosage_fname)
      smonth = df_dosage['Month'].to_numpy()
      syear = df_dosage['Year'].to_numpy()
      monthly_max_solar_BB = df_dosage['solar_I_BB(W/m2)'].to_numpy()
      monthly_max_lunar_BB = df_dosage['monthly_max_lunar_I_BB(W/m2)'].to_numpy()
      monthly_max_ALAN_BB = df_dosage['monthly_max_ALAN_I_BB(W/m2)'].to_numpy()
      monthly_max_twilight_BB = df_dosage['monthly_max_twilight_I_BB(W/m2)'].to_numpy()

      monthly_max_solar_I_Red = df_dosage['monthly_max_solar_I_Red(W/m2)'].to_numpy()
      monthly_max_solar_I_Green = df_dosage['monthly_max_solar_I_Green(W/m2)'].to_numpy()
      monthly_max_solar_I_Blue = df_dosage['monthly_max_solar_I_Blue(W/m2)'].to_numpy()
      monthly_max_lunar_I_Red = df_dosage['monthly_max_lunar_I_Red(W/m2)'].to_numpy()
      monthly_max_lunar_I_Green = df_dosage['monthly_max_lunar_I_Green(W/m2)'].to_numpy()
      monthly_max_lunar_I_Blue = df_dosage['monthly_max_lunar_I_Blue(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Red = df_dosage['monthly_max_ALAN_I_Red(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Green = df_dosage['monthly_max_ALAN_I_Green(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Blue = df_dosage['monthly_max_ALAN_I_Blue(W/m2)'].to_numpy()
      monthly_max_twilight_I_Red = df_dosage['monthly_max_twilight_I_Red(W/m2)'].to_numpy()
      monthly_max_twilight_I_Green = df_dosage['monthly_max_twilight_I_Green(W/m2)'].to_numpy()
      monthly_max_twilight_I_Blue = df_dosage['monthly_max_twilight_I_Blue(W/m2)'].to_numpy()

      monthly_max_solar_I_BB_datum = df_dosage['monthly_max_solar_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_BB_datum = df_dosage['monthly_max_lunar_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_BB_datum = df_dosage['monthly_max_ALAN_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_BB_datum = df_dosage['monthly_max_twilight_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_solar_I_Red_datum = df_dosage['monthly_max_solar_I_Red_datum(W/m2)'].to_numpy()
      monthly_max_solar_I_Green_datum = df_dosage['monthly_max_solar_I_Green_datum(W/m2)'].to_numpy()
      monthly_max_solar_I_Blue_datum = df_dosage['monthly_max_solar_I_Blue_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_Red_datum = df_dosage['monthly_max_lunar_I_Red_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_Green_datum = df_dosage['monthly_max_lunar_I_Green_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_Blue_datum = df_dosage['monthly_max_lunar_I_Blue_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Red_datum = df_dosage['monthly_max_ALAN_I_Red_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Green_datum = df_dosage['monthly_max_ALAN_I_Green_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Blue_datum = df_dosage['monthly_max_ALAN_I_Blue_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_Red_datum = df_dosage['monthly_max_twilight_I_Red_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_Green_datum = df_dosage['monthly_max_twilight_I_Green_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_Blue_datum = df_dosage['monthly_max_twilight_I_Blue_datum(W/m2)'].to_numpy()
       
      smonth = smonth[0:12]
      syear = str(syear[0])
      monthly_max_solar_BB = monthly_max_solar_BB[0:12]
      monthly_max_lunar_BB = monthly_max_lunar_BB[0:12]
      monthly_max_ALAN_BB = monthly_max_ALAN_BB[0:12]
      monthly_max_twilight_BB = monthly_max_twilight_BB[0:12]

      monthly_max_solar_I_Red = monthly_max_solar_I_Red[0:12]
      monthly_max_solar_I_Green = monthly_max_solar_I_Green[0:12]
      monthly_max_solar_I_Blue = monthly_max_solar_I_Blue[0:12]
      monthly_max_lunar_I_Red = monthly_max_lunar_I_Red[0:12]
      monthly_max_lunar_I_Green = monthly_max_lunar_I_Green[0:12]
      monthly_max_lunar_I_Blue = monthly_max_lunar_I_Blue[0:12]
      monthly_max_ALAN_I_Red = monthly_max_ALAN_I_Red[0:12]
      monthly_max_ALAN_I_Green = monthly_max_ALAN_I_Green[0:12]
      monthly_max_ALAN_I_Blue = monthly_max_ALAN_I_Blue[0:12]
      monthly_max_twilight_I_Red = monthly_max_twilight_I_Red[0:12]
      monthly_max_twilight_I_Green = monthly_max_twilight_I_Green[0:12]
      monthly_max_twilight_I_Blue = monthly_max_twilight_I_Blue[0:12]

      monthly_max_solar_I_BB_datum = monthly_max_solar_I_BB_datum[0:12]
      monthly_max_lunar_I_BB_datum = monthly_max_lunar_I_BB_datum[0:12]
      monthly_max_ALAN_I_BB_datum = monthly_max_ALAN_I_BB_datum[0:12]
      monthly_max_twilight_I_BB_datum = monthly_max_twilight_I_BB_datum[0:12]

      monthly_max_solar_I_Red_datum = monthly_max_solar_I_Red_datum[0:12]
      monthly_max_solar_I_Green_datum = monthly_max_solar_I_Green_datum[0:12]
      monthly_max_solar_I_Blue_datum = monthly_max_solar_I_Blue_datum[0:12]
      monthly_max_lunar_I_Red_datum = monthly_max_lunar_I_Red_datum[0:12]
      monthly_max_lunar_I_Green_datum = monthly_max_lunar_I_Green_datum[0:12]
      monthly_max_lunar_I_Blue_datum = monthly_max_lunar_I_Blue_datum[0:12]
      monthly_max_ALAN_I_Red_datum = monthly_max_ALAN_I_Red_datum[0:12]
      monthly_max_ALAN_I_Green_datum = monthly_max_ALAN_I_Green_datum[0:12]
      monthly_max_ALAN_I_Blue_datum = monthly_max_ALAN_I_Blue_datum[0:12]
      monthly_max_twilight_I_Red_datum = monthly_max_twilight_I_Red_datum[0:12]
      monthly_max_twilight_I_Green_datum = monthly_max_twilight_I_Green_datum[0:12]
      monthly_max_twilight_I_Blue_datum = monthly_max_twilight_I_Blue_datum[0:12]
   
      DJF_index = np.where((smonth == 'Dec') | (smonth == 'Jan') | (smonth == 'Feb'))
      MAM_index = np.where((smonth == 'Mar') | (smonth == 'Apr') | (smonth == 'May'))
      JJA_index = np.where((smonth == 'Jun') | (smonth == 'Jul') | (smonth == 'Aug'))
      SON_index = np.where((smonth == 'Sep') | (smonth == 'Oct') | (smonth == 'Nov'))

      if position == "intertidal":
         DJF_Twilight = np.append(DJF_Twilight,np.mean(monthly_max_twilight_I_BB_datum[DJF_index]))
         DJF_Lunar = np.append(DJF_Lunar,np.mean(monthly_max_lunar_I_BB_datum[DJF_index])/1e-6)
         DJF_ALAN = np.append(DJF_ALAN,np.mean(monthly_max_ALAN_I_BB_datum[DJF_index])/1e-6)
   
         MAM_Twilight = np.append(MAM_Twilight,np.mean(monthly_max_twilight_I_BB_datum[MAM_index]))
         MAM_Lunar = np.append(MAM_Lunar,np.mean(monthly_max_lunar_I_BB_datum[MAM_index])/1e-6)
         MAM_ALAN = np.append(MAM_ALAN,np.mean(monthly_max_ALAN_I_BB_datum[MAM_index])/1e-6)
   
         JJA_Twilight = np.append(JJA_Twilight,np.mean(monthly_max_twilight_I_BB_datum[JJA_index]))
         JJA_Lunar = np.append(JJA_Lunar,np.mean(monthly_max_lunar_I_BB_datum[JJA_index])/1e-6)
         JJA_ALAN = np.append(JJA_ALAN,np.mean(monthly_max_ALAN_I_BB_datum[JJA_index])/1e-6)

         SON_Twilight = np.append(SON_Twilight,np.mean(monthly_max_twilight_I_BB_datum[SON_index]))
         SON_Lunar = np.append(SON_Lunar,np.mean(monthly_max_lunar_I_BB_datum[SON_index])/1e-6)
         SON_ALAN = np.append(SON_ALAN,np.mean(monthly_max_ALAN_I_BB_datum[SON_index])/1e-6)

      if position == "surface":
         DJF_Twilight = np.append(DJF_Twilight,np.mean(monthly_max_twilight_BB[DJF_index]))
         DJF_Lunar = np.append(DJF_Lunar,np.mean(monthly_max_lunar_BB[DJF_index])/1e-6)
         DJF_ALAN = np.append(DJF_ALAN,np.mean(monthly_max_ALAN_BB[DJF_index])/1e-6)
   
         MAM_Twilight = np.append(MAM_Twilight,np.mean(monthly_max_twilight_BB[MAM_index]))
         MAM_Lunar = np.append(MAM_Lunar,np.mean(monthly_max_lunar_BB[MAM_index])/1e-6)
         MAM_ALAN = np.append(MAM_ALAN,np.mean(monthly_max_ALAN_BB[MAM_index])/1e-6)
   
         JJA_Twilight = np.append(JJA_Twilight,np.mean(monthly_max_twilight_BB[JJA_index]))
         JJA_Lunar = np.append(JJA_Lunar,np.mean(monthly_max_lunar_BB[JJA_index])/1e-6)
         JJA_ALAN = np.append(JJA_ALAN,np.mean(monthly_max_ALAN_BB[JJA_index])/1e-6)

         SON_Twilight = np.append(SON_Twilight,np.mean(monthly_max_twilight_BB[SON_index]))
         SON_Lunar = np.append(SON_Lunar,np.mean(monthly_max_lunar_BB[SON_index])/1e-6)
         SON_ALAN = np.append(SON_ALAN,np.mean(monthly_max_ALAN_BB[SON_index])/1e-6)
         
      loc_str.append(df_dosage['Location'][0])
      short_loc_str.append(df_dosage['Location'][0][0:2])

   barWidth = 0.25
   
   #pdb.set_trace()
   
   fig, axs = plt.subplots(2,2, figsize=(10,7))
   fig.suptitle("Seasonal - "+position, fontsize=16)

   # DJF - sort in (descending) order of ALAN dosage 
   short_loc_array = np.array(short_loc_str)
   alan_sort = np.argsort(DJF_ALAN)[::-1]

   bars1 = DJF_ALAN[alan_sort]
   bars2 = DJF_Lunar[alan_sort]
   bars3 = DJF_Twilight[alan_sort]/1e+3
   short_loc_array = short_loc_array[alan_sort]

   r1 = np.arange(len(DJF_Lunar))
   r2 = [x + barWidth for x in r1]
   r3 = [x + barWidth for x in r2]

   pcm = axs[0,0].bar(r1, bars1, color='#ff7f0e', width=barWidth, edgecolor='k', label='ALAN')
   pcm = axs[0,0].bar(r2, bars2, color='#C5C9C7', width=barWidth, edgecolor='k', label='Lunar')
   #pcm = axs[0,0].bar(r3, bars3, color='#00008B', width=barWidth, edgecolor='k', label='Twilight')
  
   axs[0,0].xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(short_loc_array))))
   axs[0,0].set(ylabel='Mean max irradiance [$\mu$W/m$^2$]'); axs[0,0].set(xticklabels=short_loc_array)
   axs[0,0].set(title='A) DJF')
   axs[0,0].set_ylim(top = 1500, bottom = 0);   
   axs[0,0].legend()
   
   # MAM - sort in (descending) order of ALAN dosage 
   short_loc_array = np.array(short_loc_str)
   alan_sort = np.argsort(MAM_ALAN)[::-1]

   bars1 = MAM_ALAN[alan_sort]
   bars2 = MAM_Lunar[alan_sort]
   bars3 = MAM_Twilight[alan_sort]/1e+3
   short_loc_array = short_loc_array[alan_sort]

   r1 = np.arange(len(MAM_Lunar))
   r2 = [x + barWidth for x in r1]
   r3 = [x + barWidth for x in r2]

   pcm = axs[0,1].bar(r1, bars1, color='#ff7f0e', width=barWidth, edgecolor='k', label='ALAN')
   pcm = axs[0,1].bar(r2, bars2, color='#C5C9C7', width=barWidth, edgecolor='k', label='Lunar')
   #pcm = axs[0,1].bar(r3, bars3, color='#00008B', width=barWidth, edgecolor='k', label='Twilight')
  
   axs[0,1].xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(short_loc_array))))
   axs[0,1].set(ylabel='Mean max irradiance [$\mu$W/m$^2$]'); axs[0,1].set(xticklabels=short_loc_array)
   axs[0,1].set(title='B) MAM')
   axs[0,1].set_ylim(top = 1500, bottom = 0);   
   #axs[0,1].legend()
   
   # JJA - sort in (descending) order of ALAN dosage 
   short_loc_array = np.array(short_loc_str)
   alan_sort = np.argsort(JJA_ALAN)[::-1]

   bars1 = JJA_ALAN[alan_sort]
   bars2 = JJA_Lunar[alan_sort]
   bars3 = JJA_Twilight[alan_sort]/1e+3
   short_loc_array = short_loc_array[alan_sort]

   r1 = np.arange(len(JJA_Lunar))
   r2 = [x + barWidth for x in r1]
   r3 = [x + barWidth for x in r2]

   pcm = axs[1,0].bar(r1, bars1, color='#ff7f0e', width=barWidth, edgecolor='k', label='ALAN')
   pcm = axs[1,0].bar(r2, bars2, color='#C5C9C7', width=barWidth, edgecolor='k', label='Lunar')
   #pcm = axs[1,0].bar(r3, bars3, color='#00008B', width=barWidth, edgecolor='k', label='Twilight')
  
   axs[1,0].xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(short_loc_array))))
   axs[1,0].set(ylabel='Mean max irradiance [$\mu$W/m$^2$]'); axs[1,0].set(xlabel='City', xticklabels=short_loc_array)
   axs[1,0].set(title='C) JJA')
   axs[1,0].set_ylim(top = 1500, bottom = 0);   
   #axs[1,0].legend()

   # SON - sort in (descending) order of ALAN dosage 
   short_loc_array = np.array(short_loc_str)
   alan_sort = np.argsort(SON_ALAN)[::-1]

   bars1 = SON_ALAN[alan_sort]
   bars2 = SON_Lunar[alan_sort]
   bars3 = SON_Twilight[alan_sort]/1e+3
   short_loc_array = short_loc_array[alan_sort]

   r1 = np.arange(len(SON_Lunar))
   r2 = [x + barWidth for x in r1]
   r3 = [x + barWidth for x in r2]

   pcm = axs[1,1].bar(r1, bars1, color='#ff7f0e', width=barWidth, edgecolor='k', label='ALAN')
   pcm = axs[1,1].bar(r2, bars2, color='#C5C9C7', width=barWidth, edgecolor='k', label='Lunar')
   #pcm = axs[1,1].bar(r3, bars3, color='#00008B', width=barWidth, edgecolor='k', label='Twilight')
  
   axs[1,1].xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(short_loc_array))))
   axs[1,1].set(ylabel='Mean max irradiance [$\mu$W/m$^2$]'); axs[1,1].set(xlabel='City', xticklabels=short_loc_array)
   axs[1,1].set(title='D) SON')
   axs[1,1].set_ylim(top = 1500, bottom = 0);   
   #axs[1,1].legend()

   fig.savefig(outdir+'/max_seasonal_'+position+'_'+syear+'.png')
   plt.close(fig)
   print('Produced image: '+outdir+'/max_seasonal_'+position+'_'+syear+'.png')
   
   return   


   
def dosage_plot(location):

   if location:
      print("Location: ", location)

   for dosage_fname in glob.glob(outdir+'/Dosage_'+location+'_2020.csv'):
      df_dosage = pd.read_csv(dosage_fname)
      smonth = df_dosage['Month'].to_numpy()
      syear = df_dosage['Year'].to_numpy()
      monthly_solar_BB = df_dosage['solar_I_BB(J/m2)'].to_numpy()
      monthly_lunar_BB = df_dosage['monthly_lunar_I_BB(J/m2)'].to_numpy()
      monthly_ALAN_BB = df_dosage['monthly_ALAN_I_BB(J/m2)'].to_numpy()
      monthly_twilight_BB = df_dosage['monthly_twilight_I_BB(J/m2)'].to_numpy()

      monthly_solar_I_Red = df_dosage['monthly_solar_I_Red(J/m2)'].to_numpy()
      monthly_solar_I_Green = df_dosage['monthly_solar_I_Green(J/m2)'].to_numpy()
      monthly_solar_I_Blue = df_dosage['monthly_solar_I_Blue(J/m2)'].to_numpy()
      monthly_lunar_I_Red = df_dosage['monthly_lunar_I_Red(J/m2)'].to_numpy()
      monthly_lunar_I_Green = df_dosage['monthly_lunar_I_Green(J/m2)'].to_numpy()
      monthly_lunar_I_Blue = df_dosage['monthly_lunar_I_Blue(J/m2)'].to_numpy()
      monthly_ALAN_I_Red = df_dosage['monthly_ALAN_I_Red(J/m2)'].to_numpy()
      monthly_ALAN_I_Green = df_dosage['monthly_ALAN_I_Green(J/m2)'].to_numpy()
      monthly_ALAN_I_Blue = df_dosage['monthly_ALAN_I_Blue(J/m2)'].to_numpy()
      monthly_twilight_I_Red = df_dosage['monthly_twilight_I_Red(J/m2)'].to_numpy()
      monthly_twilight_I_Green = df_dosage['monthly_twilight_I_Green(J/m2)'].to_numpy()
      monthly_twilight_I_Blue = df_dosage['monthly_twilight_I_Blue(J/m2)'].to_numpy()

      monthly_solar_I_BB_datum = df_dosage['monthly_solar_I_BB_datum(J/m2)'].to_numpy()
      monthly_lunar_I_BB_datum = df_dosage['monthly_lunar_I_BB_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_BB_datum = df_dosage['monthly_ALAN_I_BB_datum(J/m2)'].to_numpy()
      monthly_twilight_I_BB_datum = df_dosage['monthly_twilight_I_BB_datum(J/m2)'].to_numpy()
      monthly_solar_I_Red_datum = df_dosage['monthly_solar_I_Red_datum(J/m2)'].to_numpy()
      monthly_solar_I_Green_datum = df_dosage['monthly_solar_I_Green_datum(J/m2)'].to_numpy()
      monthly_solar_I_Blue_datum = df_dosage['monthly_solar_I_Blue_datum(J/m2)'].to_numpy()
      monthly_lunar_I_Red_datum = df_dosage['monthly_lunar_I_Red_datum(J/m2)'].to_numpy()
      monthly_lunar_I_Green_datum = df_dosage['monthly_lunar_I_Green_datum(J/m2)'].to_numpy()
      monthly_lunar_I_Blue_datum = df_dosage['monthly_lunar_I_Blue_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_Red_datum = df_dosage['monthly_ALAN_I_Red_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_Green_datum = df_dosage['monthly_ALAN_I_Green_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_Blue_datum = df_dosage['monthly_ALAN_I_Blue_datum(J/m2)'].to_numpy()
      monthly_twilight_I_Red_datum = df_dosage['monthly_twilight_I_Red_datum(J/m2)'].to_numpy()
      monthly_twilight_I_Green_datum = df_dosage['monthly_twilight_I_Green_datum(J/m2)'].to_numpy()
      monthly_twilight_I_Blue_datum = df_dosage['monthly_twilight_I_Blue_datum(J/m2)'].to_numpy()
       
      smonth = smonth[0:12]
      syear = str(syear[0])
      monthly_solar_BB = monthly_solar_BB[0:12]
      monthly_lunar_BB = monthly_lunar_BB[0:12]
      monthly_ALAN_BB = monthly_ALAN_BB[0:12]
      monthly_twilight_BB = monthly_twilight_BB[0:12]

      monthly_solar_I_Red = monthly_solar_I_Red[0:12]
      monthly_solar_I_Green = monthly_solar_I_Green[0:12]
      monthly_solar_I_Blue = monthly_solar_I_Blue[0:12]
      monthly_lunar_I_Red = monthly_lunar_I_Red[0:12]
      monthly_lunar_I_Green = monthly_lunar_I_Green[0:12]
      monthly_lunar_I_Blue = monthly_lunar_I_Blue[0:12]
      monthly_ALAN_I_Red = monthly_ALAN_I_Red[0:12]
      monthly_ALAN_I_Green = monthly_ALAN_I_Green[0:12]
      monthly_ALAN_I_Blue = monthly_ALAN_I_Blue[0:12]
      monthly_twilight_I_Red = monthly_twilight_I_Red[0:12]
      monthly_twilight_I_Green = monthly_twilight_I_Green[0:12]
      monthly_twilight_I_Blue = monthly_twilight_I_Blue[0:12]

      monthly_solar_I_BB_datum = monthly_solar_I_BB_datum[0:12]
      monthly_lunar_I_BB_datum = monthly_lunar_I_BB_datum[0:12]
      monthly_ALAN_I_BB_datum = monthly_ALAN_I_BB_datum[0:12]
      monthly_twilight_I_BB_datum = monthly_twilight_I_BB_datum[0:12]

      monthly_solar_I_Red_datum = monthly_solar_I_Red_datum[0:12]
      monthly_solar_I_Green_datum = monthly_solar_I_Green_datum[0:12]
      monthly_solar_I_Blue_datum = monthly_solar_I_Blue_datum[0:12]
      monthly_lunar_I_Red_datum = monthly_lunar_I_Red_datum[0:12]
      monthly_lunar_I_Green_datum = monthly_lunar_I_Green_datum[0:12]
      monthly_lunar_I_Blue_datum = monthly_lunar_I_Blue_datum[0:12]
      monthly_ALAN_I_Red_datum = monthly_ALAN_I_Red_datum[0:12]
      monthly_ALAN_I_Green_datum = monthly_ALAN_I_Green_datum[0:12]
      monthly_ALAN_I_Blue_datum = monthly_ALAN_I_Blue_datum[0:12]
      monthly_twilight_I_Red_datum = monthly_twilight_I_Red_datum[0:12]
      monthly_twilight_I_Green_datum = monthly_twilight_I_Green_datum[0:12]
      monthly_twilight_I_Blue_datum = monthly_twilight_I_Blue_datum[0:12]

      fig, axs = plt.subplots(2,2, figsize=(10,7))
      fig.suptitle(location+" "+syear, fontsize=16)

      pcm = axs[0,0].plot(smonth, monthly_solar_BB/1e+6, 'k', label='Broadband')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Red/1e+6, 'r', label='Red(620-740nm)')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Green/1e+6, 'g', label='Green(495-560nm)')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Blue/1e+6, 'b', label='Blue(400-500nm)')
      #pcm = axs[0,0].plot(smonth, (monthly_solar_I_Blue+monthly_solar_I_Red+monthly_solar_I_Green)/1e+6, 'gray')
      axs[0,0].set_title('Solar - surface'); axs[0,0].set(ylabel='Dosage [MJ]'); #axs[0,0].set(xlabel='Month')
      axs[0,0].set_ylim(top = 500, bottom = 0);
     
      pcm = axs[0,1].plot(smonth, monthly_twilight_BB/1e+3, 'k', label='Broadband')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Red/1e+3, 'r', label='Red(620-740nm)')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Green/1e+3, 'g', label='Green(495-560nm)')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Blue/1e+3, 'b', label='Blue(400-500nm)')
      #pcm = axs[0,1].plot(smonth, (monthly_twilight_I_Blue+monthly_twilight_I_Red+monthly_twilight_I_Green)/1e+3, 'gray')
      axs[0,1].set_title('Twilight - surface'); axs[0,1].set(ylabel='Dosage [kJ]'); #axs[0,1].set(xlabel='Month')
      axs[0,1].legend()
      axs[0,1].set_ylim(top = 500, bottom = 0);
      
      pcm = axs[1,0].plot(smonth, monthly_lunar_BB, 'k', label='Broadband')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Red, 'r', label='Red(620-740nm)')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Green, 'g', label='Green(495-560nm)')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Blue, 'b', label='Blue(400-500nm)')
      #pcm = axs[1,0].plot(smonth, (monthly_lunar_I_Blue+monthly_lunar_I_Red+monthly_lunar_I_Green), 'gray')
      axs[1,0].set_title('Lunar - surface'); axs[1,0].set(ylabel='Dosage [J]'); axs[1,0].set(xlabel='Month')
      axs[1,0].set_ylim(top = 550, bottom = 0);
   
      pcm = axs[1,1].plot(smonth, monthly_ALAN_BB, 'k', label='Broadband')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Red, 'r', label='Red(620-740nm)')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Green, 'g', label='Green(495-560nm)')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Blue, 'b', label='Blue(400-500nm)')
      #pcm = axs[1,1].plot(smonth, (monthly_ALAN_I_Blue+monthly_ALAN_I_Red+monthly_ALAN_I_Green), 'gray')
      axs[1,1].set_title('ALAN - surface'); axs[1,1].set(ylabel='Dosage [J]'); axs[1,1].set(xlabel='Month')
      axs[1,1].set_ylim(top = 550, bottom = 0);
   
      fig.savefig(outdir+'/Dosage_surface_'+location+'_'+syear+'.png')
      print('Produced image: '+outdir+'/Dosage_surface_'+location+'_'+syear+'.png')
      plt.close(fig)
      
      fig, axs = plt.subplots(2,2, figsize=(10,7))
      fig.suptitle(location+" "+syear, fontsize=16)

      pcm = axs[0,0].plot(smonth, monthly_solar_I_BB_datum/1e+6, 'k', label='Broadband')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Red_datum/1e+6, 'r', label='Red(620-740nm)')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Green_datum/1e+6, 'g', label='Green(495-560nm)')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Blue_datum/1e+6, 'b', label='Blue(400-500nm)')
      axs[0,0].set_title('Solar - intertidal'); axs[0,0].set(ylabel='Dosage [MJ]'); #axs[0,0].set(xlabel='Month')
      axs[0,0].set_ylim(top = 500, bottom = 0);
     
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_BB_datum/1e+3, 'k', label='Broadband')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Red_datum/1e+3, 'r', label='Red(620-740nm)')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Green_datum/1e+3, 'g', label='Green(495-560nm)')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Blue_datum/1e+3, 'b', label='Blue(400-500nm)')
      axs[0,1].set_title('Twilight - intertidal'); axs[0,1].set(ylabel='Dosage [kJ]'); #axs[0,1].set(xlabel='Month')
      axs[0,1].legend()
      axs[0,1].set_ylim(top = 500, bottom = 0);
      
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_BB_datum, 'k', label='Broadband')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Red_datum, 'r', label='Red(620-740nm)')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Green_datum, 'g', label='Green(495-560nm)')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Blue_datum, 'b', label='Blue(400-500nm)')
      axs[1,0].set_title('Lunar - intertidal'); axs[1,0].set(ylabel='Dosage [J]'); axs[1,0].set(xlabel='Month')
      axs[1,0].set_ylim(top = 550, bottom = 0);
   
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_BB_datum, 'k', label='Broadband')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Red_datum, 'r', label='Red(620-740nm)')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Green_datum, 'g', label='Green(495-560nm)')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Blue_datum, 'b', label='Blue(400-500nm)')
      axs[1,1].set_title('ALAN - intertidal'); axs[1,1].set(ylabel='Dosage [J]'); axs[1,1].set(xlabel='Month')
      axs[1,1].set_ylim(top = 550, bottom = 0);
   
      fig.savefig(outdir+'/Dosage_intertidal_'+location+'_'+syear+'.png')
      print('Produced image: '+outdir+'/Dosage_intertidal_'+location+'_'+syear+'.png')
      plt.close(fig)

   return 0

def max_plot(location):

   if location:
      print("Location: ", location)

   for dosage_fname in glob.glob(outdir+'/Dosage_'+location+'_2020.csv'):
      df_dosage = pd.read_csv(dosage_fname)
      smonth = df_dosage['Month'].to_numpy()
      syear = df_dosage['Year'].to_numpy()
      monthly_max_solar_BB = df_dosage['solar_I_BB(W/m2)'].to_numpy()
      monthly_max_lunar_BB = df_dosage['monthly_max_lunar_I_BB(W/m2)'].to_numpy()
      monthly_max_ALAN_BB = df_dosage['monthly_max_ALAN_I_BB(W/m2)'].to_numpy()
      monthly_max_twilight_BB = df_dosage['monthly_max_twilight_I_BB(W/m2)'].to_numpy()

      monthly_max_solar_I_Red = df_dosage['monthly_max_solar_I_Red(W/m2)'].to_numpy()
      monthly_max_solar_I_Green = df_dosage['monthly_max_solar_I_Green(W/m2)'].to_numpy()
      monthly_max_solar_I_Blue = df_dosage['monthly_max_solar_I_Blue(W/m2)'].to_numpy()
      monthly_max_lunar_I_Red = df_dosage['monthly_max_lunar_I_Red(W/m2)'].to_numpy()
      monthly_max_lunar_I_Green = df_dosage['monthly_max_lunar_I_Green(W/m2)'].to_numpy()
      monthly_max_lunar_I_Blue = df_dosage['monthly_max_lunar_I_Blue(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Red = df_dosage['monthly_max_ALAN_I_Red(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Green = df_dosage['monthly_max_ALAN_I_Green(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Blue = df_dosage['monthly_max_ALAN_I_Blue(W/m2)'].to_numpy()
      monthly_max_twilight_I_Red = df_dosage['monthly_max_twilight_I_Red(W/m2)'].to_numpy()
      monthly_max_twilight_I_Green = df_dosage['monthly_max_twilight_I_Green(W/m2)'].to_numpy()
      monthly_max_twilight_I_Blue = df_dosage['monthly_max_twilight_I_Blue(W/m2)'].to_numpy()

      monthly_max_solar_I_BB_datum = df_dosage['monthly_max_solar_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_BB_datum = df_dosage['monthly_max_lunar_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_BB_datum = df_dosage['monthly_max_ALAN_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_BB_datum = df_dosage['monthly_max_twilight_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_solar_I_Red_datum = df_dosage['monthly_max_solar_I_Red_datum(W/m2)'].to_numpy()
      monthly_max_solar_I_Green_datum = df_dosage['monthly_max_solar_I_Green_datum(W/m2)'].to_numpy()
      monthly_max_solar_I_Blue_datum = df_dosage['monthly_max_solar_I_Blue_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_Red_datum = df_dosage['monthly_max_lunar_I_Red_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_Green_datum = df_dosage['monthly_max_lunar_I_Green_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_Blue_datum = df_dosage['monthly_max_lunar_I_Blue_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Red_datum = df_dosage['monthly_max_ALAN_I_Red_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Green_datum = df_dosage['monthly_max_ALAN_I_Green_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_Blue_datum = df_dosage['monthly_max_ALAN_I_Blue_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_Red_datum = df_dosage['monthly_max_twilight_I_Red_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_Green_datum = df_dosage['monthly_max_twilight_I_Green_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_Blue_datum = df_dosage['monthly_max_twilight_I_Blue_datum(W/m2)'].to_numpy()
       
      smonth = smonth[0:12]
      syear = str(syear[0])
      monthly_max_solar_BB = monthly_max_solar_BB[0:12]
      monthly_max_lunar_BB = monthly_max_lunar_BB[0:12]
      monthly_max_ALAN_BB = monthly_max_ALAN_BB[0:12]
      monthly_max_twilight_BB = monthly_max_twilight_BB[0:12]

      monthly_max_solar_I_Red = monthly_max_solar_I_Red[0:12]
      monthly_max_solar_I_Green = monthly_max_solar_I_Green[0:12]
      monthly_max_solar_I_Blue = monthly_max_solar_I_Blue[0:12]
      monthly_max_lunar_I_Red = monthly_max_lunar_I_Red[0:12]
      monthly_max_lunar_I_Green = monthly_max_lunar_I_Green[0:12]
      monthly_max_lunar_I_Blue = monthly_max_lunar_I_Blue[0:12]
      monthly_max_ALAN_I_Red = monthly_max_ALAN_I_Red[0:12]
      monthly_max_ALAN_I_Green = monthly_max_ALAN_I_Green[0:12]
      monthly_max_ALAN_I_Blue = monthly_max_ALAN_I_Blue[0:12]
      monthly_max_twilight_I_Red = monthly_max_twilight_I_Red[0:12]
      monthly_max_twilight_I_Green = monthly_max_twilight_I_Green[0:12]
      monthly_max_twilight_I_Blue = monthly_max_twilight_I_Blue[0:12]

      monthly_max_solar_I_BB_datum = monthly_max_solar_I_BB_datum[0:12]
      monthly_max_lunar_I_BB_datum = monthly_max_lunar_I_BB_datum[0:12]
      monthly_max_ALAN_I_BB_datum = monthly_max_ALAN_I_BB_datum[0:12]
      monthly_max_twilight_I_BB_datum = monthly_max_twilight_I_BB_datum[0:12]

      monthly_max_solar_I_Red_datum = monthly_max_solar_I_Red_datum[0:12]
      monthly_max_solar_I_Green_datum = monthly_max_solar_I_Green_datum[0:12]
      monthly_max_solar_I_Blue_datum = monthly_max_solar_I_Blue_datum[0:12]
      monthly_max_lunar_I_Red_datum = monthly_max_lunar_I_Red_datum[0:12]
      monthly_max_lunar_I_Green_datum = monthly_max_lunar_I_Green_datum[0:12]
      monthly_max_lunar_I_Blue_datum = monthly_max_lunar_I_Blue_datum[0:12]
      monthly_max_ALAN_I_Red_datum = monthly_max_ALAN_I_Red_datum[0:12]
      monthly_max_ALAN_I_Green_datum = monthly_max_ALAN_I_Green_datum[0:12]
      monthly_max_ALAN_I_Blue_datum = monthly_max_ALAN_I_Blue_datum[0:12]
      monthly_max_twilight_I_Red_datum = monthly_max_twilight_I_Red_datum[0:12]
      monthly_max_twilight_I_Green_datum = monthly_max_twilight_I_Green_datum[0:12]
      monthly_max_twilight_I_Blue_datum = monthly_max_twilight_I_Blue_datum[0:12]

      fig, axs = plt.subplots(2,2, figsize=(10,7))
      fig.suptitle(location+" "+syear, fontsize=16)

      pcm = axs[0,0].plot(smonth, monthly_max_solar_BB, 'k', label='Broadband')
      pcm = axs[0,0].plot(smonth, monthly_max_solar_I_Red, 'r', label='Red(620-740nm)')
      pcm = axs[0,0].plot(smonth, monthly_max_solar_I_Green, 'g', label='Green(495-560nm)')
      pcm = axs[0,0].plot(smonth, monthly_max_solar_I_Blue, 'b', label='Blue(400-500nm)')
      #pcm = axs[0,0].plot(smonth, (monthly_max_solar_I_Blue+monthly_max_solar_I_Red+monthly_max_solar_I_Green)/1e+6, 'gray')
      axs[0,0].set_title('Solar - surface'); axs[0,0].set(ylabel='Max irradiance [W/m$^2$]'); #axs[0,0].set(xlabel='Month')
      axs[0,0].set_ylim(top = 500, bottom = 0);
     
      pcm = axs[0,1].plot(smonth, monthly_max_twilight_BB, 'k', label='Broadband')
      pcm = axs[0,1].plot(smonth, monthly_max_twilight_I_Red, 'r', label='Red(620-740nm)')
      pcm = axs[0,1].plot(smonth, monthly_max_twilight_I_Green, 'g', label='Green(495-560nm)')
      pcm = axs[0,1].plot(smonth, monthly_max_twilight_I_Blue, 'b', label='Blue(400-500nm)')
      #pcm = axs[0,1].plot(smonth, (monthly_max_twilight_I_Blue+monthly_max_twilight_I_Red+monthly_max_twilight_I_Green)/1e+3, 'gray')
      axs[0,1].set_title('Twilight - surface'); axs[0,1].set(ylabel='Max irradiance [W/m$^2$]'); #axs[0,1].set(xlabel='Month')
      axs[0,1].legend()
      axs[0,1].set_ylim(top = 10, bottom = 0);
      
      pcm = axs[1,0].plot(smonth, monthly_max_lunar_BB/1e-6, 'k', label='Broadband')
      pcm = axs[1,0].plot(smonth, monthly_max_lunar_I_Red/1e-6, 'r', label='Red(620-740nm)')
      pcm = axs[1,0].plot(smonth, monthly_max_lunar_I_Green/1e-6, 'g', label='Green(495-560nm)')
      pcm = axs[1,0].plot(smonth, monthly_max_lunar_I_Blue/1e-6, 'b', label='Blue(400-500nm)')
      #pcm = axs[1,0].plot(smonth, (monthly_max_lunar_I_Blue+monthly_max_lunar_I_Red+monthly_max_lunar_I_Green), 'gray')
      axs[1,0].set_title('Lunar - surface'); axs[1,0].set(ylabel='Max irradiance [$\mu$W/m$^2$]'); axs[1,0].set(xlabel='Month')
      axs[1,0].set_ylim(top = 1500, bottom = 0);
   
      pcm = axs[1,1].plot(smonth, monthly_max_ALAN_BB/1e-6, 'k', label='Broadband')
      pcm = axs[1,1].plot(smonth, monthly_max_ALAN_I_Red/1e-6, 'r', label='Red(620-740nm)')
      pcm = axs[1,1].plot(smonth, monthly_max_ALAN_I_Green/1e-6, 'g', label='Green(495-560nm)')
      pcm = axs[1,1].plot(smonth, monthly_max_ALAN_I_Blue/1e-6, 'b', label='Blue(400-500nm)')
      #pcm = axs[1,1].plot(smonth, (monthly_max_ALAN_I_Blue+monthly_max_ALAN_I_Red+monthly_max_ALAN_I_Green), 'gray')
      axs[1,1].set_title('ALAN - surface'); axs[1,1].set(ylabel='Max irradiance [$\mu$W/m$^2$]'); axs[1,1].set(xlabel='Month')
      axs[1,1].set_ylim(top = 1500, bottom = 0);
   
      fig.savefig(outdir+'/max_irradiance_surface_'+location+'_'+syear+'.png')
      print('Produced image: '+outdir+'/max_irradiance_surface_'+location+'_'+syear+'.png')
      plt.close(fig)
      
      fig, axs = plt.subplots(2,2, figsize=(10,7))
      fig.suptitle(location+" "+syear, fontsize=16)

      pcm = axs[0,0].plot(smonth, monthly_max_solar_I_BB_datum, 'k', label='Broadband')
      pcm = axs[0,0].plot(smonth, monthly_max_solar_I_Red_datum, 'r', label='Red(620-740nm)')
      pcm = axs[0,0].plot(smonth, monthly_max_solar_I_Green_datum, 'g', label='Green(495-560nm)')
      pcm = axs[0,0].plot(smonth, monthly_max_solar_I_Blue_datum, 'b', label='Blue(400-500nm)')
      axs[0,0].set_title('Solar - intertidal'); axs[0,0].set(ylabel='Max irradiance [W/m$^2$]'); #axs[0,0].set(xlabel='Month')
      axs[0,0].set_ylim(top = 500, bottom = 0);
     
      pcm = axs[0,1].plot(smonth, monthly_max_twilight_I_BB_datum, 'k', label='Broadband')
      pcm = axs[0,1].plot(smonth, monthly_max_twilight_I_Red_datum, 'r', label='Red(620-740nm)')
      pcm = axs[0,1].plot(smonth, monthly_max_twilight_I_Green_datum, 'g', label='Green(495-560nm)')
      pcm = axs[0,1].plot(smonth, monthly_max_twilight_I_Blue_datum, 'b', label='Blue(400-500nm)')
      axs[0,1].set_title('Twilight - intertidal'); axs[0,1].set(ylabel='Max irradiance [W/m$^2$]'); #axs[0,1].set(xlabel='Month')
      axs[0,1].legend()
      axs[0,1].set_ylim(top = 10, bottom = 0);
      
      pcm = axs[1,0].plot(smonth, monthly_max_lunar_I_BB_datum/1e-6, 'k', label='Broadband')
      pcm = axs[1,0].plot(smonth, monthly_max_lunar_I_Red_datum/1e-6, 'r', label='Red(620-740nm)')
      pcm = axs[1,0].plot(smonth, monthly_max_lunar_I_Green_datum/1e-6, 'g', label='Green(495-560nm)')
      pcm = axs[1,0].plot(smonth, monthly_max_lunar_I_Blue_datum/1e-6, 'b', label='Blue(400-500nm)')
      axs[1,0].set_title('Lunar - intertidal'); axs[1,0].set(ylabel='Max irradiance [$\mu$W/m$^2$]'); axs[1,0].set(xlabel='Month')
      axs[1,0].set_ylim(top = 1500, bottom = 0);
   
      pcm = axs[1,1].plot(smonth, monthly_max_ALAN_I_BB_datum/1e-6, 'k', label='Broadband')
      pcm = axs[1,1].plot(smonth, monthly_max_ALAN_I_Red_datum/1e-6, 'r', label='Red(620-740nm)')
      pcm = axs[1,1].plot(smonth, monthly_max_ALAN_I_Green_datum/1e-6, 'g', label='Green(495-560nm)')
      pcm = axs[1,1].plot(smonth, monthly_max_ALAN_I_Blue_datum/1e-6, 'b', label='Blue(400-500nm)')
      axs[1,1].set_title('ALAN - intertidal'); axs[1,1].set(ylabel='Max irradiance [$\mu$W/m$^2$]'); axs[1,1].set(xlabel='Month')
      axs[1,1].set_ylim(top = 1500, bottom = 0);
   
      fig.savefig(outdir+'/max_irradiance_intertidal_'+location+'_'+syear+'.png')
      print('Produced image: '+outdir+'/max_irradiance_intertidal_'+location+'_'+syear+'.png')
      plt.close(fig)

   return 0

def metonic_max_plot():

   for dosage_fname in glob.glob(outdir+'/Dosage_Plymouth_Dockyard_2001-2020.csv'):
      df_dosage = pd.read_csv(dosage_fname)
      smonth = df_dosage['Month'].to_numpy()
      syear = df_dosage['Year'].to_numpy()
      monthly_max_solar_BB = df_dosage['solar_I_BB(W/m2)'].to_numpy()
      monthly_max_lunar_BB = df_dosage['monthly_max_lunar_I_BB(W/m2)'].to_numpy()
      monthly_max_ALAN_BB = df_dosage['monthly_max_ALAN_I_BB(W/m2)'].to_numpy()
      monthly_max_twilight_BB = df_dosage['monthly_max_twilight_I_BB(W/m2)'].to_numpy()

      monthly_max_solar_I_BB_datum = df_dosage['monthly_max_solar_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_lunar_I_BB_datum = df_dosage['monthly_max_lunar_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_ALAN_I_BB_datum = df_dosage['monthly_max_ALAN_I_BB_datum(W/m2)'].to_numpy()
      monthly_max_twilight_I_BB_datum = df_dosage['monthly_max_twilight_I_BB_datum(W/m2)'].to_numpy()

      monthly_solar_BB = df_dosage['solar_I_BB(J/m2)'].to_numpy()
      monthly_lunar_BB = df_dosage['monthly_lunar_I_BB(J/m2)'].to_numpy()
      monthly_ALAN_BB = df_dosage['monthly_ALAN_I_BB(J/m2)'].to_numpy()
      monthly_twilight_BB = df_dosage['monthly_twilight_I_BB(J/m2)'].to_numpy()

      monthly_solar_I_BB_datum = df_dosage['monthly_solar_I_BB_datum(J/m2)'].to_numpy()
      monthly_lunar_I_BB_datum = df_dosage['monthly_lunar_I_BB_datum(J/m2)'].to_numpy()
      monthly_ALAN_I_BB_datum = df_dosage['monthly_ALAN_I_BB_datum(J/m2)'].to_numpy()
      monthly_twilight_I_BB_datum = df_dosage['monthly_twilight_I_BB_datum(J/m2)'].to_numpy()

       
      jan_index = np.where(smonth == 'Jan')
      feb_index = np.where(smonth == 'Feb')
      mar_index = np.where(smonth == 'Mar')
      apr_index = np.where(smonth == 'Apr')
      may_index = np.where(smonth == 'May')
      jun_index = np.where(smonth == 'Jun')
      jul_index = np.where(smonth == 'Jul')
      aug_index = np.where(smonth == 'Aug')
      sep_index = np.where(smonth == 'Sep')
      oct_index = np.where(smonth == 'Oct')
      nov_index = np.where(smonth == 'Nov')
      dec_index = np.where(smonth == 'Dec')

      mean_max_lunar = []; std_max_lunar = []
      mean_max_ALAN = []; std_max_ALAN = []
      mean_max_lunar_datum = []; std_max_lunar_datum = []
      mean_max_ALAN_datum = []; std_max_ALAN_datum = []

      # irradiances
      # Jan
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[jan_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[jan_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[jan_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[jan_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[jan_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[jan_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[jan_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[jan_index]));

      # Feb
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[feb_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[feb_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[feb_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[feb_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[feb_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[feb_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[feb_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[feb_index]));

      # Mar
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[mar_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[mar_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[mar_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[mar_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[mar_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[mar_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[mar_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[mar_index]));

      # Apr
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[apr_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[apr_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[apr_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[apr_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[apr_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[apr_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[apr_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[apr_index]));

      # May
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[may_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[may_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[may_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[may_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[may_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[may_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[may_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[may_index]));

      # Jun
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[jun_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[jun_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[jun_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[jun_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[jun_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[jun_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[jun_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[jun_index]));

      # Jul
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[jul_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[jul_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[jul_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[jul_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[jul_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[jul_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[jul_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[jul_index]));

      # Aug
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[aug_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[aug_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[aug_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[aug_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[aug_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[aug_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[aug_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[aug_index]));

      # Sep
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[sep_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[sep_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[sep_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[sep_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[sep_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[sep_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[sep_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[sep_index]));

      # Oct
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[oct_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[oct_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[oct_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[oct_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[oct_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[oct_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[oct_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[oct_index]));

      # Nov
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[jan_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[jan_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[jan_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[jan_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[jan_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[jan_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[jan_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[jan_index]));

      # Dec
      mean_max_lunar = np.append(mean_max_lunar,np.mean(monthly_max_lunar_BB[dec_index])); 
      std_max_lunar = np.append(std_max_lunar,np.std(monthly_max_lunar_BB[dec_index]));
      mean_max_ALAN = np.append(mean_max_ALAN,np.mean(monthly_max_ALAN_BB[dec_index])); 
      std_max_ALAN = np.append(std_max_ALAN,np.std(monthly_max_ALAN_BB[dec_index]));
      mean_max_lunar_datum = np.append(mean_max_lunar_datum,np.mean(monthly_max_lunar_I_BB_datum[dec_index])); 
      std_max_lunar_datum = np.append(std_max_lunar_datum,np.std(monthly_max_lunar_I_BB_datum[dec_index]));
      mean_max_ALAN_datum = np.append(mean_max_ALAN_datum,np.mean(monthly_max_ALAN_I_BB_datum[dec_index])); 
      std_max_ALAN_datum = np.append(std_max_ALAN_datum,np.std(monthly_max_ALAN_I_BB_datum[dec_index]));

      # Dosages
      mean_lunar = []; std_lunar = []
      mean_ALAN = []; std_ALAN = []
      mean_lunar_datum = []; std_lunar_datum = []
      mean_ALAN_datum = []; std_ALAN_datum = []

      # Jan
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[jan_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[jan_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[jan_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[jan_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[jan_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[jan_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[jan_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[jan_index]));

      # Feb
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[feb_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[feb_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[feb_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[feb_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[feb_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[feb_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[feb_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[feb_index]));

      # Mar
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[mar_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[mar_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[mar_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[mar_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[mar_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[mar_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[mar_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[mar_index]));

      # Apr
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[apr_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[apr_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[apr_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[apr_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[apr_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[apr_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[apr_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[apr_index]));

      # May
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[may_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[may_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[may_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[may_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[may_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[may_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[may_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[may_index]));

      # Jun
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[jun_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[jun_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[jun_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[jun_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[jun_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[jun_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[jun_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[jun_index]));

      # Jul
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[jul_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[jul_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[jul_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[jul_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[jul_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[jul_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[jul_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[jul_index]));

      # Aug
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[aug_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[aug_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[aug_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[aug_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[aug_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[aug_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[aug_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[aug_index]));

      # Sep
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[sep_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[sep_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[sep_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[sep_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[sep_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[sep_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[sep_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[sep_index]));

      # Oct
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[oct_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[oct_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[oct_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[oct_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[oct_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[oct_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[oct_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[oct_index]));

      # Nov
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[jan_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[jan_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[jan_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[jan_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[jan_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[jan_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[jan_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[jan_index]));

      # Dec
      mean_lunar = np.append(mean_lunar,np.mean(monthly_lunar_BB[dec_index])); 
      std_lunar = np.append(std_lunar,np.std(monthly_lunar_BB[dec_index]));
      mean_ALAN = np.append(mean_ALAN,np.mean(monthly_ALAN_BB[dec_index])); 
      std_ALAN = np.append(std_ALAN,np.std(monthly_ALAN_BB[dec_index]));
      mean_lunar_datum = np.append(mean_lunar_datum,np.mean(monthly_lunar_I_BB_datum[dec_index])); 
      std_lunar_datum = np.append(std_lunar_datum,np.std(monthly_lunar_I_BB_datum[dec_index]));
      mean_ALAN_datum = np.append(mean_ALAN_datum,np.mean(monthly_ALAN_I_BB_datum[dec_index])); 
      std_ALAN_datum = np.append(std_ALAN_datum,np.std(monthly_ALAN_I_BB_datum[dec_index]));

      
      outdf = pd.DataFrame()
      outdf['mean_lunar_dosage(J/m2)'] = mean_lunar
      outdf['std_lunar_dosage(J/m2)'] = std_lunar
      outdf['mean_ALAN_dosage(J/m2)'] = mean_ALAN
      outdf['std_ALAN_dosage(J/m2)'] = std_ALAN 
      outdf['mean_lunar_dosage_datum(J/m2)'] = mean_lunar_datum
      outdf['std_lunar_dosage_datum(J/m2)'] = std_lunar_datum
      outdf['mean_ALAN_dosage_datum(J/m2)'] = mean_ALAN_datum
      outdf['std_ALAN_dosage_datum(J/m2)'] = std_ALAN_datum

      outdf['mean_max_lunar(uW/m2)'] = mean_max_lunar/1e-6
      outdf['std_max_lunar(uW/m2)'] = std_max_lunar/1e-6
      outdf['mean_max_ALAN(uW/m2)'] = mean_max_ALAN/1e-6
      outdf['std_max_ALAN(uW/m2)'] = std_max_ALAN/1e-6
      outdf['mean_max_lunar_datum(uW/m2)'] = mean_max_lunar_datum/1e-6
      outdf['std_max_lunar_datum(uW/m2)'] = std_max_lunar_datum/1e-6
      outdf['mean_max_ALAN_datum(uW/m2)'] = mean_max_ALAN_datum/1e-6
      outdf['std_max_ALAN_datum(uW/m2)'] = std_max_ALAN_datum/1e-6

      print("Writing out metonic variabilty file")
      print(outdir+'/metonic_variability_plymouth_2001-2020.csv')
      outdf.to_csv(outdir+'/metonic_variability_plymouth_2001-2020.csv', index=False)

   return 0

def main():
   aparser = argparse.ArgumentParser()
   aparser.add_argument("-loc", "--location", action="store", type=str, help="geographical location")
   aparser.add_argument("-p", "--plot", action="store_true", help="plots of the different light source contributions")
   aparser.add_argument("-m", "--max", action="store_true", help="maximum plots of the different light source contributions")
   aparser.add_argument("-pos", "--position", action="store", type=str, help="position < surface | intertidal >")
   aparser.add_argument("-me", "--metonic", action="store_true", help="calculation of variability over 20 years")
   args = aparser.parse_args()
   
   if args.position:
      print("Plotting City seasonal intercomparison")
      print(" "+args.position)
      dosage_seasons(args.position)
      max_seasons(args.position)
      sys.exit(0)
   
   if args.plot and args.location:
      print("Plotting option set")
      if args.max:
         print("Plotting option set for maximum")
         plot_flag = max_plot(args.location)
      else:
         plot_flag = dosage_plot(args.location)
      sys.exit(0)
   
   if args.location:
      dosage_opf(args.location)
      sys.exit(0)
   
   if args.metonic:
      metonic_max_plot()
      sys.exit(0)
   
   return

# Run script if called from command line.   
if __name__=='__main__':
    main()

