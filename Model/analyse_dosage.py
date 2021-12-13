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
       
      # Output dataframe
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
      
      outdf.to_csv(outfname, index=False)
      print("============================================")
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
      axs[0,0].set_title('Surface: Solar'); axs[0,0].set(ylabel='Dosage [MJ]'); #axs[0,0].set(xlabel='Month')
      axs[0,0].set_ylim(top = 500, bottom = 0);
     
      pcm = axs[0,1].plot(smonth, monthly_twilight_BB/1e+3, 'k', label='Broadband')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Red/1e+3, 'r', label='Red(620-740nm)')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Green/1e+3, 'g', label='Green(495-560nm)')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Blue/1e+3, 'b', label='Blue(400-500nm)')
      #pcm = axs[0,1].plot(smonth, (monthly_twilight_I_Blue+monthly_twilight_I_Red+monthly_twilight_I_Green)/1e+3, 'gray')
      axs[0,1].set_title('Surface: Twilight'); axs[0,1].set(ylabel='Dosage [kJ]'); #axs[0,1].set(xlabel='Month')
      axs[0,1].legend()
      axs[0,1].set_ylim(top = 500, bottom = 0);
      
      pcm = axs[1,0].plot(smonth, monthly_lunar_BB, 'k', label='Broadband')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Red, 'r', label='Red(620-740nm)')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Green, 'g', label='Green(495-560nm)')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Blue, 'b', label='Blue(400-500nm)')
      #pcm = axs[1,0].plot(smonth, (monthly_lunar_I_Blue+monthly_lunar_I_Red+monthly_lunar_I_Green), 'gray')
      axs[1,0].set_title('Surface: Lunar'); axs[1,0].set(ylabel='Dosage [J]'); axs[1,0].set(xlabel='Month')
      axs[1,0].set_ylim(top = 500, bottom = 0);
   
      pcm = axs[1,1].plot(smonth, monthly_ALAN_BB, 'k', label='Broadband')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Red, 'r', label='Red(620-740nm)')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Green, 'g', label='Green(495-560nm)')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Blue, 'b', label='Blue(400-500nm)')
      #pcm = axs[1,1].plot(smonth, (monthly_ALAN_I_Blue+monthly_ALAN_I_Red+monthly_ALAN_I_Green), 'gray')
      axs[1,1].set_title('Surface: ALAN'); axs[1,1].set(ylabel='Dosage [J]'); axs[1,1].set(xlabel='Month')
      axs[1,1].set_ylim(top = 500, bottom = 0);
   
      fig.savefig(outdir+'/Dosage_Surface_'+location+'_'+syear+'.png')
      plt.close(fig)
      
      fig, axs = plt.subplots(2,2, figsize=(10,7))
      fig.suptitle(location+" "+syear, fontsize=16)

      pcm = axs[0,0].plot(smonth, monthly_solar_I_BB_datum/1e+6, 'k', label='Broadband')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Red_datum/1e+6, 'r', label='Red(620-740nm)')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Green_datum/1e+6, 'g', label='Green(495-560nm)')
      pcm = axs[0,0].plot(smonth, monthly_solar_I_Blue_datum/1e+6, 'b', label='Blue(400-500nm)')
      axs[0,0].set_title('Datum: Solar'); axs[0,0].set(ylabel='Dosage [MJ]'); #axs[0,0].set(xlabel='Month')
      axs[0,0].set_ylim(top = 500, bottom = 0);
     
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_BB_datum/1e+3, 'k', label='Broadband')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Red_datum/1e+3, 'r', label='Red(620-740nm)')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Green_datum/1e+3, 'g', label='Green(495-560nm)')
      pcm = axs[0,1].plot(smonth, monthly_twilight_I_Blue_datum/1e+3, 'b', label='Blue(400-500nm)')
      axs[0,1].set_title('Datum: Twilight'); axs[0,1].set(ylabel='Dosage [kJ]'); #axs[0,1].set(xlabel='Month')
      axs[0,1].legend()
      axs[0,1].set_ylim(top = 500, bottom = 0);
      
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_BB_datum, 'k', label='Broadband')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Red_datum, 'r', label='Red(620-740nm)')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Green_datum, 'g', label='Green(495-560nm)')
      pcm = axs[1,0].plot(smonth, monthly_lunar_I_Blue_datum, 'b', label='Blue(400-500nm)')
      axs[1,0].set_title('Datum: Lunar'); axs[1,0].set(ylabel='Dosage [J]'); axs[1,0].set(xlabel='Month')
      axs[1,0].set_ylim(top = 500, bottom = 0);
   
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_BB_datum, 'k', label='Broadband')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Red_datum, 'r', label='Red(620-740nm)')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Green_datum, 'g', label='Green(495-560nm)')
      pcm = axs[1,1].plot(smonth, monthly_ALAN_I_Blue_datum, 'b', label='Blue(400-500nm)')
      axs[1,1].set_title('Datum: ALAN'); axs[1,1].set(ylabel='Dosage [J]'); axs[1,1].set(xlabel='Month')
      axs[1,1].set_ylim(top = 500, bottom = 0);
   
      fig.savefig(outdir+'/Dosage_Datum_'+location+'_'+syear+'.png')
      plt.close(fig)
      
      
      #pdb.set_trace()   

   return 0


def main():
   aparser = argparse.ArgumentParser()
   aparser.add_argument("-loc", "--location", action="store", type=str, help="geographical location")
   aparser.add_argument("-p", "--plot", action="store_true", help="plots of the different light source contributions")
   args = aparser.parse_args()
   
   if args.location:
      dosage_opf(args.location)
   
   if args.plot:
      print("Plotting option set")
      plot_flag = dosage_plot(args.location)
   
   
   
   return

# Run script if called from command line.   
if __name__=='__main__':
    main()

