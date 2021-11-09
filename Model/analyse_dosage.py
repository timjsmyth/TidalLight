#!/bin/python3
import os
import glob
import numpy as np
import pandas as pd
import pdb
import calendar
import argparse

uW_W = 1.0e-6 # scaling from uW to W
indir = "/users/rsg/tjsm/ALAN/tidal/TidalLight/Model/Output/annual"
outdir = "/users/rsg/tjsm/ALAN/tidal/TidalLight/Model/Output/annual/analysis"

def main():
   aparser = argparse.ArgumentParser()
   aparser.add_argument("-loc", "--location", action="store", type=str, help="geographical location")
   args = aparser.parse_args()
   
   location = args.location
   print("Location: ", location)

   for fname in glob.glob(indir+'/*'+location+'*.csv'):
      #print(fname)
      if "Lunar_" in fname:
         print("Lunar file", fname)
         df_lunar_irradiances = pd.read_csv(fname)
      if "Solar_" in fname:
         print(" Solar file", fname)
         df_solar_irradiances = pd.read_csv(fname)
      if "ALAN_" in fname:
         print("  ALAN file", fname)
         df_ALAN_irradiances = pd.read_csv(fname)
   
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
      if (m == 12):
         # Final run through is total of year
         print("Summing between days: ", mm[0],mm[12])
         mm_index = np.where((df_solar_irradiances['Jday(decimal)'] >= mm[0]) & (df_solar_irradiances['Jday(decimal)'] <= mm[12]))
         mm_twilight_index = np.where((df_solar_irradiances['Jday(decimal)'] >= mm[0]) & (df_solar_irradiances['Jday(decimal)'] <= mm[12]) & (df_solar_irradiances['night(0)_day(1)'] > 0) & (df_solar_irradiances['night(0)_day(1)'] < 1.0))
      else:
         print("Summing between days: ", mm[m],mm[m+1]-1)
         mm_index = np.where((df_solar_irradiances['Jday(decimal)'] >= mm[m]) & (df_solar_irradiances['Jday(decimal)'] <= mm[m+1])) 
         mm_twilight_index = np.where((df_solar_irradiances['Jday(decimal)'] >= mm[m]) & (df_solar_irradiances['Jday(decimal)'] <= mm[m+1]) & (df_solar_irradiances['night(0)_day(1)'] > 0) & (df_solar_irradiances['night(0)_day(1)'] < 1.0))
      # Broadband
      monthly_solar_I_BB.append(np.sum(solar_I_BB[mm_index[0]]))
      monthly_lunar_I_BB.append(np.sum(lunar_I_BB[mm_index[0]]))
      monthly_ALAN_I_BB.append(np.sum(ALAN_I_BB[mm_index[0]]))
      monthly_twilight_I_BB.append(np.sum(solar_I_BB[mm_twilight_index[0]]))
      monthly_solar_I_BB_datum.append(np.sum(solar_I_BB_datum[mm_index[0]]))
      monthly_lunar_I_BB_datum.append(np.sum(lunar_I_BB_datum[mm_index[0]]))
      monthly_ALAN_I_BB_datum.append(np.sum(ALAN_I_BB_datum[mm_index[0]]))
      monthly_twilight_I_BB_datum.append(np.sum(solar_I_BB_datum[mm_twilight_index[0]]))
      
      # Spectral
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


# Run script if called from command line.   
if __name__=='__main__':
    main()

