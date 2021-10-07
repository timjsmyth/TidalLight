#!/usr/bin/env python
import os
import sys
from timeit import default_timer as timer
import numpy as np
import netCDF4 as nc
from PIL import Image
from itertools import product
from random import randrange, uniform
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import RegularGridInterpolator
from PIL import Image
from osgeo import gdal
from numpy import zeros, newaxis
from scipy.ndimage import zoom
from matplotlib.colors import LogNorm
# from mpl_toolkits.basemap import Basemap
import argparse
import pandas as pd

def find_nearest(array, value):
   array = np.asarray(array)
   idx = (np.abs(array - value)).argmin()
   val = array[idx]
   return val, idx

def find_next_nearest_above(array, value):
   array = np.asarray(array)
   idx = (np.abs(array - value)).argmin()
   above_idx = idx + 1 
   val_above = array[above_idx]
   return val_above, above_idx

def find_next_nearest_below(array, value):
   array = np.asarray(array)
   idx = (np.abs(array - value)).argmin()
   below_idx = idx - 1
   val_below = array[below_idx]
   return val_below, below_idx

def check_val(kd_val, lat, LatIndice, lon, LonIndice):
   kd_ = kd_val[LatIndice, LonIndice]
   
   latlon = [lat, LatIndice, lon, LonIndice]
   if isinstance(kd_, float) == True:
      # lat_indice = LatIndice; lon_indice = LonIndice
      # ext_lats = lat; ext_lons = lon
      kd = float(kd_)
      print("found one at",lat, lon,
            kd)
      return kd, latlon
   else:
      kd_ = []
      return kd_, latlon
      

def read_falchi(latitude, longitude, sky_condition, input_flag):
   
   if sky_condition=="cloudy":
      # "CLOUDY" sky
      m_blue = 17.97854676
      m_green = 35.80029137
      m_red = 29.40875316
      c_blue = 57.22738163
      c_green = 18.58249682
      c_red = 23.73604658
   else:
      # "CLEAR" sky
      m_blue = 4.530851872
      m_green = 7.272195079
      m_red = 6.372282777
      c_blue = 59.57581038
      c_green = 25.091669
      c_red = 26.19770495

   # From Batnes et al. (2015) - table 5
   thresh_irr_total_uW_m2 = 0.102 # in uW/m2 converted from uE/m2/s  

   # This needs changing depending on where the 12.8GB of HYDROLIGHT data files and Falchi Atlas are stored...
   #datadir = '../../ALAN_Map+Kd/' # Relative path to directory that contains Kd files 
   datadir = '/home/scratch/data/workspace/ALAN/Kd' # Relative path to directory that contains Kd files 

   # FALCHI map runs between 85N -> 60S (145 degrees)
   Image.MAX_IMAGE_PIXELS = 751939200
   im = gdal.Open(datadir + "/../World_Atlas_2015.tif")
   falchi_map = im.ReadAsArray()
   width = im.RasterXSize
   height = im.RasterYSize
   gt = im.GetGeoTransform()
   minx = gt[0]
   miny = gt[3] + width*gt[4] + height*gt[5] 
   maxx = gt[0] + width*gt[1] + height*gt[2]
   maxy = gt[3]
   
   falchi_lats=np.zeros(int(height))*0.
   falchi_lons=np.zeros(int(width))*0.

   falchi_lats_orig = falchi_lats
   falchi_lons_orig = falchi_lons

   for ii in range(width):
      falchi_lons[ii] = ii*gt[1]+gt[0]

   for jj in range(height):
      falchi_lats[jj] = gt[3]+jj*gt[5]
   
   ##### Find parameters of input location on Falchi map for checking against Kd locations later in the script #####
   falchi_lat_pos, falchi_lat_indice = find_nearest(falchi_lats, latitude) # falchi_lat_pos = position in degrees nearest to the input value "latitude" on the falchi map; falchi_lat_indice = row index of falchi_lat_pos in falchi_lats array 
   falchi_lon_pos, falchi_lon_indice = find_nearest(falchi_lons, longitude) # falchi_lon_pos = position in degrees nearest to the input value "longitude" on the falchi map; falchi_lon_indice = row index of falchi_lon_pos in falchi_lons array
   # print(falchi_map.shape, type(falchi_map))

   ##### Create header lists and DataFrames #####
   column_r = ["Lat", "Lon", "Month", "ALAN_R (uW/m^2)", "Kd_R"]
   column_g = ["Lat", "Lon", "Month", "ALAN_G (uW/m^2)", "Kd_G"]
   column_b = ["Lat", "Lon", "Month", "ALAN_B (uW/m^2)", "Kd_B"]
   df_r = pd.DataFrame(columns=column_r)
   df_g = pd.DataFrame(columns=column_g)
   df_b = pd.DataFrame(columns=column_b)
   
   
   
   # IF Working from home using external harddrive as data store (mount D: drive int the "mnt" folder using sudo mkdir /mnt/d; sudo mount -t drvfs D: /mnt/d)
   # datadir = "/mnt/d/PML_LastDay/ALICE_Project/ALAN_Map+Kd/"
   month=0 # Create a counter variable to represent Months in conjunction with the Kd files
   for filename in sorted(os.listdir(datadir)): # Loop through all files in the "datadir"
      if filename.endswith(".nc"): # Only loop through NetCDF4 files
         ##### Determine the filename for the netCDF and read in the file ######
         fname = os.path.join(datadir, filename)
         month=month+1 # Month counter
         infile=nc.Dataset(fname, 'r') # read NetCDF4 files

         ##### Split NetCDF4 file into variables according to headers #####
         kd_blue=infile.variables['kd_blue'][:] 
         kd_green=infile.variables['kd_green'][:]
         kd_red=infile.variables['kd_red'][:]
         nc_dims = [dim for dim in infile.dimensions]  # list of nc dimensions # NOT Relevant 
         time=infile.variables['time'][:]
         lats=infile.variables['lat'][:]
         lons=infile.variables['lon'][:]
         
         lats_orig = lats # Store lats parameter for safe keeping # NOT Relevant 
         lons_orig = lons # Store lats parameter for safe keeping # NOT Relevant
         ext_lats, lat_indices = find_nearest(lats, latitude) # Extracted latitude in degrees and row index (ext_lats, lat_indices) nearest to the input lat within the latitudes of the Kd input file (lats)
         ext_lons, lon_indices = find_nearest(lons, longitude) # Extracted longitude in degrees and row index (ext_lons, lon_indices) nearest to the input lon within the longitudes of the Kd input file (lons)
         # print("lat/lon nearest to input pos in Kd coords \n", ext_lats, ext_lons)

         ##### Find Kd for given Latitude and Longitude #####
         ## BLUE ##
         kd_blue_ext = kd_blue[0][lat_indices, lon_indices]
         # print("Blue Kd\n", kd_blue_ext)
         kd_blue_shape = kd_blue_ext.shape
         ## GREEN ## 
         kd_green_ext = kd_green[0][lat_indices, lon_indices]
         # print("Green Kd\n", kd_green_ext)
         kd_green_shape = kd_green_ext.shape
         ## RED ##
         kd_red_ext = kd_red[0][lat_indices, lon_indices]

         kd_red_shape = kd_red_ext.shape
         # print("Red Kd\n", kd_red_ext)
         empty = np.ma.masked_all(kd_red_shape)

         if type(kd_red_ext) != float:
            Next_lat = ext_lats; Next_lon = ext_lons
            next_lat = ext_lats; next_lon = ext_lons
            count = 0
            while type(kd_red_ext) != float:
               count = count+1
               # print(count)
               Next_lat, Nlat_indice = find_next_nearest_above(lats, Next_lat)
               Next_lon, Nlon_indice = find_next_nearest_above(lons, Next_lon)
               # print('above', Next_lat, Next_lon)
               next_lat, nlat_indice = find_next_nearest_below(lats, next_lat)
               next_lon, nlon_indice = find_next_nearest_below(lons, next_lon)
               # print('below', next_lat, next_lon)
               kd_red_ext, latlon = check_val(kd_red[0], Next_lat, Nlat_indice, Next_lon, Nlon_indice)
               if type(kd_red_ext) == float:
                  # print("exiting loop")
                  break
               kd_red_ext, latlon = check_val(kd_red[0], next_lat, nlat_indice, next_lon, nlon_indice)
               if type(kd_red_ext) == float:
                  # print("exiting loop")
                  break
               kd_red_ext, latlon = check_val(kd_red[0], Next_lat, Nlat_indice, next_lon, nlon_indice)
               if type(kd_red_ext) == float:
                  # print("exiting loop")
                  break
               kd_red_ext, latlon = check_val(kd_red[0], next_lat, nlat_indice, Next_lon, Nlon_indice)
               if type(kd_red_ext) == float:
                  # print("exiting loop")
                  break
               elif count > 10:
                  break
               else:
                  continue
               # kd_red_ext = kd_red[0][Nlat_indice, Nlon_indice]
               # if isinstance(kd_red_ext, float) == True:
               #    lat_indice = Nlat_indice; lon_indice = Nlon_indice
               #    ext_lats = Next_lat; ext_lons = Next_lon
               #    print(Next_lat, Next_lon, "In the loop", kd_red_ext)
               #    break
               # kd_red_ext = kd_red[0][Nlat_indice, nlon_indice]
               # if isinstance(kd_red_ext, float) == True:
               #    lat_indice = Nlat_indice; lon_indice = Nlon_indice
               #    ext_lats = Next_lat; ext_lons = Next_lon
               #    print(Next_lat, Next_lon, "In the loop", kd_red_ext)
               #    break
               
               # kd_red_ext = kd_red[0][nlat_indice, nlon_indice]
               # print("below", next_lat, next_lon)
               # if isinstance(kd_red_ext, float) == True:
               #    lat_indice = nlat_indice; lon_indice = nlon_indice
               #    ext_lats = next_lat; ext_lons = next_lon
               #    print(next_lat, next_lon, "In the loop", kd_red_ext)
               #    break
            ext_lats = latlon[0]; ext_lons = latlon[2]; lat_indice = latlon[1]; lon_indice = latlon[3]      
               
            ## BLUE ##
            kd_blue_ext = kd_blue[0][lat_indice, lon_indice]
            kd_blue_shape = kd_blue_ext.shape
            ## GREEN ## 
            kd_green_ext = kd_green[0][lat_indice, lon_indice]
            kd_green_shape = kd_green_ext.shape
            ## RED ##
            kd_red_ext = kd_red[0][lat_indice, lon_indice]
            kd_red_shape = kd_red_ext.shape
            

         ##### Determine the nearest position to the extracted Kd position on the falchi map #####
         Kd_lat_pos, Kd_lat_pos_idx = find_nearest(falchi_lats, ext_lats) # Nearest position & row index on the falchi map to the value of the extracted latitude from the Kd file parameter 'lats'
         Kd_lon_pos, Kd_lon_pos_idx = find_nearest(falchi_lons, ext_lons) # Nearest position & row index on the falchi map to the value of the extracted longitude from the Kd file parameter 'lons'
         # print("lat/Lon nearest to Falchi pos in Kd coords \n", Kd_lat_pos, Kd_lon_pos)
         # print("lat/lon nearest to input pos in Kd coords \n", ext_lats, ext_lons)

         Kd_pos_ALAN = falchi_map[Kd_lat_pos_idx, Kd_lon_pos_idx]
         # print("Kd pos ALAN value from Falchi map\n", Kd_pos_ALAN)
         input_pos_ALAN = falchi_map[falchi_lat_indice, falchi_lon_indice]
         # print("Input pos ALAN value from Falchi map\n", input_pos_ALAN)
         if month < 2:
            if Kd_pos_ALAN != input_pos_ALAN:
               print("WARNING:\n    Position on Falchi map is not equal for coordinates found in Kd files.\n     Review input latitude and longitude.")
               print("Kd Latitude:", Kd_lat_pos, "      Falchi Latitude:", falchi_lat_pos)
               print("Kd Longitude:", Kd_lon_pos, "     Falchi Longitude:", falchi_lon_pos)
               print("Kd ALAN (uCd/m^2):", Kd_pos_ALAN, "      Input ALAN (uCd/m^2):", input_pos_ALAN) 
               if input_flag == 1:
                  to_be_or_not_to_be = input("Do you wish to continue despite this difference?\n"
                                             "  Y: yes continue, I am aware the resulting Kd and ALAN values will not originate from the same location\n" 
                                             "  n: no, stop! I must change my inputs immediately\n")
                  if to_be_or_not_to_be == "n":
                     exit(1)
               else: 
                  pass
            else:
               # print("Kd location matches input location on Falchi map - Hooray")
               print("Kd ALAN (uCd/m^2):", Kd_pos_ALAN, "      Input ALAN (uCd/m^2):", input_pos_ALAN) 

         ######### Calculations ###################
         print("USING REGRESSION INTERCEPT")
         ## Calculate surface irradiance in uW/m^2 ##
         sfce_irr_blue_uW_m2 = m_blue*input_pos_ALAN + c_blue
         # print("Blue Irr uW/m2", sfce_irr_blue_uW_m2, "Kd_Blue", kd_blue_ext)
         # print("Surface Irradiance: Blue \n", sfce_irr_blue_uW_m2)
         sfce_irr_green_uW_m2 = m_green*input_pos_ALAN + c_green
         # print("Green Irr uW/m2", sfce_irr_green_uW_m2, "Kd_Green", kd_green_ext)
         # print("Surface Irradiance: Green \n", sfce_irr_green_uW_m2)
         sfce_irr_red_uW_m2 = m_red*input_pos_ALAN + c_red
         # print("Red Irr uW/m2", sfce_irr_red_uW_m2, "Kd_Red", kd_red_ext)
         # print("Surface Irradiance: Red \n", sfce_irr_blue_uW_m2)
         sfce_irr_total_uW_m2 = sfce_irr_blue_uW_m2 + sfce_irr_green_uW_m2 + sfce_irr_red_uW_m2
         
         ALAN_mCdm2 = input_pos_ALAN
         input_pos_ALAN = sfce_irr_total_uW_m2
         print("ALAN in uW/m^2 = ", input_pos_ALAN)
         R_Red = sfce_irr_red_uW_m2/sfce_irr_total_uW_m2
         ALAN_Red = R_Red*input_pos_ALAN # ALAN in uW/m^2
         R_Green = sfce_irr_green_uW_m2/sfce_irr_total_uW_m2
         ALAN_Green = R_Green*input_pos_ALAN # ALAN in uW/m^2
         R_Blue = sfce_irr_blue_uW_m2/sfce_irr_total_uW_m2
         ALAN_Blue = R_Blue*input_pos_ALAN # ALAN in uW/m^2
         ###### Collate and append to DataFrames ######
         _r = [falchi_lat_pos, falchi_lon_pos, month, ALAN_Red, kd_red_ext]
         _g = [falchi_lat_pos, falchi_lon_pos, month, ALAN_Green, kd_green_ext]
         _b = [falchi_lat_pos, falchi_lon_pos, month, ALAN_Blue, kd_blue_ext]
   
         df_r = df_r.append(pd.Series(_r, index=df_r.columns), ignore_index=True)
         df_g = df_g.append(pd.Series(_g, index=df_g.columns), ignore_index=True)
         df_b = df_b.append(pd.Series(_b, index=df_b.columns), ignore_index=True)
         print("Month: ", month)
         ALAN_total = pd.Series([input_pos_ALAN])
         ALAN_mCd = pd.Series([ALAN_mCdm2])
   # print("Total ALAN", sfce_irr_total_uW_m2)
   return df_r, df_g, df_b, ALAN_total, ALAN_mCd

if __name__=='__main__':
   parser = argparse.ArgumentParser(description='Produce ALAN maps')
   parser.add_argument("-sky",  type=str, action='store',   help="Input either [cloudy], [clear] for scaling factors")
   parser.add_argument("-lat", "--latitude",type=float, action="store", help="Longitude in decimal degrees")
   parser.add_argument("-lon", "--longitude", type=float, action="store", help="Latitude in decimal degrees")
   args = parser.parse_args() 
   

   latitude = 50.27 # latitude_deg
   longitude = -4.13 # longitude_deg
   sky = "clear"
   if args.sky:
      sky = args.sky
   if args.latitude:
      latitude = args.latitude
   else:
      print("No location input, default location: Plymouth, United Kingdom")
   if args.longitude:
      longitude = args.longitude
   
   Red, Green, Blue, ALAN_Total, ALAN_mCd = read_falchi(latitude, longitude, sky)
   print("Red\n", Red)
   print("Green\n", Green)
   print("Blue\n", Blue)
   
