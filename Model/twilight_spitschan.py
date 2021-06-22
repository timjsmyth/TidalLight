#!/bin/python3
import os
import sys
import numpy as np
import pandas as pd
import argparse
import pdb
import warnings
warnings.simplefilter('ignore', np.RankWarning)

def main(solar_elevation, atmosphere_type, spectral_group):
      
   # 1. Read in hyperspectral intensity as a function of solar zenith angle
   #    Read these in from the Spitzchan et al. (2016) Sci. Rep. doi:10.1038/srep26756 outputs
   #    LUT/Rural_Twilight.csv
   #    LUT/Urban_Twilight.csv
   #    Intensities for Civil (0 - 6 degrees); Nautical (6 - 12 degrees) and 
   #    Astronomical (12 - 18 degrees) below horizon twilights
   #    Each of these distinct parts of the Twilight will have their own spectral functions
   if atmosphere_type == "rural":
      df = pd.read_csv("Required/LUT/Rural_Twilight.csv", skiprows=1)
   if atmosphere_type == "urban":
      df = pd.read_csv("Required/LUT/Urban_Twilight.csv", skiprows=1)

   df.rename(columns={'Elev=0.0(log(W/m2/nm))':'0.0',
                      'Elev=-1.0(log(W/m2/nm))':'-1.0',
                      'Elev=-2.0(log(W/m2/nm))':'-2.0',
                      'Elev=-3.0(log(W/m2/nm))':'-3.0',
                      'Elev=-4.0(log(W/m2/nm))':'-4.0',
                      'Elev=-5.0(log(W/m2/nm))':'-5.0',
                      'Elev=-6.0(log(W/m2/nm))':'-6.0',
                      'Elev=-7.0(log(W/m2/nm))':'-7.0',
                      'Elev=-8.0(log(W/m2/nm))':'-8.0',
                      'Elev=-9.0(log(W/m2/nm))':'-9.0',
                      'Elev=-10.0(log(W/m2/nm))':'-10.0',
                      'Elev=-11.0(log(W/m2/nm))':'-11.0',
                      'Elev=-12.0(log(W/m2/nm))':'-12.0',
                      'Elev=-13.0(log(W/m2/nm))':'-13.0',
                      'Elev=-14.0(log(W/m2/nm))':'-14.0',
                      'Elev=-15.0(log(W/m2/nm))':'-15.0',
                      'Elev=-16.0(log(W/m2/nm))':'-16.0',
                      'Elev=-17.0(log(W/m2/nm))':'-17.0',
                      'Elev=-18.0(log(W/m2/nm))':'-18.0'},
                      inplace=True)   

   # calculate the spectral intensity at the required angle
   # Find which two columns the angle is between.
   elevation_angles = len(df.columns) - 1      
   elevation_angles = np.array(-(np.arange(elevation_angles)))
   
   out_df = pd.DataFrame()

   if (solar_elevation >= np.min(elevation_angles) and solar_elevation <= np.max(elevation_angles)):
      print("Sun is below horizon and twilight spectrum can be determined")
      # determine which LUT elevation angles the spectrum lies between
      index = np.where(np.float(solar_elevation) == elevation_angles)

      # Need to interpolate between two angles within the LUT (integer values between 0 and -18)
      if not index[0] and np.float(solar_elevation) != 0:
         indexu = np.max(np.where(np.float(solar_elevation) < elevation_angles))
         indexl = np.min(np.where(np.float(solar_elevation) > elevation_angles))
         
         # extract the two columns needed for the interpolation
         upper_spectrum = df.iloc[:,indexu+1].to_numpy() # bug fix within the df to add index+1 [wavelength is column 0]
         lower_spectrum = df.iloc[:,indexl+1].to_numpy() # bug fix within the df to add index+1 [wavelength is column 0]
         
         # setup empty array for the interpolated values
         actual_spectrum = np.zeros(len(lower_spectrum))
         
         # loop over the length of the spectrum for interpolation
         for i in range(len(actual_spectrum)):
            x = [elevation_angles[indexl], elevation_angles[indexu]]
            y = [lower_spectrum[i], upper_spectrum[i]]
            actual_spectrum[i] = np.interp(np.float(solar_elevation),x,y) 

      # No need for interpolation within the LUT if integer value of solar_elevation (0 -> -18) specified
      else:
         actual_spectrum = df.iloc[:,index[0]+1].to_numpy() # bug fix within the df to add index+1 [wavelength is column 0]
      
      # 3. Output as hyper, multi or broadband spectral (RGB)
      # i) hyperspectral
      if spectral_group == "hyperspectral":
         out_df['Start_Wavelength(nm)'] = df['Wavelength(nm)'] - 0.5
         out_df['End_Wavelength(nm)'] = df['Wavelength(nm)'] + 0.5
         out_df['Twilight_spectrum(log(W/m2/nm))'] = actual_spectrum
         out_df['Elevation(degrees)'] = solar_elevation
   
      # ii) broadband (400 - 700 nm)
      if spectral_group == "broadband":
         wavelength = df['Wavelength(nm)'].to_numpy()
         wavelength_range = np.where((df['Wavelength(nm)'] >= 400) & (df['Wavelength(nm)'] <= 700))
         wavelength = wavelength[wavelength_range[0]]
         actual_spectrum = actual_spectrum[wavelength_range[0]]
         broadband_integral = np.log10(np.sum(10**(actual_spectrum)))
         out_df['Start_Wavelength(nm)'] = [np.min(wavelength)]
         out_df['End_Wavelength(nm)'] = [np.max(wavelength)]
         out_df['Twilight_broadband(log(W/m2))'] = [broadband_integral]
         out_df['Elevation(degrees)'] = [solar_elevation]
         
      # iii) Skye Meter wavelengths Blue (400 - 500); Green (495 - 560nm); Red (620 - 740)
      if spectral_group == "skye":
         wavelength = df['Wavelength(nm)'].to_numpy()
         
         b_wavelength_range = np.where((df['Wavelength(nm)'] >= 400) & (df['Wavelength(nm)'] <= 500))
         b_wavelength = wavelength[b_wavelength_range[0]]
         b_actual_spectrum = actual_spectrum[b_wavelength_range[0]]
         b_integral = np.log10(np.sum(10**(b_actual_spectrum)))
         
         g_wavelength_range = np.where((df['Wavelength(nm)'] >= 495) & (df['Wavelength(nm)'] <= 560))
         g_wavelength = wavelength[g_wavelength_range[0]]
         g_actual_spectrum = actual_spectrum[g_wavelength_range[0]]
         g_integral = np.log10(np.sum(10**(g_actual_spectrum)))
   
         r_wavelength_range = np.where((df['Wavelength(nm)'] >= 620) & (df['Wavelength(nm)'] <= 740))
         r_wavelength = wavelength[r_wavelength_range[0]]
         r_actual_spectrum = actual_spectrum[r_wavelength_range[0]]
         r_integral = np.log10(np.sum(10**(r_actual_spectrum)))
         
         out_df['Start_Wavelength(nm)'] = [np.min(r_wavelength),np.min(g_wavelength),np.min(b_wavelength)]
         out_df['End_Wavelength(nm)'] = [np.max(r_wavelength),np.max(g_wavelength),np.max(b_wavelength),]
         out_df['Twilight_broadband(log(W/m2))'] = [r_integral,b_integral,g_integral]
         out_df['Elevation(degrees)'] = solar_elevation
         print(out_df) 
   return out_df

if __name__=='__main__':
   parser = argparse.ArgumentParser(
   description=__doc__,
   formatter_class=argparse.RawDescriptionHelpFormatter
   )

   # Require the solar elevation to be specified (between 0 and -18 degrees)
   parser.add_argument("-e", "--solar_elevation", type=float, required=True, help="--solar_elevation <degrees>")

   # Require either a rural or urban atmosphere to be specified
   atm_group = parser.add_mutually_exclusive_group(required=True)
   atm_group.add_argument("-u", "--urban", action='store_true', help="Read in Urban Twilight spectrum from LUT/Urban_Twilight.csv")
   atm_group.add_argument("-r", "--rural", action='store_true', help="Read in Rural Twilight spectrum from LUT/Rural_Twilight.csv")

   # Require one of hyperspectral, Skye or broadband spectrum options for output
   spec_group = parser.add_mutually_exclusive_group(required=True)
   spec_group.add_argument("-hy", "--hyperspectral", action='store_true', default=False, help="Return Dataframe of Hyperspectral irradiances")
   spec_group.add_argument("-s", "--skye", action='store_true', default=False, help="Return Dataframe of Skye meter irradiances (W/m2)")
   spec_group.add_argument("-b", "--broadband", action='store_true', default=False, help="Return Dataframe of broadband (PAR) irradiances (W/m2)")

   # Output filename (optional)
   parser.add_argument("-o", "--ofile", type=str, default='', help="--ofile <Output Filename>")

   args = parser.parse_args()
   if ((args.solar_elevation <= 0) and (args.solar_elevation >= -18.0)):
         solar_elevation = args.solar_elevation
   else:
      print("Solar elevation must be between 0 and -18 degrees")
      # sys.exit(0)
   main(args.solar_elevation, args.atm_group, args.spec_group)
   # Write out dataframe as a csv
   if args.ofile:
      filename = args.ofile
      print("Output filename: ", filename)
      out_df.to_csv(filename,index=False)
