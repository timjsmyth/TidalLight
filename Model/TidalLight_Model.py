#!/usr/bin/env python3

# Authour: Adam Wright
# Purpose: Calculate intensity of lunar, solar and ALAN brightness at the seasurface and seabed
# This script only contains a broadband Irradiance model and does not split spectrally
# a number of assumptions are made in the calculation process such as: Both the Sun and Moon are point sources, Lunar albedo is constant throughout cycle.
# Sky condition is constant for the entirety of the modelled time, scattering is constant through the atmosphere
# twilight is not factored in.
# Tidal data must be in the form of 'BODCTideData' .csv file 

# Import packages
from pysolar.solar import get_altitude 
from pysolar.radiation import get_air_mass_ratio
import utide
import matplotlib.dates as mdates
from mpl_toolkits.mplot3d import Axes3D

import get_TidalCoef
import Solplot_RGB
import Lunplot_RGB
import Aplot_RGB
import Falchi_Kd_Position
import twilight_spitschan
import gcirrad

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
warnings.simplefilter(action='ignore', category=FutureWarning)
import spectres

from tpxo_tide_prediction import (
    read_parameter_file,
    tide_predict,
    write_tides
    )

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w') # prevent function output from appearing in terminal

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__ # resume function output

def daterange(date_start, date_end):
    for n in range(int((date_end - date_start).days)):
        yield date_start + datetime.timedelta(n)

def get_air_mass_kasten_young(altitude_deg):
    # Original C
    # return 1.0 / (sin(altitude_deg * DtoR) + 0.50572 * pow(altitude + 6.07995, -1.6364));
    return 1.0 / (np.sin(np.radians(altitude_deg)) + 0.50572 * (altitude_deg + 6.07995)**(-1.6364))

# Lumme-Bowell (1981) 
# phase_angle in radians 
def lumme_bowell(phase_angle):
    Q = 0.108    # multiple-scattering fraction 

    #TO DO: find more precise angle for smooth curve? 
    if (np.degrees(phase_angle) < 25.07):
       Phi1 = 1.0 - np.sin(phase_angle) / (0.124 + 1.407 * np.sin(phase_angle) - 0.758 * np.sin(phase_angle) * np.sin(phase_angle))
    else:
       Phi1 = np.exp(-3.343 * (np.tan(phase_angle / 2.0))**0.632)

    # multiple scattering 
    Phi_m = 1.0 / np.pi * (np.sin(phase_angle) + (np.pi - phase_angle) * np.cos(phase_angle))

    phase_factor = (1.0 - Q) * Phi1 + Q * Phi_m

    return phase_factor   
       
def main():

    
##  Create Parser instructions
    aparser = argparse.ArgumentParser()
    # TIDE
    aparser.add_argument("-t", "--tidal", action="store_true", help="Tidal model - REQUIRED: TEXT FILE CONTAINING TIDAL DATA FROM COASTAL OBSERVATORY")
    aparser.add_argument("-d", "--datum", action="store", type=float, help="Depth in metres to datum below Mean Spring tide")
    aparser.add_argument("-dp", "--datum_percentage", action="store", type=float, help="Depth percentage of tidal range for Mean Spring tide")

    # LIGHT SOURCES
    aparser.add_argument("-s", "--solar", action="store_true", help="Solar Irradiance at surface and seabed")
    aparser.add_argument("-A", "--ALAN", action="store", type=float, help="2 = clear condition (2.88 PFFD), 1 = cloudy condition (6.24 PFFD), Irradiance at surface and seabed")
    aparser.add_argument("-l", "--lunar", action="store_true", help="Lunar Irradiance at surface and seabed, Alt, Az, distance") 
    # MODEL PARAMETERS
    aparser.add_argument("-T", "--time", action="store", type=float, help="Time increment in decimal hours i.e. 0.25 = 15 minutes")
    aparser.add_argument("-lat", "--latitude",type=float, action="store", help="Longitude in decimal degrees")
    aparser.add_argument("-lon", "--longitude", type=float, action="store", help="Latitude in decimal degrees")
    # MODEL OUTPUTS
    aparser.add_argument("-p", "--plots", action="store_true", help="Activate plot module")
    aparser.add_argument("-o", "--output", action="store_true", help="generate output file(s)")
    aparser.add_argument("-pr", "--prompt", action="store_true", help="Activate the prompt segments of the script for added customisation")

    aparser.add_argument("-cph", "--chlorophyll", action="store", type=float, help="Input Chlorophyll measure (UNITS)")
    aparser.add_argument("-bsc", "--backscatter", action="store", type=float, help="Input backscatter measure (UNITS)")
    aparser.add_argument("-fc", "--CDOM", action="store", type=float, help="Input fCDOM measure (UNITS)")
    aparser.add_argument("-loc", "--location", action="store", type=str, help="geographical location")
    aparser.add_argument("-thr", "--threshold", action="store", type=float, help="Threshold Intensity")
    aparser.add_argument("-stn", "--station", action="store", type=int, help="Station number - i#")
    aparser.add_argument("-TC","--TideCoef", action="store_true", help="Use TXPO tidal coefficients for location") 
    aparser.add_argument("-TPXO_d","--TPXO_dir", action="store", type=str, help="TPXO netCDF directory")

    aparser.add_argument("-start", "--start", action="store", type=str, help="Start Date (yyyy-mm-dd)")
    aparser.add_argument("-end", "--end", action="store", type=str, help="End Date (yyyy-mm-dd)")

    args = aparser.parse_args()
    sta = timeit.default_timer()  
    # if args.chlorophyll == None:
    #     print("forgotten input values remember for next time!")
    #     args.chlorophyll = 0.1
    #     args.CDOM = 0.1
    #     args.backscatter = 0.01
    if args.prompt:
        timeline = input("Do you wish to enter a date range [Y/n]:\n"
                 "Y: enter date range for model\n"
                 "n: default to last 7 days \n")
        if timeline == "Y":
            # ask the user to input date range
            start_date = input("Please enter start date: ")
            end_date = input("Please enter end date: ")
        else:
            # otherwise just run past 7 days
            start_date = datetime.datetime.utcnow() - datetime.timedelta(days=7)
            end_date = datetime.datetime.utcnow()
    else:
        if args.start:
           start_date = args.start; print("Start date: ", start_date)
        else:
           start_date = "2001-01-01"; print("Start date: ", start_date)
        if args.end:
           end_date = args.end; print("End date: ", end_date)
        else:  
           end_date = "2001-01-14"; print("End date: ", end_date)
           
    data_start_date = subprocess.getoutput(f"date --date '{start_date}' +'%F %H:%M:%S'")
    #print("Start date", data_start_date)
    data_end_date = subprocess.getoutput(f"date --date '{end_date}' +'%F %H:%M:%S'")
    #print("End date", data_end_date)

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

    # Chlorophyll = args.chlorophyll; print('Chlorophyll = ', Chlorophyll)
    # fCDOM = args.CDOM; print('fCDOM = ', fCDOM)
    # backscatter = args.backscatter; print('backscatter = ', backscatter)
    # Tide_fname = "DEV19.txt"
    # plot_index = args.plots
    
##  Compute Kd
    # kd_red, kd_green, kd_blue, r_wavenm, g_wavenm, b_wavenm, SS_wavenm = LUT_kd.kd_compute(Chlorophyll, fCDOM, backscatter)
    # kd_SS = kd_red + kd_green + kd_blue
    # print("Diffusion coefficient for Spectra \n", kd_SS)

##  Split Spectra of Sun and Moon
    # blockPrint()
    path_name = os.getcwd() + "/Required/" 
    fname_Sol = path_name + 'Solarspectra.csv'
    fname_Lun = path_name + 'Moon_spectra.csv'
    fname_ALAN = path_name + 'Lightspectra.csv'
    fname_gcirrad = path_name + 'gcirrad.dat'
    # Table 2 from Velikodsky et al., (2011) New Earth-based absolute photometry of the Moon, Icarus, 214, 30 - 45
    fname_lunar_albedo = path_name + 'lunar_albedo.dat' 
    # enablePrint()
    

    #blockPrint()
    # print("Only interested in ALAN from Falchi Map?\nY/n")
    # Falchi_inp = input()
    # ALAN_TYPE = " "
    # if Falchi_inp =="n":
    #     print("using spectral data from spreadsheets")
    print("No customisation of light source")
    if args.prompt:
        input_flag = 1
    else:
        input_flag = 0
    # blockPrint()
    SolSpec, LunSpec, ASpec, ALAN_TYPE = SpectralSplit.spectra_run(fname_Sol, fname_Lun, fname_ALAN, input_flag)
    SolBB = SolSpec['PAR'].to_numpy()[0]
    LunBB = LunSpec['PAR'].to_numpy()[0]
    SolSpec = SolSpec[['Red','Green','Blue']]
    LunSpec = LunSpec[['Red','Green','Blue']]
    enablePrint()

##  Define constants, components and arrays    
    # Isc = 521.74 # (W/m^2) A solar constant of 1373 W/m^2 (total solar irradiance) is assumed, 
    # of which approximately 38% (521.74 W/m^2) is PAR (Kirk 1994) (from Roberts paper)
    Isc = SolBB # Replaced above approximation with the integration from Solarspectra.csv 
    # k_atmos = 0.01 # Value from the Roberts paper
    # k_atmos = 0.21 # Value from Jeff Conrad 
    # Atmospheric transmission
    df_trans_atmos = pd.read_csv(fname_gcirrad, delim_whitespace=True)
    df_k_atmos = SpectralSplit.k_spectral_split('Wavelength','Trans',df_trans_atmos)
    k_atmos_bb = df_k_atmos['Broadband'].to_numpy()[0]
    df_k_atmos_spec = df_k_atmos[['Red','Green','Blue']]
    
    # Lunar spectral albedo
    df_spectral_lunar_albedo = pd.read_csv(fname_lunar_albedo, delim_whitespace=True)
    df_albedo_lunar = SpectralSplit.lunar_albedo_spectral_split('Wavelength','AverageMoonAlbedo',df_spectral_lunar_albedo)

    # interpolate the lunar spectral albedo onto the atmospheric transmission wavelengths
    gcirrad_wavs = df_trans_atmos['Wavelength'].to_numpy()
    lunar_albedo_wavs = df_spectral_lunar_albedo['Wavelength'].to_numpy()
    spec_albedo = df_spectral_lunar_albedo['AverageMoonAlbedo'].to_numpy()
    lunar_UV_albedo_fill = 0.0728 # for wavelengths shorter than 350 nm (305 - 350 nm)
    lunar_albedo_gcirrad = spectres.spectres(gcirrad_wavs, lunar_albedo_wavs, spec_albedo, fill=lunar_UV_albedo_fill, verbose=False)    
    # add column into the gcirrad dataframe
    df_trans_atmos['LunarAlbedo'] = lunar_albedo_gcirrad

    # 0.000064692 = sr value for subtending angle of the lunar disc (0.26 degree semi-diameter --> sr)
    lun_sva = 0.000064692/np.pi 

    ###################################### Select location ####################################################
    # Locations
    ##Tokyo: 35.6762, 139.6503
    ##Mumbai: 19.0760, 72.8777
    ##New York: 40.7128, -74.0060
    ##Shanghai: 31.2304, 121.4737
    ##Lagos: 6.5244, 3.3792
    ##Los Angeles: 34.0522, -118.2437
    ##Buenos Aires: -34.6037, -58.3816
    
    geo_location = str(args.location)
    if geo_location == 'Eilat':
        latitude_deg =  29.526; longitude_deg = 34.968  # lat; lon # Eilat South Tide Gauge
        Tide_fname = "TideEilatGulf_L1.csv"
    elif geo_location =='Plymouth_L4':
        latitude_deg = 50.27; longitude_deg = -4.13
        Tide_fname = "TidePlymouth_L1.csv"
    elif geo_location == 'Plymouth_Dockyard':
        latitude_deg = 50.3819; longitude_deg = -4.1927 #  lat; lon # Plymouth (L4)
        Tide_fname = "TidePlymouth_L1.csv"
    elif geo_location =='Tokyo':
        latitude_deg = 35.575; longitude_deg = 139.908
        Tide_fname = "TideTokyo_L1.csv"
    elif geo_location == 'Lagos':
        latitude_deg = 6.5244; longitude_deg = 3.3792
        Tide_fname = "TideLagos_L1.csv"
    elif geo_location == 'NewYork':
        latitude_deg = 40.7128; longitude_deg = -74.0060
        Tide_fname = "TideNewYork_L1.csv"
    elif geo_location == 'LosAngeles':
        latitude_deg = 34.0522; longitude_deg = -118.2437
        Tide_fname = "TideLosAngeles_L1.csv"
    elif geo_location == 'BuenosAires':
        latitude_deg = -34.6037; longitude_deg = -58.3816
        Tide_fname = "TideBuenosAires_L1.csv"
    elif geo_location == 'Shanghai':
        latitude_deg = 31.2304; longitude_deg = 121.4737
        Tide_fname = "TideShanghai_L1.csv" 
    elif geo_location == 'Mumbai':
        latitude_deg = 19.0760; longitude_deg = 72.8777
        Tide_fname = "TideMumbai_L1.csv"
    elif geo_location == 'Site_85':
        latitude_deg = 32.167; longitude_deg = 130.03
        Tide_fname = "TideTokyo_L1.csv" # Dummy tide file
    elif geo_location == 'Site_78':
        latitude_deg = -5.00; longitude_deg = 119.00
        Tide_fname = "TideTokyo_L1.csv" # Dummy tide file
    elif geo_location == 'Site_183':
        latitude_deg = 29.50; longitude_deg = 34.92
        Tide_fname = "TideEilatGulf_L1.csv"
    elif geo_location == 'Site_80':
        latitude_deg = -5.06; longitude_deg = 119.329
        Tide_fname = "TideTokyo_L1.csv" # Dummy tide file

    #elif geo_location == 'Ruben_1':
    #    latitude_deg = -14.279126735499915; longitude_deg = -170.70000372000004 #
    #    Tide_fname = "TideTokyo_L1.csv" # Dummy tide file

    #elif geo_location == 'GEOTAG':
    #    latitude_deg = LATLATLAT; longitude_deg = LONLONLON
    #    Tide_fname = "TideGEOTAG_L1.csv"
    
    if args.latitude:
        latitude_deg = args.latitude
        geo_location = ''
    if args.longitude:
        longitude_deg = args.longitude
    if (args.latitude and args.longitude):
        location = EarthLocation(lat=latitude_deg, lon=longitude_deg, height=1)

    if args.station:
        locations = pd.read_csv("../TamirLocDoc.csv")
        # pdb.set_trace()
        # geo_location = ''
        stationID = "i" + str(args.station)
        row_index = locations[locations["ID"]==stationID].index.values

        latitude_deg = float(locations['LatDeci'].loc[row_index])
        longitude_deg = float(locations['Lon'].loc[row_index])
        print("location (Lat/Lon)", latitude_deg, longitude_deg)
        location = EarthLocation(lat=latitude_deg, lon=longitude_deg, height=1)

    t_incr = args.time
    condition = "clear"
    ##  ALAN average broadband intensity for clear and cloudy sky conditions
    if args.ALAN==1:
        condition = 'cloudy'
    elif args.ALAN==2:
        condition = 'clear'

###################################################################################################################################
##                                 Retrieve Falchi Atlas & Kd values for Location
################################################################################################################################### 
#    
##  Read in the Falchi map of global ALAN to retrieve accurate ALAN values for the input location and Kd for R, G, B spectra. 
    if args.tidal:
        print("Accessing Falchi Map")
        print("\nUsing values taken directly from Falchi Atlas\n")
        ALAN_TYPE = "OriginalFalchi"
        sky_condition = condition

        if os.path.isfile(f"Kd_Falchi_Output_{geo_location}.csv") == False:
            df_R, df_G, df_B, ALAN_total, ALAN_mCd = Falchi_Kd_Position.read_falchi(latitude_deg, longitude_deg, sky_condition, input_flag)
            Kd_Falchi_Output = pd.concat([df_R["Lat"], df_R["Lon"], df_R["Month"], df_R["ALAN_R (uW/m^2)"], df_R["Kd_R"], df_G["ALAN_G (uW/m^2)"], df_G["Kd_G"], df_B["ALAN_B (uW/m^2)"], df_B["Kd_B"], ALAN_total, ALAN_mCd], keys = ["Latitude", "Longitude", "Month", "ALAN_R_(uW/m^2)", "Kd_R", "ALAN_G_(uW/m^2)", "Kd_G", "ALAN_B_(uW/m^2)", "Kd_B", "Irr_(uW/m^2)", "Lum_(uCd/m^2)"], axis=1)
            pd.DataFrame(Kd_Falchi_Output).to_csv(f"Kd_Falchi_Output_{geo_location}.csv") 
            # Label values for use throughout script
            ALAN_R = Kd_Falchi_Output["ALAN_R_(uW/m^2)"]; ALAN_G = Kd_Falchi_Output["ALAN_G_(uW/m^2)"]; ALAN_B = Kd_Falchi_Output["ALAN_B_(uW/m^2)"]; Kd_R = Kd_Falchi_Output["Kd_R"]; Kd_G = Kd_Falchi_Output["Kd_G"]; Kd_B = Kd_Falchi_Output["Kd_B"]; ALAN_total = Kd_Falchi_Output["Irr_(uW/m^2)"][0]; ALAN_mCd = Kd_Falchi_Output["Lum_(uCd/m^2)"][0]; Month_Kd = Kd_Falchi_Output["Month"] 
        else:
            Kd_Falchi_Output = pd.read_csv(f"Kd_Falchi_Output_{geo_location}.csv")
            # Label values for use throughout script
            ALAN_R = Kd_Falchi_Output["ALAN_R_(uW/m^2)"]; ALAN_G = Kd_Falchi_Output["ALAN_G_(uW/m^2)"]; ALAN_B = Kd_Falchi_Output["ALAN_B_(uW/m^2)"]; Kd_R = Kd_Falchi_Output["Kd_R"]; Kd_G = Kd_Falchi_Output["Kd_G"]; Kd_B = Kd_Falchi_Output["Kd_B"]; ALAN_total = Kd_Falchi_Output["Irr_(uW/m^2)"][0]; ALAN_mCd = Kd_Falchi_Output["Lum_(uCd/m^2)"][0]; Month_Kd = Kd_Falchi_Output["Month"] 
            latitude_deg = Kd_Falchi_Output['Latitude'][0]
            longitude_deg = Kd_Falchi_Output['Longitude'][0]
            location = EarthLocation(lat=latitude_deg, lon=longitude_deg, height=1)
        
        R = 0; G = 0; B=0
        try:
            data = {"Red": [ALAN_R[0]], "Green": [ALAN_G[0]], "Blue": [ALAN_B[0]]} # 0 = month counter from Falchi output, but Falchi Atlas does not very with time. Hence 0
        except:
            data = {"Red": [R], "Green": [G], "Blue": [B]}
        ASpec = pd.DataFrame(data,columns=["Red", "Green", "Blue"])

    # Assign empty variables
    sol = []; night = []; dec_day = []; Io = []; IBT= []; solar_surface = []; solar_datum = []; solar_alt = [];
    alt = []; az = []; phase = []; I = []; lIBT = []; I_atmos = []; lunar_surface= []; lunar_datum = []
    A = []; aIBT = []; ALAN_surface = []; ALAN_datum = []
    TL = []; t = []; waterdepth = []; tide_h = []; ftide_h = []; Zc = []
    # kd as a function of time
    kd_red_t = []; kd_green_t = []; kd_blue_t = []; kd_BB_t = [];
    frames = []; Lun_frames = []; Sol_frames = []; A_frames = []
    date_record = []; day_summed =[]; fullmoon_mask = []
    datum_intensity_Lunar = []; datum_intensity_ALAN = []; datum_intensity_Solar = []


    # col_names_r = ['620', '630', '640', '650', '660', '670', '680', '690', '700', '710', '720', '730', '740']
    # col_names_g = ['500', '510', '520', '530', '540', '550', '560']
    # col_names_b = ['400', '410', '420', '430', '440', '450', '460', '470', '480', '490', '500']
    col_names_SS = ['Red', 'Green', 'Blue'] # To this end SS Solar Spec and LS Lunar Spec are equal in index
    col_names_SS_sol = ['I_Red(W/m2)','I_Green(W/m2)','I_Blue(W/m2)']
    col_names_SS_lun = ['I_Red(uW/m2)','I_Green(uW/m2)','I_Blue(uW/m2)']
    col_names_SS_A = ['I_Red(uW/m2)','I_Green(uW/m2)','I_Blue(uW/m2)']
    
    col_names_SS_sol_datum = ['I_Red_datum(W/m2)','I_Green_datum(W/m2)','I_Blue_datum(W/m2)']
    col_names_SS_lun_datum = ['I_Red_datum(uW/m2)','I_Green_datum(uW/m2)','I_Blue_datum(uW/m2)']
    col_names_SS_A_datum = ['I_Red_datum(uW/m2)','I_Green_datum(uW/m2)','I_Blue_datum(uW/m2)']
    
    # Spectral Intensity Dataframes
    SolI_SSb = pd.DataFrame(columns=col_names_SS_sol_datum); LunI_LSb = pd.DataFrame(columns=col_names_SS_lun_datum); AI_ASb = pd.DataFrame(columns=col_names_SS_A_datum)
    SolI_SS = pd.DataFrame(columns=col_names_SS_sol); LunI_LS = pd.DataFrame(columns=col_names_SS_lun); AI_AS = pd.DataFrame(columns=col_names_SS_A)
    SolI_SS_TOA = pd.DataFrame(columns=col_names_SS_sol); LunI_LS_atm = pd.DataFrame(columns=col_names_SS_lun)
    SolI_SSRes = pd.DataFrame(columns=col_names_SS); LunI_LSRes = pd.DataFrame(columns=col_names_SS); AI_ASRes = pd.DataFrame(columns=col_names_SS)
    crit_depth = np.arange(0, 200, 0.25).tolist()

###################################################################################################################################
##                                 Retrieve Tidal parameters for Location 
################################################################################################################################### 
# Using tidal data from a corrected scientific tide measurement station run the tidal model UTide to solve for harmonic coefficients
    
    if (args.tidal):
        
        print('Running tidal model...')

        if args.datum_percentage:
            datum_percentage = args.datum_percentage
        else:
            datum_percentage = 0.1
        
        if args.TPXO_dir:
           t_incr_m = np.int64(t_incr*60)
           step = np.timedelta64(t_incr_m, 'm') # set stepsize
           # set time range
           times = np.arange(np.datetime64(start_date), np.datetime64(end_date) + step , step)
           TL = tide_predict(args.TPXO_dir, latitude_deg,longitude_deg, times)
           TL_TPXO = TL
        else:
           ## Read in tidal data 
           ## Tidal data from a variety of sources 
           ## 1. PSMSL
           ## https://www.psmsl.org/data/
           ## 2. University of Hawaii Sea Level Centre
           ## http://uhslc.soest.hawaii.edu/data/
           ## 3. Generate data from the TPXO data files 
           ## https://tpxows.azurewebsites.net/  
           ## ==> generate weekly tidal data, multiple times and then create a ?_L1.csv file from about a month's worth of data
           ## Use routines found in the Tide_PrePro directory and create a suitable ?_L1.csv file
           ## Put data in Model/Required/TideGauge directory

           ##https://www.tide-forecast.com/ - good site for eyeballing how accurately the tide is recreated.

           print('     using tidegauge data for tidal coefficients')
           tidepath = os.getcwd() + "/Required/TideGauge/" + Tide_fname # path of data file output
           tides_ = pd.read_csv(tidepath, delimiter=',', engine='python') #usecols=np.arange(16,48), engine='python')
           df = pd.DataFrame(tides_)
           for i in range(len(tides_)):
           ## Convert/alter DataFrame Variables to type datetime and float
              tide_date = df.iloc[i,0]
              tide_level =df.iloc[i,1]
              tide_level = float(tide_level)

           ## Convert date to datetime format
              T = datetime.datetime.strptime(tide_date, '%Y/%m/%d %H:%M:%S') # All stations to a standard timestamp
           ## append variables to lists
              TL.append(tide_level)
              t.append(T)

           ## Convert time back to DataFrame
           t_df = pd.to_datetime(t)
           ## Convert time to format of UTide ##'date2num' function only seems to work with pandas DF ##
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
           dts = (pd.DataFrame(columns=['NULL'],index=pd.date_range(start_date, end_date,freq='60T')).between_time('00:00','23:59').index.strftime('%Y-%m-%d %H:%M:%S').tolist())
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
        
    #                                         TidalLight - MODEL BREAKDOWN
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                                           SOLAR MODULE: args.solar
    #
    # Solar module: Contains intensity at seasurface and tidally modulated results (args.seabed)
    # Time loop is created for a given year (input at command line -y) in increments specified by -T 
    # within loop pysolar modules utilise ephemeris calculations for azimuth, altitude and airmass to calculate Intensity at the top of the atmosphere
    # Intensity at sea level follows the method described in Roberts et al., (2018) Limnol. Oceangr. 63, 2018, 91-106
    # Tidal model is produced for the given time inputs.
    # Intensity attenuation through the water column is modulated by the tidal model (args.seabed) 
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                                           LUNAR MODULE: args.lunar
    #
    # Contains forward modelling functions from astropy and astroplan packages
    # predict azimuth, altitude, and distance from ICRS reference frame (International Celestial Reference Frame)
    # Lunar phase calculated as decimal where 1 = full-moon
    # 
    # NOTE: Lunar albedo is wavelength dependent the albedo for specific wavelengths will be important for accuracy and validity of future spectral decomp of intensity
    # Plan is to use pysolar airmass function with the output of astropy altitude for consistency

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                                           TIDAL MODULE: args.seabed
    #
    # Models tidal range for a given location using the UTide function which creates the relevant harmonic coefficients
    # for a data set from a coastal observatory with time-series data. The data for the UK is collected from the BODC website.
    # Data file must be in the same location as the script*
    # Roberts et al., (2018) describes the proccess through which intensity is attenuated through the water columnn.

    ##############################################
    #      START TidalLight MODEL LOOP 
    ##############################################
    month_counter = 0 
    TPXO_count = 0    
    # for delt in range(date_start, date_end):
    for model_date in daterange(date_start, date_end):
        year = model_date.year # Year in datetime format
        day = model_date.day # Day in datetime format
        month = model_date.month # Month in datetime format
        # Print first date of each month as a progress counter
        if month != month_counter: 
            month_counter = month
            print(model_date)
        else:
            pass
        if args.tidal:
            monthly_Kd_idx = Kd_Falchi_Output[Kd_Falchi_Output['Month'] == month].index[0] # Match the month of the loop to the correct Kd 
            kd_red = Kd_R[monthly_Kd_idx]; kd_green = Kd_G[monthly_Kd_idx]; kd_blue = Kd_B[monthly_Kd_idx] # Assign monthly Kd to variable
            KD = [kd_red, kd_green, kd_blue] # create a Kd array
            #KD_Bb = scipy.integrate.simps(KD, dx=1) # Calculate Broadband Kd
            KD_Bb = np.mean(KD) # Calculate Broadband Kd
            
        enablePrint()
        print(model_date)
        blockPrint()
        for hour in range(0,24):
            
            minutes = np.arange(start=0, stop=1, step=t_incr) # create minute intervals
            for Minute in minutes:
                H = ((hour+Minute)/24) # create Decimal component of the day 
                minute = int(Minute*60) # reduce to a 60-base value 
                if minute>=60:
                    minute = 0 
                date = datetime.datetime(year, month, day, hour, minute, tzinfo=datetime.timezone.utc) # Create date - accurate to minute incr.
                
                date_record.append(date)
                jday = date.timetuple().tm_yday # Retrieve date as a tuple
                dday = int(jday)+H # Combine times as a decimal day
                Ecc = (1.0+1.67E-2*np.cos(2.0*np.pi*float(dday-3)/365.0))**2 # correction for eccentricity of Earth's orbit
                dec_day.append(dday) 
                altitude_deg = get_altitude(latitude_deg, longitude_deg, date) # Retrieve altitude of the sun
                solar_alt.append(altitude_deg)

                # Print Sunset/Sunrise times to terminal
                # if 0.1>altitude_deg>-0.1:
                #     if hour<12: 
                #         print("Sunrise: ", hour,":",minute,"\nSolar Alt. (degrees)", altitude_deg) # disable blockPrint for output
                #     elif hour>12:
                #         print("Sunset: ", hour,":",minute,"\nSolar Alt. (degrees)", altitude_deg)
                #
                # In the United Kingdom, there is a legally enforced lighting-up time, 
                # defined as from half an hour after sunset to half an hour before sunrise
                # https://en.wikipedia.org/wiki/Lighting-up_time
                # Civil twilight when the Sun's centre is 6 degrees below the horizon, 
                # is roughly equivalent to Lighting-up Time.
                # https://www.rmg.co.uk/stories/topics/dawn-dusk-twilight

                # Playing Night time.. Day time: https://www.youtube.com/watch?v=Ln2Xq8fCNI8
                #enablePrint()
                if -18<=altitude_deg<=0:
                    twilight_df = twilight_spitschan.main(altitude_deg, "rural", "skye")
                    if -6<=altitude_deg<=0:
                       solar_day = float(0.75)
                    else:
                       solar_day = float(0.5)
                elif altitude_deg < -18:
                    solar_day = float(0)
                    Night = 2700
                    night.append(Night) # Night values are assigned arbitrarily for plotting effect
                else:
                    solar_day = float(1)
                    Night = -600
                    night.append(Night) # Night values are assigned arbitrarily for plotting effect
                #blockPrint()
                
                sol.append(solar_day)
######################################################## 
                # IRRADIANCE CALCULATIONS
########################################################
                df_TOA = gcirrad.TOA_irradiance(df_trans_atmos, Ecc)
                # Top of Atmosphere values of the SolSpec
                SolSpecTOA = SpectralSplit.spectral_split('Wavelength','Fo',df_TOA)
                PARTOA = SolSpecTOA[['PAR']].to_numpy()[0]
                SolSpecTOA = SolSpecTOA[['Red','Green','Blue']]

                # SOLAR CALCULATION 
                if (args.solar):
##                   Calculate variable using Pysolar
                    #airmass = get_air_mass_ratio(altitude_deg)  ###### NOTE ####### diffusivity through airmass is wavelength specific - currently using broadband attenuation for each wvlnth
                    if (altitude_deg > 0.):
                       airmass = get_air_mass_kasten_young(altitude_deg)
                       # Run a more sophisticated atmospheric model
                       df_dir_dif = gcirrad.dir_dif_irradiance(df_trans_atmos, altitude_deg, Ecc)
                       # Surface values of the SolSpec
                       SolSpecSurface = SpectralSplit.spectral_split('Wavelength','Ed',df_dir_dif)
                       PARSurface = SolSpecSurface[['PAR']].to_numpy()[0]
                       SolSpecSurface = SolSpecSurface[['Red','Green','Blue']]
#                    else:
#                       airmass = 30.0
##                  Calculate Intensity based on the Earth's orbital eccentricity
                    #Iatmos = Isc*Ecc
##                  Calculate Iatmos for each wavelength
                    SSurface = np.array([])
                    #SAtmos = np.array([])
                    STOA = np.array([])
                    
                    for ss in range(len(SolSpecTOA.columns)):
                        #Iatmos_SS = SolSpec.iloc[:,ss]*Ecc # Irradiance at atmosphere (Masters, 2004): Renewable and Efficient Power Systems Ch7
                        
                        ITOA_SS = SolSpecTOA.iloc[:,ss]
                        STOA = np.append(STOA, ITOA_SS)
                        #k_atmos_spec = df_k_atmos_spec.iloc[:,ss].to_numpy()[0]
                        if solar_day==0:
                            IO_SS = float(0)
                        elif (solar_day >= 0.5 and solar_day < 1.0):
                            # Bugfix: twilight_spitzchan data are in log base 10
                            #IO_SS = np.exp(twilight_df["Twilight_broadband(log(W/m2))"][ss]) #Twilight component
                            IO_SS = 10**(twilight_df["Twilight_broadband(log(W/m2))"][ss]) #Twilight component
                            if IO_SS > 15:
                               IO_SS=0
                        else:
                            IO_SS = SolSpecSurface.iloc[:,ss]
                            #IO_SS = Iatmos_SS*np.sin(np.deg2rad(altitude_deg))*np.exp(-1.*airmass*k_atmos_spec) # Irradiance at surface
                        SSurface = np.append(SSurface, IO_SS)

                    SolI_SS_TOA = SolI_SS_TOA.append(pd.Series(STOA, index=SolI_SS_TOA.columns),ignore_index=True)
                    SolI_SS = SolI_SS.append(pd.Series(SSurface, index=SolI_SS.columns),ignore_index=True)  
                    ## Need to introduce pd.concat in future, but this particular coding currently not producing expected result
                    #SolI_SS_TOA = pd.concat((pd.DataFrame(STOA, index=SolI_SS_TOA.columns).T, SolI_SS_TOA))
                    #SolI_SS = pd.concat((pd.DataFrame(SSurface, index=SolI_SS.columns).T, SolI_SS))  
                    
                    if solar_day==0:
                        IO = float(0)
                    elif (solar_day >= 0.5 and solar_day < 1.0):
                        # Bugfix: twilight_spitzchan data are in log base 10
                        #IO = np.exp(twilight_df["Twilight_PAR(log(W/m2))"][0])
                        IO = 10**(twilight_df["Twilight_PAR(log(W/m2))"][0])
                        if IO > 15:
                           IO = float(0)
                    else:
                        IO = PARSurface # Broadband
                        #IO = Iatmos*np.sin(np.deg2rad(altitude_deg))*np.exp(-1.*airmass*k_atmos_bb) # Broadband
                    Io.append(float(IO))
                
######################################################## 

                # ALAN CALCULATION  
                ## Calculate Iatmos for each wavelength
                ASurface = np.array([])
                if args.ALAN:
                    for aa in range(len(ASpec.columns)):
                        #if solar_day==1:
                        # Only "switch on" ALAN after "Lighting up time" - which is end of Civic Twilight
                        if solar_day>=0.75:
                            Aint = 0
                        else:
                            Aint = ASpec.iloc[:,aa]
                        ASurface = np.append(ASurface, Aint)
                    
                    AI_AS = AI_AS.append(pd.Series(ASurface, index=AI_AS.columns),ignore_index=True)  
                    ## Need to introduce pd.concat in future, but this particular coding currently not producing expected result
                    #AI_AS = pd.concat((pd.DataFrame(ASurface, index=AI_AS.columns).T, AI_AS)) 
                    #if solar_day==1:
                    # Only "switch on" ALAN after "Lighting up time" - which is end of Civic Twilight
                    if solar_day>=0.75:
                        Aint = 0
                    else:
                        Aint = ALAN_total
                    A.append(Aint) # ALAN array at surface

######################################################## 
                # LUNAR CALCULATIONS
                if (args.lunar):
                    # Lunar albedo: 
                    # Lane & Irvine (1973) ~12% for full moon DOI: 10.1086/111414
                    # Buratti et al., (1996) ~16% (zero phase value)
                    # Lunar_albedo = 0.16
                    # Use spectral albedo from Table 2 from Velikodsky et al., (2011)
                    # For broadband using df_albedo_lunar['Broadband'] - which is approx 12%
                    dt = Time(date) # Time within model
                    EarthMoon = get_moon(dt, location) # Relative positioning              
                    moon_icrs = EarthMoon.transform_to('icrs') # Relative Ephemeris calculation    
                    moonaltaz = moon_icrs.transform_to(AltAz(obstime=dt, location=location)) # Transform positioning to Alt, Az
                    alt_ = float(Angle(moonaltaz.alt).degree) # Lunar Altitude 
                    az_ = float(Angle(moonaltaz.az).degree) # Lunar Azimuth                     
                    Phase = float(moon.moon_illumination(dt)) # Lunar Phase 
                    Phase_lb = lumme_bowell(moon.moon_phase_angle(dt).value) # Lunar Phase function - Lumme-Bowell
                    albedo_phased = df_albedo_lunar['Broadband'].to_numpy()[0]*Phase_lb
                    # Iatmos in place of Isc - to correct for eccentricity of Earth (moon) orbit 
                    # Broadband top of atmosphere lunar 
                    Lunar_refl = PARTOA*lun_sva*albedo_phased*1000000. # convert from W/m^2 to uW/m^2  
                    airmass = get_air_mass_kasten_young(alt_)
                    if alt_<0:
                        Iol = 0
                    elif altitude_deg>0:
                        Iol = 0
                    else:
                        # Use the Gregg and Carder atmospheric model, with correction for the lunar spectral albedo (lunar=True)
                        # and the lunar phase (multiplied by the lunar solid view angle)
                        df_dir_dif = gcirrad.dir_dif_irradiance(df_trans_atmos, alt_, Ecc, lunar=True, phase_sva=Phase_lb*lun_sva)
                        LunSpecSurface = SpectralSplit.spectral_split('Wavelength','Ed',df_dir_dif)
                        LunPARSurface = LunSpecSurface[['PAR']].to_numpy()[0]
                        LunSpecSurface = LunSpecSurface[['Red','Green','Blue']]
                        Iol = LunPARSurface*1000000. # convert from W/m^2 to uW/m^2
                        #Iol = (Lunar_refl*np.sin(np.deg2rad(alt_))*np.exp(-1.*airmass*k_atmos_bb)) # Lunar Irr at surface
                    alt.append(alt_)
                    phase.append(Phase)
                    az.append(az_)
                    I.append(float(Iol))
                    I_atmos.append(float(Lunar_refl))

                    LSurface = np.array([])
                    LAtmos = np.array([])
                    #for ll in range(len(LunSpec.columns)):
                    #pdb.set_trace()
                    for ll in range(len(SolSpec.columns)):
                    #for ll in range(len(LunSpecSurface.columns)):
                        #Lunar_refl_LS = (LunSpec.iloc[:,ll]) # spectral data is a measurement which has been extrapolated to top of atmosphere. 
                        # CHANGE TJS: Need to add in the Phase information here
                        #Lunar_refl_LS = (LunSpec.iloc[:,ll]) * Phase_lb # spectral data is a measurement which has been extrapolated to top of atmosphere. 
                        #Lunar_refl_LS = SolSpecTOA.iloc[:,ll]*Ecc*lun_sva*albedo_phased # spectral data is a measurement which has been extrapolated to top of atmosphere. 
                        Lunar_refl_LS = SolSpecTOA.iloc[:,ll]*lun_sva*albedo_phased # spectral data is a measurement which has been extrapolated to top of atmosphere. 
                        Lunar_refl_LS = Lunar_refl_LS*1000000 # convert from W/m^2 to uW/m^2
                        LAtmos = np.append(LAtmos, Lunar_refl_LS)
                        #k_atmos_spec = df_k_atmos_spec.iloc[:,ss].to_numpy()[0]                      
                        if alt_<0:
                            Iol_LS = 0
                        # Lunar irradiance only meaningful after sunset
                        elif altitude_deg>0:
                            Iol_LS = 0
                        else:
                            Iol_LS = LunSpecSurface.iloc[:,ll]*1000000. # convert from W/m^2 to uW/m^2
                            #Iol_LS = (Lunar_refl_LS*np.sin(np.deg2rad(alt_))*np.exp(-1.*airmass*k_atmos_spec))
                            #Iol_LS = Iol_LS*1000000 # convert from W/m^2 to uW/m^2 
                        LSurface = np.append(LSurface, Iol_LS)
                    
                    LunI_LS = LunI_LS.append(pd.Series(LSurface, index=LunI_LS.columns),ignore_index=True)
                    LunI_LS_atm = LunI_LS_atm.append(pd.Series(LAtmos, index=LunI_LS_atm.columns),ignore_index=True)
                    ## Need to introduce pd.concat in future, but this particular coding currently not producing expected result
                    #LunI_LS = pd.concat((pd.DataFrame(LSurface, index=LunI_LS.columns).T, LunI_LS))
                    #LunI_LS_atm = pd.concat((pd.DataFrame(LAtmos, index=LunI_LS_atm.columns).T, LunI_LS_atm))
                    
                    #if (alt_ > 5.0 and altitude_deg < 0):
                    #   enablePrint()
                    #   pdb.set_trace()
                    #   blockPrint()
                    
                    if Phase >=0.9:
                        fullmoon_mask.append(1)
                    else: 
                        fullmoon_mask.append(0)
######################################################## 
                # TIDAL CALCULATIONS  
########################################################
                if args.tidal:
                    # The tide at the precise time points already calculated earlier in the model
                    # for the TPXO data
                    if args.TPXO_dir:
                       tide_h.append(TL_TPXO[TPXO_count])
                       depth_to_datum = float(TL_TPXO[TPXO_count] - datum) # depth_to_datum = height of water column above datum
                       if TL_TPXO[TPXO_count]<=datum:
                          depth_to_datum = float(0)
                       ftide_h.append(TL_TPXO[TPXO_count])
                       waterdepth.append(depth_to_datum) 
                       TPXO_count += 1
                    # Where tide gauge data has been used to back calculate the tidal harmonics
                    # calculate the tide for the precise time points here
                    else:
                       date2 = datetime.datetime(year, month, day, hour, minute) # create new date variable
                       date2 = str(date2) # convert to string
                       fmt_date = datetime.datetime.strptime(date2, '%Y-%m-%d %H:%M:%S') # for some reason UTide requires it be converted back to a datetime before running
                       tidetime = pd.to_datetime(fmt_date) # convert to dataframe
                       tide_time = mdates.date2num(tidetime.to_pydatetime()) # mdates is the required format for UTide, this is the only way I managed to get the function to work
                       reconst = utide.reconstruct(tide_time, c) # Reconstruct a tidal model for the selected date range using the constants 'c' determined by the tidal module
                       tide_h.append(reconst.h) # retrieve tide height (sea level) and store in a list
                       depth_to_datum = float(reconst.h - datum) # depth_to_datum = height of water column above datum
                       if reconst.h<=datum:     
                          depth_to_datum = float(0) 
                       ftide_h.append(float(reconst.h))
                       waterdepth.append(depth_to_datum)
                    
                    #kPAR = ((0.5+0.5*np.cos(2*np.pi*dday)/365))*0.1*reconst.h+0.4 # Old computation of diffusivity in atmosphere for PAR (Masters, 2004)
                    # need to update / append the value of Kd(s) for the output file
                    kd_red_t.append(kd_red); kd_green_t.append(kd_green); kd_blue_t.append(kd_blue); kd_BB_t.append(KD_Bb)

########################################################
        # ATTENUATION OF INTENSITY TO SEABED                 
######################################################## 
                    # SOLAR
                    if args.solar:
                        # Broadband 
                        if depth_to_datum==0: # when datum level is above the wter line Intensity at seabed (datum) is not modulated by the water
                            iBT = float(IO)
                        else:
                            iBT = float(IO*np.exp(-KD_Bb*depth_to_datum)) # Roberts et al., 2018 ## Substitute ion KD_Bb for crit depth
                            
                        IBT.append(iBT) # IBT: Intensity Below Tide
                        # Spectra
                        SS_array = np.array([]) # Solar Surface
                        SSRes_array = np.array([]) # Solar Residuals: Surface - Seabed
                        for ss in range(len(SolSpec.columns)):
                            iBT_SS = float(SSurface[ss]*np.exp(-KD[ss]*depth_to_datum)) # Beers law - Io.exp(-Kd*z)
                            
                            if depth_to_datum==0:
                                iBT_SS = float(SSurface[ss])
                            SS_array = np.append(SS_array, iBT_SS)
                            SSRes_array = np.append(SSRes_array, (SSurface[ss]-iBT_SS))
                        SolI_SSb = SolI_SSb.append(pd.Series(SS_array, index=SolI_SSb.columns), ignore_index=True)
                        SolI_SSRes = SolI_SSRes.append(pd.Series(SSRes_array, index=SolI_SSRes.columns), ignore_index=True)
                        ## Need to introduce pd.concat in future, but this particular coding currently not producing expected result
                        #SolI_SSb = pd.concat((pd.DataFrame(SS_array, index=SolI_SSb.columns).T, SolI_SSb))
                        #SolI_SSRes = pd.concat((pd.DataFrame(SSRes_array, index=SolI_SSRes.columns).T, SolI_SSRes))

########################################################  
                          
                    # LUNAR 
                    if args.lunar:
                        # Broadband
                        if depth_to_datum==0:
                            LIBT = float(Iol)
                        else:
                            LIBT = float(Iol*np.exp(-KD_Bb*depth_to_datum))
                        lIBT.append(LIBT)
                        # Spectra
                        LSb_array = np.array([])
                        LRes_array = np.array([])
                        for ll in range(len(LunSpec.columns)):
                            LIBT_LS = float(LSurface[ll]*np.exp(-KD[ll]*depth_to_datum))
                            if depth_to_datum==0:
                                LIBT_LS = float(LSurface[ll])
                            LSb_array = np.append(LSb_array, LIBT_LS)
                            LRes_array = np.append(LRes_array, (LSurface[ll]-LIBT_LS))
                        LunI_LSb = LunI_LSb.append(pd.Series(LSb_array, index=LunI_LSb.columns), ignore_index=True)
                        LunI_LSRes = LunI_LSRes.append(pd.Series(LRes_array, index=LunI_LSRes.columns), ignore_index=True)
                        ## Need to introduce pd.concat in future, but this particular coding currently not producing expected result
                        #LunI_LSb = pd.concat((pd.DataFrame(LSb_array, index=LunI_LSb.columns).T, LunI_LSb))
                        #LunI_LSRes = pd.concat((pd.DataFrame(LRes_array, index=LunI_LSRes.columns).T, LunI_LSRes))
                        

########################################################               

                    # ALAN
                    if args.ALAN:
                        # Broadband
                        if depth_to_datum==0:
                            AIBT = float(Aint)
                        else:
                            AIBT = float(Aint*np.exp(-KD_Bb*depth_to_datum))
                        aIBT.append(AIBT)
                        # Spectra
                        ASb_array = np.array([])
                        ARes_array = np.array([])
                        for aa in range(len(ASpec.columns)):
                            AIBT_AS = float(ASurface[aa]*np.exp(-KD[aa]*depth_to_datum))
                            if depth_to_datum==0:
                                AIBT_AS = float(ASurface[aa])
                            
                            ASb_array = np.append(ASb_array, AIBT_AS)
                            ARes_array = np.append(ARes_array, (ASurface[aa]-AIBT_AS))
                        AI_ASb = AI_ASb.append(pd.Series(ASb_array, index=AI_ASb.columns), ignore_index=True)
                        AI_ASRes = AI_ASRes.append(pd.Series(ARes_array, index=AI_ASRes.columns), ignore_index=True)
                        ## Need to introduce pd.concat in future, but this particular coding currently not producing expected result
                        #AI_ASb = pd.concat((pd.DataFrame(ASb_array, index=AI_ASb.columns).T, AI_ASb))
                        #AI_ASRes = pd.concat((pd.DataFrame(ARes_array, index=AI_ASRes.columns).T, AI_ASRes))
########################################################
#             END TidalLight MODEL LOOP               
######################################################## 


    enablePrint()
##################################################################################################
#                                Broadband Crit. Depth     
################################################################################################## 
    Bb_Zc = [];R_Zc = [];G_Zc = [];B_Zc = []
    for i in range(len(KD)):
        if args.ALAN:
            datum_intensity_ALAN.append(sum(AI_ASb.iloc[:,i])) # Dosage ALAN
        if args.lunar:
            datum_intensity_Lunar.append(sum(LunI_LSb.iloc[:,i])) # Dosage Lunar
        if args.solar:
            datum_intensity_Solar.append(sum(SolI_SSb.iloc[:,i])) # Dosage Solar
            
        A_zBb = []
        zBb = []
        critical_depth_Bb = []
        KD = [Kd_R[i], Kd_G[i], Kd_B[i]] 
        #KD_Bb = scipy.integrate.simps(KD, dx=1) # Calculate Broadband Kd
        KD_Bb = np.mean(KD) # Calculate Broadband Kd
        if args.ALAN:
            if args.threshold:
                threshold_intensity = args.threshold
            else:
                threshold_intensity = 0.102 # Batnes et al., 2015 uW/m^2 # my attempt at converting uW/cm2 to uW/m2 = 0.03... # 4.1x10^-5 uW/cm^2/nm (from Tamir et al., 2017)  
            
            for zz in range(len(crit_depth)):
                AIBT_z = float(ALAN_total*np.exp(-KD_Bb*crit_depth[zz]))
                A_zBb.append(AIBT_z)
                zBb.append(crit_depth[zz])
                if AIBT_z < threshold_intensity: # this value is in uW/m^2, 
                    critical_depth_Bb.append(crit_depth[zz])
                    break
                else:
                    continue
    
##################################################################################################
#                               Spectral Crit. Depth 
##################################################################################################
            A_zR = []; A_zG = []; A_zB = []
            zR = []; zG = []; zB = []
            critical_depth_R = []; critical_depth_G = []; critical_depth_B = []
            for zz in range(len(crit_depth)):
                AIBT_zR = float(ALAN_R[i]*np.exp(-Kd_R[i]*crit_depth[zz]))
                A_zR.append(AIBT_zR) # append the value of Irradiance at the corresponding depth "zz"
                zR.append(crit_depth[zz]) # record the value of depth. 
                if AIBT_zR < 150: # uW/m^2, value taken from Davies et al., 2020 (Biologically Relevant light)
                    critical_depth_R.append(crit_depth[zz])
                else:
                    
                    continue
            
            for zz in range(len(crit_depth)):
                AIBT_zG = float(ALAN_G[i]*np.exp(-Kd_G[i]*crit_depth[zz]))
                A_zG.append(AIBT_zG) # append the value of Irradiance at the corresponding depth "zz"
                zG.append(crit_depth[zz]) # record the value of depth. 
                if AIBT_zG < 0.75: # uW/m^2, value taken from Davies et al., 2020 (Biologically Relevant light)

                    critical_depth_G.append(crit_depth[zz])
                else:

                    continue
            
            for zz in range(len(crit_depth)):
                AIBT_zB = float(ALAN_B[i]*np.exp(-Kd_B[i]*crit_depth[zz]))
                A_zB.append(AIBT_zB) # append the value of Irradiance at the corresponding depth "zz"
                zB.append(crit_depth[zz]) # record the value of depth. 
                if AIBT_zB < 0.19: # uW/m^2, value taken from Davies et al., 2020 (Biologically Relevant light)

                    critical_depth_B.append(crit_depth[zz])
                else:
                    continue
            # Clunky Safeguarding - if thresholds are too high.
            if critical_depth_G == []:
                critical_depth_G = [0]
            if critical_depth_B == []:
                critical_depth_B = [0]
            if critical_depth_R == []:
                critical_depth_R = [0]
            Bb_Zc.append(critical_depth_Bb[0])
            R_Zc.append(critical_depth_R[0]); G_Zc.append(critical_depth_G[0]); B_Zc.append(critical_depth_B[0])
    print(f"Dosage of Solar Irr. at datum\n{date_start} - {date_end}\nRGB:  {datum_intensity_Solar} J/m^2")
    print(f"Dosage of Lunar Irr. at datum\n{date_start} - {date_end}\nRGB:  {datum_intensity_Lunar} \u03bcJ/m^2")
    print(f"Dosage of ALAN at datum\n{date_start} - {date_end}\nRGB:  {datum_intensity_ALAN} \u03bcJ/m^2")

    Sol_datum_dosage_df = pd.DataFrame(np.array(datum_intensity_Solar).reshape(-1,3), columns = list("RGB"))
    Lun_datum_dosage_df = pd.DataFrame(np.array(datum_intensity_Lunar).reshape(-1,3), columns = list("RGB"))
    ALAN_datum_dosage_df = pd.DataFrame(np.array(datum_intensity_ALAN).reshape(-1,3), columns = list("RGB"))


    #### FEEL LIKE THIS COULD BE CONDENSED INTO A FUNCTION..
################################################################################
#                               Sums and checks 
################################################################################
    if args.solar:
        # Solar Total: Spectral check
        Total_Bb = max(Io)
        Tot = np.array([])
        Tot_atm = np.array([])
        for ii in range(len(SolSpec.columns)):
            MAX = max(SolI_SS.iloc[:,ii]) # At the surface
            MAX_atm = max(SolI_SS_TOA.iloc[:,ii])
            Tot = np.append(Tot, MAX)
            Tot_atm = np.append(Tot_atm, MAX_atm)
        Total_ss = np.sum(Tot)
        Total_ss_atm = np.sum(Tot_atm)
        
        print("SOLAR:\n     Max PAR[W/m2] at surface =", Total_Bb)
        print("     Max PAR[W/m2] at TOA =", Isc)
        print("     Max summed RGB[W/m2] at surface =", Total_ss)
        print("     Max summed RGB[W/m2] at TOA =", Total_ss_atm)
        
        #enablePrint()
        #pdb.set_trace()
        #blockPrint()
        
    if args.lunar:
        # Lunar Total: Spectral check 
        Total_Bb = max(I)
        Total_Bb_atm = max(I_atmos)
        Tot = np.array([])
        Tot_atm = np.array([])
        for ii in range(len(LunSpec.columns)):
            MAX = max(LunI_LS.iloc[:,ii]) # At the surface
            MAX_atm = max(LunI_LS_atm.iloc[:,ii])
            Tot = np.append(Tot, MAX)
            Tot_atm = np.append(Tot_atm, MAX_atm)
        Total_ls = np.sum(Tot)
        Total_ls_atm = np.sum(Tot_atm)

        print("LUNAR: \n     Max PAR[uW/m2] at surface =", Total_Bb)
        print('     Max PAR[uW/m2] at TOA =', Total_Bb_atm)
        print("     Max summed RGB[uW/m2] at surface =", Total_ls)
        print("     Max summed RGB[uW/m2] at TOA =", Total_ls_atm)

    if args.ALAN:
        # ALAN Total: Spectral check
        Total_Bb = max(A)
        Tot = np.array([])
        for ii in range(len(ASpec.columns)):
            MAX = max(AI_AS.iloc[:,ii]) # At the surface
            Tot = np.append(Tot, MAX)
        Total_As = np.sum(Tot)
        print("ALAN:\n     Total Broadband signature at surface =", Total_Bb)
        print("     Total Spectral signature (summed) at surface =", Total_As)

#############################################################
#                        Plots 
#############################################################
    location = geo_location
    skycondition = condition
    figurepath = os.getcwd() + f"/Output/{location}_{start_date}_{end_date}_Datum_{datum_percentage}MST"
    if args.plots:
        # Broadband.Broadband(dec_day, A, tide_h, waterdepth, night, sol, aIBT, datum, I, lIBT, phase, Io, IBT, location, skycondition)
        # Tamir_2017.figure_5(AIBT_zR, AIBT_zG, AIBT_zB, zR, zG, zB)
        
        if args.solar:
            Solplot_RGB.SolOverlay(dec_day, SolSpec, SolI_SS, SolI_SSb, Io, tide_h, waterdepth, sol, IBT, datum,datum_percentage, location, figurepath)
            
            # Solplot_RGB.Sol3d(dec_day, SolSpec, SolI_SS, SolI_SSb, SolI_SSRes, datum,datum_percentage)
            # Solplot_RGB.SolRes(dec_day, SolSpec, SolI_SS, SolI_SSb, SolI_SSRes, tide_h, waterdepth, sol, datum, datum_percentage)
            
        if args.lunar:
            Lunplot_RGB.LunOverlay(dec_day, LunSpec, LunI_LS, LunI_LSb, I, tide_h, waterdepth, sol, lIBT, datum,datum_percentage, phase, location, figurepath)

            # Lunplot_RGB.Lun3d(dec_day, LunSpec, LunI_LSb, LunI_LSRes, LunI_LS, datum)
            # Lunplot_RGB.LunRes(dec_day, LunSpec, LunI_LS, LunI_LSb, LunI_LSRes, tide_h, waterdepth, sol, datum)
            
        
        if args.ALAN:
            Aplot_RGB.AOverlay(dec_day, ASpec, AI_AS, AI_ASb, A, tide_h, waterdepth, night, sol, aIBT, col_names_SS, datum, datum_percentage, location, skycondition, ALAN_TYPE, figurepath)
            # Aplot_RGB.A3d(dec_day, ASpec, AI_ASb, AI_ASRes, AI_AS, datum)
            # Aplot_RGB.ARes(dec_day, ASpec, AI_AS, AI_ASb, AI_ASRes, tide_h, waterdepth, sol, datum)
        
        # plt.plot(dec_day, Zc)
        #plt.show()
    ###############################################
    #                   Outputs 
    ###############################################
    #print("Output section still requires updating don't know what the new standard 'TidalLight' model output should be")
    # pdb.set_trace()
    #Dosage = pd.concat([ALAN_datum_dosage_df, Lun_datum_dosage_df, Sol_datum_dosage_df], keys=["ALAN_(uW/m2)", "Lunar_(uW/m2)", "Solar_(uW/m2)"], axis=1)
    #Dosage.to_csv(f"Output/Dosage_{year}_{location}.csv")
    if args.output:
        directory = os.getcwd()+"/Output/"
        if not os.path.exists(directory):
            os.mkdir(directory)
        Output_fname = f"{geo_location}_{start_date}_{end_date}_Datum-{datum_percentage}MST.csv" # OLD DEPTH DESCRIPTION: Depth-{datum}m-above-max-low-tide.csv" # filename of data output
        datapath = os.getcwd() + "/Output/" + Output_fname # path of data file output
        sta = timeit.default_timer()
        df0 = pd.DataFrame({'Location_Lat(degN)' : latitude_deg, 'Location_Lon(degE)' : longitude_deg, 'time_increment(hr)' : t_incr, 'depth_to_datum(m)' : waterdepth, 'modelled_tidal_depth(m)' : ftide_h, 'date' : date_record, 'Jday(decimal)' : dec_day, 'Kd_Blue(1/m)' : kd_blue_t, 'Kd_Green(1/m)' : kd_green_t, 'Kd_Red(1/m)' : kd_red_t, 'Kd_Bb(1/m)' : kd_BB_t})
        frames.append(df0)
        # print("THIS IS A QUESTION!!\n Do you want to add additional detail to output database and ouptut a database of the critical depths?\nY/n?")
        # change_output_name = input()
        # if change_output_name=="Y":
        #     print("type unique identifier for output data")
        #     identifier = input()
        
        if args.solar: 
            # Normalise and store data
            NSolI_SSb = pd.DataFrame(columns=col_names_SS)
            NSolI_SS = pd.DataFrame(columns=col_names_SS)
            NSSb = np.array([])
            NSS = np.array([])

            df1 = pd.DataFrame({'Solar_Alt(deg)' : solar_alt, 'night(0)_day(1)' : sol, 'I_BB(W/m2)' : Io, 'I_BB_Datum(W/m2)' : IBT})
          
            for ss in range(len(col_names_SS)):
                
                NSolI_SS.iloc[:,ss] = (SolI_SS.iloc[:,ss]/max(SolI_SS.iloc[:,ss]))
                NSolI_SSb.iloc[:,ss] = (SolI_SSb.iloc[:,ss]/max(SolI_SSb.iloc[:,ss]))
            #Sresult = pd.concat([df0, df1, NSolI_SS, NSolI_SSb, SolI_SS, SolI_SSb, SolI_SSRes], keys =['','Solar position', 'Surface Normalised', 'Seabed Normalised', 'Surface', 'Seabed', 'Residuals'], axis=1)
            #Sresult = pd.concat([df0, df1, NSolI_SS, NSolI_SSb, SolI_SS, SolI_SSb, SolI_SSRes], axis=1)
            Sresult = pd.concat([df0, df1, SolI_SS.reset_index(drop=True), SolI_SSb.reset_index(drop=True)], axis=1)
            Sresult.to_csv('Output/Solar_' + Output_fname, index=False)
            
        if args.lunar:
            # Normalise and store data
            NLunI_LSb = pd.DataFrame(columns=col_names_SS)
            NLunI_LS = pd.DataFrame(columns=col_names_SS)
            for ll in range(len(col_names_SS)):
                NLunI_LSb.iloc[:,ll] = (LunI_LSb.iloc[:,ll]/max(LunI_LSb.iloc[:,ll]))
                NLunI_LS.iloc[:,ll] = (LunI_LS.iloc[:,ll]/max(LunI_LS.iloc[:,ll]))
            df2 = pd.DataFrame({'Lunar_Alt (deg)' : alt, 'Lunar_Az (deg)': az, 'phase' : phase, 'I_TOA(uW/m2)' : I_atmos, 'I_BB(uW/m2)': I, 'I_BB_datum(uW/m2)' : lIBT})
            #Lresult = pd.concat([df0, df2, NLunI_LS, NLunI_LSb, LunI_LS, LunI_LSb, LunI_LSRes, ], keys= ['', 'Lunar position', 'Surface Normalised', 'Seabed Normalised', 'Surface', 'Seabed', 'Residuals'], axis=1)
            Lresult = pd.concat([df0, df2, LunI_LS.reset_index(drop=True), LunI_LSb.reset_index(drop=True)], axis=1)
            # if change_output_name == "Y":
            #     Output_fname = identifier + Output_fname
            Lresult.to_csv('Output/Lunar_' + Output_fname, index=False)

        if args.ALAN:
            # Normalise and store data
            NAI_ASb = pd.DataFrame(columns=col_names_SS)
            NAI_AS = pd.DataFrame(columns=col_names_SS)
            for aa in range(len(col_names_SS)):
                NAI_ASb.iloc[:,aa] = (AI_ASb.iloc[:,aa]/max(AI_ASb.iloc[:,aa]))
                NAI_AS.iloc[:,aa] = (AI_AS.iloc[:,aa]/max(AI_AS.iloc[:,aa]))   
            #df3 = pd.DataFrame({'night(0)_day(1)' : sol, 'sky_condition': condition, 'Falchi_ALAN(mCd/m^2)': ALAN_mCd, 'ALAN_BB(uW/m^2)' : ALAN_total, 'ALAN_BB_datum(uW/m^2)': aIBT}) 
            df3 = pd.DataFrame({'night(0)_day(1)' : sol, 'sky_condition': condition, 'Falchi_ALAN(mCd/m^2)': ALAN_mCd, 'ALAN_BB(uW/m^2)' : A, 'ALAN_BB_datum(uW/m^2)': aIBT}) 
            #Aresult = pd.concat([df0, df3, NAI_AS, NAI_ASb, AI_AS, AI_ASb, AI_ASRes], keys= ['','ALAN', 'Surface Normalised', 'Seabed Normalised', 'Surface', 'Seabed', 'Residuals'], axis=1)
            Aresult = pd.concat([df0, df3, AI_AS.reset_index(drop=True), AI_ASb.reset_index(drop=True)], axis=1)
            Output_fname = f"{geo_location}_{sky_condition}_{start_date}_{end_date}_Datum-{datum_percentage}MST.csv" # filename of data output (MST = Mean Spring Tide)
            # print("THIS IS A QUESTION!!\n Do you want to add additional detail to output database and ouptut a database of the critical depths?\nY/n?")
            if args.station:
                Output_fname = str(args.station) + "_" + Output_fname
                d = {"ALAN_R_Depth[620-740nm](uW/m2)": A_zR, "R_Depth(m)": zR, "ALAN_G_Depth[495-560nm](uW/m2)": A_zG, "G_Depth(m)": zG, "ALAN_B_Depth[400-500nm](uW/m2)": A_zB, "B_Depth(m)": zB,  "ALAN_Bb_Depth(uW/m2)" : A_zBb, "Kd_Red": KD[0], "Kd_Green": KD[1],"Kd_Blue": KD[2], "crit_depth_R(m)" : critical_depth_R[0], "crit_depth_G(m)" : critical_depth_G[0], "crit_depth_B(m)" : critical_depth_B[0], "crit_depth_Bb(m)" : critical_depth_Bb[0]} # nm range taken from (Davies et al., 2020) the source of the Hydrolight dataset
                Tamir_comparison = pd.DataFrame({k:pd.Series(v) for k,v in d.items() })
                Tamir_comparison.to_csv(f'Output/ALAN_TAMIR_COMPARISON_i' + Output_fname)
            Aresult.to_csv(f'Output/ALAN_' + Output_fname, index=False)
        sto = timeit.default_timer()
        print('Data output completed in', (sto-sta)/60, 'minutes')
        #plt.show()
    
###############################################
#      Thank you for using TidalLight
#   Authors: Adam E Wright, Tim Smyth, ...
###############################################  
  
# Run script if called from command line.   
if __name__=='__main__':
    main()
    #plt.show()
  
        
       
                     





