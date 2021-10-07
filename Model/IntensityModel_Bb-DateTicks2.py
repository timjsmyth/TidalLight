#!/usr/bin/env python3

# Authour: Adam Wright
# Purpose: Calculate intensity of lunar, solar and ALAN brightness at the seasurface and seabed
# This script only contains a broadband Irradiance model and does not split spectrally
# a number of assumtions are made in the calculation process such as: Both the Sun and Moon are point sources, Lunar albedo is constant throughout cycle.
# Sky condition is constant for the entirety of the modelled time, scattering is constant through the atmosphere
# twilight is not factored in.
# Tidal data must be in the form of 'BODCTideData'.csv file 

# Import packages
from pysolar.solar import *
import pdb
from pysolar.radiation import *
import utide
import matplotlib.dates as mdates

import numpy as np
import datetime
import matplotlib
# matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import argparse

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_moon, Angle
from astroplan import moon, Observer
import timeit
import calendar

import os
import subprocess
from dateutil import parser





timeline = input("Do you wish to enter a date range [Y/n]:\n"
                 "Y: enter date range for model\n"
                 "n: default to last 7 days \n")

if timeline == "Y":
    # ask the user to input date range
    start_date_input = input("Please enter start date: ")
    end_date_input = input("Please enter end date: ")
else:
    # otherwise just run past 7 days
    start_date_input = datetime.datetime.utcnow() - datetime.timedelta(days=7)
    end_date_input = datetime.datetime.utcnow()

data_start_date = subprocess.getoutput(f"date --date '{start_date_input}' +'%F %H:%M:%S'")
#print("Start date", data_start_date)
data_end_date = subprocess.getoutput(f"date --date '{end_date_input}' +'%F %H:%M:%S'")
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
tdiff_day = tt_e-tt_s
year = date_start.year


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
    sta = timeit.default_timer() # Start timer for entire function
##  Create Parser instructions
    parser = argparse.ArgumentParser()
    # TIDE
    parser.add_argument("-t", "--tidal", action="store_true", help="Tidal model - REQUIRED: TEXT FILE CONTAINING TIDAL DATA FROM COASTAL OBSERVATORY")
    parser.add_argument("-d", "--datum", action="store", type=float, help="Depth in metres to a datum of choice, this is a reference position in the water column below the lowest astronomical tide of the given year")
    # LIGHT SOURCES
    parser.add_argument("-s", "--solar", action="store_true", help="Solar Irradiance at sea level")
    parser.add_argument("-A", "--ALAN", action="store", type=float, help="0 = clear condition (2.88 PFFD), 1 = cloudy condition (6.24 PFFD)")
    parser.add_argument("-l", "--lunar", action="store_true", help="Lunar information dday, Alt, Az, distance") 
    # MODEL PARAMETERS
    parser.add_argument("-T", "--time", action="store", type=float, help="Time increment in decimal hours i.e. 0.25 = 15mins")
    parser.add_argument("-lat", "--latitude",type=float, action="store", help="Longitude in decimal degrees")
    parser.add_argument("-lon", "--longitude", type=float, action="store", help="Latitude in decimal degrees")
    # MODEL OUTPUTS
    parser.add_argument("-p", "--plots", action="store",type=int, help="plot: 1=Solar, 2=Lunar, 3=ALAN, 4=ALL ")
    parser.add_argument("-o", "--output", action="store_true", help="generate output file(s)")
    parser.add_argument("-S", "--save", action="store_true", help="add annotation to plot and save figure")
    args = parser.parse_args()
    if args.save:
        print('annotation to be added to figure path:') 
        annotation = input()
##  Define constants, components and arrays    
    Isc = 2400 # Solar constant (umol/m^2/s)
    k_atmos = 0.01 # Atmospheric diffusion coefficient set from Roberts et al. (2018)
    latitude_deg = 50.25 # latitude_deg
    longitude_deg = -4.13 # longitude_deg
    
    t_incr = args.time
    #year = args.year
    location = EarthLocation(lat=50*u.deg, lon=longitude_deg*u.deg, height=1*u.m)
    location_eq = EarthLocation(lat=0*u.deg, lon=longitude_deg*u.deg, height=1*u.m)
    sol = []; night = []; z = []; dec_day = []; Io = []; IBT = []; sol_eq = []; I_atmos = []
    alt = []; az = []; phase = []; I = []; lIBT = []; PhaseAng_EQ = []; PhaseAng = []; Phase_lit = []
    RightAscension = []; Declination = []; A = []; aIBT = []
    TL = []; t = []; waterdepth = []; tide_h = []
    frames = []; Dates = []
    alt_EQ = []; az_EQ = []; phase_EQ = []; I_EQ = []
       
##  ALAN average broadband intensity for clear and cloudy sky conditions
    if args.ALAN==1:
        condition = 'cloudy'
        avgPFFD = 6.24 # $\mu$mol m^-2 s^1
    elif args.ALAN==2:
        condition = 'clear'
        avgPFFD = 2.88 # $\mu$mol m$^{-2}$ s$^{-1}$


    if (args.tidal):
        datum = args.datum

##      Read in tidal data from BODC
        # This will change given the format and stored style of each tidal data file
        print('Running tidal model...')
        Tide_data_fname = "DEV19.txt" # filename of tidal data
        tidepath = os.getcwd() + "/Required/" + Tide_data_fname # file path of tidal data
        tides_ = pd.read_csv(tidepath, skiprows=11, delimiter='    ', engine='python') #usecols=np.arange(16,48), engine='python')
        df = pd.DataFrame(tides_)
        for i in range(len(tides_)):

##          Convert/alter DataFrame Variables to type datetime and float
            # Strip variables of additional Letters numbers and blankspaces NOTE this is explicitly for the BODC Tidal data for Plymouth
            string = str(df.iloc[i,1])
            string0 = str(df.iloc[i,0])
            string0 = string0[-19:]
            string0 = string0.lstrip()
            string = string.replace("M", "")
            string = string.replace("T", "")
            string2 = str(df.iloc[i,2])
            string2 = string2.replace("M", "")
            string2 = string2.replace("T", "")
            tidelevel = float(string)

            # Convert date to datetime format
            T = datetime.datetime.strptime(string0, '%Y/%m/%d %H:%M:%S')

            # Append variables to lists
            TL.append(tidelevel)
            t.append(T)

        # Convert time back to DataFrame
        t_df = pd.to_datetime(t)
        # Convert time to format of UTide ##'date2num' function only seems to work with pandas DF ##
        time = mdates.date2num(t_df.to_pydatetime())
        tide = np.array(TL, dtype=float) # Tidal data array for model input
        # Create Tidal model 'c'
        c = utide.solve(time, u=tide, v=None, lat=latitude_deg,
                        nodal=False,
                        trend=False,
                        method='ols',
                        conf_int='linear',
                        Rayleigh_min=0.95)
        print('finished tidal model...')
    
    
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
    # predict azimuth, altitude, and distance from ICRS reference frame (International Celectial Reference Frame)
    # Lunar phase calculated as decimal where 1 = full-moon
    # 
    # TO BE DONE: calculate solar reflection intensity employing a mean lunar albedo of 0.12 
    # NOTE: Lunar albedo is wavelength dependent the albedo for specific wavelengths will be important for accuracy and validity of future spectral decomp of intensity
    # Plan is to use pysolar airmass function with the output of astropy altitude for consistency

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                                           TIDAL MODULE: args.seabed
    #
    # Models tidal range for a given location using the UTide function which creates the relevant harmonic coefficients
    # for a data set from a coastal observatory with time-series data. The data for the UK is collected from the BODC website.
    # Data file must be in the same location as the scirpt*
    # Roberts et al., (2018) describes the proccess through which intensity is attenuated through the water columnn.

    #
    # Create time-loop
    #
    
    sta = timeit.default_timer() # Function requires ~55mins-95mins
        
    for j in range(0, tdiff_day): #Establish for loop for days 
        jday = tt_s+j # set counter of days between input date range 
        Js = str(jday) # Convert to string for parsing
        Ys = str(year) # Convert to string for parsing
        YsJs = Ys + Js # Adjoin strings
        ddt = datetime.datetime.strptime(YsJs,'%Y%j').date() # Create Python datetime in preparation for for loop
        day = ddt.day # Get Python datetime Day
        month = ddt.month # Get Python datetime month

        Sta = timeit.default_timer() # Start Timer
        for h in range(0,24): # Establish for loop for hours
            for m in np.arange(0, 1, t_incr): # Establish for loop for minutes
                H = ((h+m)/24) # Create Decimal component of the day 
                m = int(m*60) # Convert decimal input to minutes for datetime processing
                if m>=60: # IF input is put in incorrectly then default to hourly iterations
                    m = 0 
                date = datetime.datetime(year, month, day, h, m, tzinfo=datetime.timezone.utc) # Create Python datetime variable for use in iterable calculations, All in UTC
                Dates.append(date) # ppend dates to list for plotting secondary x-axis opposed to Serial day 
                jday = date.timetuple().tm_yday # Convert to serial day
                dday = int(jday)+H # Add decimal component to serial day
                dec_day.append(dday) # Append decimal day to list for plotting
                altitude_deg = get_altitude(50, longitude_deg, date) # Get solar altitude based on geodetic position 
                altitude_deg_eq = get_altitude(0.1, longitude_deg, date) # Get Solar altitude for the equator
                solar_day = dday*float(1)
                z.append(-299)
                if altitude_deg<=0: # Define Nighttime
                    solar_day = float(0) # Set Nighttime value
                    Night = 2700
                elif -6<=altitude_deg<0: # Define Civil Twilight
                    solar_day = float(0.75) # Set Civil Twilight value
                elif -12<=altitude_deg<-6: # Define Nautical Twilight
                    solar_day = float(0.5) # Set Nautical Twilight value
                elif -18<=altitude_deg<-12: # Define Astronomical Twilight
                    solar_day = 0.25 # Set Astronomical Twilight value
                else:  # Define DayTime 
                    solar_day = float(1) # Set DayTime value
                    Night = -600
                sol.append(solar_day) # Append Solar day for DATA output 
                night.append(Night) # Night values are assigned arbitrarily for plotting effect
                ############# REPEAT ABOVE STEPS FOR EQ #######################
                if altitude_deg_eq<=0:
                    solar_day_eq = float(0)
                    Night = 2700
                # elif -6<=altitude_deg_eq<0: ################ THIS IS UNCOMMENTED IN THE SCRIPT
                #     solar_day_eq = float(0.75) ############# IntensityModel_Bb_Twilight.py  
                # elif -12<=altitude_deg_eq<-6:
                #     solar_day_eq = float(0.5)
                # elif -18<=altitude_deg_eq<-12:
                #     solar_day_eq = 0.25
                else:
                    solar_day_eq = float(1)
                sol_eq.append(solar_day_eq)
                ################################################################
            
                # SOLAR CALCULATION 
                if (args.solar):
                    ##  Calculate variable using Pysolar
                    airmass = get_air_mass_ratio(altitude_deg) # Airmass using altitude = 1/sin(altitude(radians))
                    ##  Calculate Intensity                 
                    Iatmos = Isc*(1+0.0344*np.cos(np.deg2rad(360.*float(dday)/365)))
                    if solar_day==0:
                        IO = float(0)
                    else:
                        IO = Iatmos*np.sin(np.deg2rad(altitude_deg))*np.exp(-1.*airmass*k_atmos)
                    Io.append(IO)

                # ALAN CALCULATION    
                if args.ALAN:
                    if solar_day==1:
                        Aint = 0
                    else:
                        Aint = avgPFFD
                    A.append(Aint)
                    
                # LUNAR CALCULATIONS
                if args.lunar:
                    Lunar_albedo = 0.16
                    dt = Time(date) # Time within model
                    EarthMoon = get_moon(dt, location) # Relative positioning              
                    moon_icrs = EarthMoon.transform_to('icrs') # Relative Ephemeris calculation    
                    moonaltaz = moon_icrs.transform_to(AltAz(obstime=dt, location=location)) # Transform positioning to Alt, Az
                    alt_ = float(Angle(moonaltaz.alt).degree) # Lunar Altitude 
                    az_ = float(Angle(moonaltaz.az).degree) # Lunar Azimuth                     
                    Phase = float(moon.moon_illumination(dt)) # Lunar Phase 
                    Phase_lb = lumme_bowell(moon.moon_phase_angle(dt).value) # Lunar Phase function - Lumme-Bowell
                    #albedo_phased = Lunar_albedo*Phase 
                    albedo_phased = Lunar_albedo*Phase_lb 
                    # LUNAR REFLECTION CALC. 
                    Lunar_refl = (Isc/np.pi)*albedo_phased*0.000064692 # 0.000064692 = sr value for subtending angle of the lunar disc (0.26 degree semi-diameter --> sr)
                    #airmass = get_air_mass_ratio(alt_) ########## diffusivity through airmass is wavelength specific - currently using broadband attenuation for each wvlnth
                    airmass = get_air_mass_kasten_young(alt_)
                                  


                    if Phase_lb<7:
                        Lunar_refl = Lunar_refl*1.35 # Lunar reflectance is brighter during a full moon when the phase angle < 7(Krisciunas & Schaefer, 1991), A Model of the Brightness of the moon

                    if altitude_deg>0:
                        Iol = 0
                    elif alt_<0: # Define Moonset and Moonrise
                        Iol = 0
                    
                    # elif 0<solar_day<1: # Twilight modulation of Intensity
                    #     Iol = (Lunar_refl*np.sin(np.deg2rad(alt_))*np.exp(-1.*airmass*k_atmos))*(1-solar_day) # Calculate Lunar Intensity at surface using Beers law. Modulated by twilight factor (solar_day)              
                    else: # Night
                        Iol = (Lunar_refl*np.sin(np.deg2rad(alt_))*np.exp(-1.*airmass*k_atmos)) # Beers Law (Roberts et al., 2018 [Eq.9]) 
                    alt.append(alt_)
                    PhaseAng.append(Phase_lb) # Append Phase angle in degrees
                    phase.append(Phase) # Append Noramlised phase/ % of moon Illuminated 
                    az.append(az_) # Append Lunar Azimuth
                    I.append(Iol) # Append Lunar Intensity at sea level

                    ################## REPEAT FOR EQ #######################

                    Lunar_albedo = 0.16
                    dt = Time(date) # Time within model
                    EarthMoon_eq = get_moon(dt, location_eq) # Relative positioning              
                    moon_icrs_eq = EarthMoon_eq.transform_to('icrs') # Relative Ephemeris calculation    
                    moonaltaz_eq = moon_icrs_eq.transform_to(AltAz(obstime=dt, location=location_eq)) # Transform positioning to Alt, Az
                    alt_eq = float(Angle(moonaltaz_eq.alt).degree) # Lunar Altitude 
                    az_eq = float(Angle(moonaltaz_eq.az).degree) # Lunar Azimuth                     
                    Phase_eq = float(moon.moon_illumination(dt)) # Lunar Phase 
                    Phase_lb_eq = lumme_bowell(moon.moon_phase_angle(dt).value) # Lunar Phase function - Lumme-Bowell
                    #albedo_phased = Lunar_albedo*Phase 
                    albedo_phased_eq = Lunar_albedo*Phase_lb_eq 
                    # LUNAR REFLECTION CALC. 
                    Lunar_refl_eq = (Isc/np.pi)*albedo_phased_eq*0.000064692 # 0.000064692 = sr value for subtending angle of the lunar disc (0.26 degree semi-diameter --> sr)
                    #airmass = get_air_mass_ratio(alt_) ########## diffusivity through airmass is wavelength specific - currently using broadband attenuation for each wvlnth
                    airmass_eq = get_air_mass_kasten_young(alt_eq)
                   
                    if Phase_lb_eq<7:
                        Lunar_refl_eq = Lunar_refl_eq*1.35
                    if altitude_deg_eq>0:
                        Iol_eq =0
                    elif alt_eq<0:
                        Iol_eq = 0
                    # elif 0<solar_day_eq<1:
                    #     Iol_eq = (Lunar_refl_eq*np.sin(np.deg2rad(alt_eq))*np.exp(-1*airmass_eq*k_atmos))*(1-solar_day_eq) ### Solar_day_eq = Daytime, Twilight, Nightime multiplication factor
                    else:
                        Iol_eq = (Lunar_refl_eq*np.sin(np.deg2rad(alt_eq))*np.exp(-1*airmass_eq*k_atmos))   

                    PhaseAng_EQ.append(Phase_lb_eq)
                    alt_EQ.append(alt_eq)
                    phase_EQ.append(Phase_eq)
                    az_EQ.append(az_eq)
                    I_EQ.append(Iol_eq)
                    
                #####################################################################

                # TIDAL CALCULATIONS    
                if args.tidal:
                    date2 = datetime.datetime(year, month, day, h, m)
                    date2 = str(date2)
                    fmt_date = datetime.datetime.strptime(date2, '%Y-%m-%d %H:%M:%S')
                    tidetime = pd.to_datetime(fmt_date)
                    tide_time = mdates.date2num(tidetime.to_pydatetime())
                    reconst = utide.reconstruct(tide_time, c)
                    tide_h.append(reconst.h[0])
                    depth_to_datum = float(reconst.h-datum) # depth = height of water column above datum
                    if reconst.h<=datum:
                        depth_to_datum = float(0)
                    waterdepth.append(depth_to_datum)
                    kPAR = ((0.5+0.5*np.cos(2*np.pi*dday)/365))*0.1*reconst.h+0.4
                    # ATTENUATION OF INTENSITY TO SEABED
                    # SOLAR
                    if args.solar:
                        if depth_to_datum==0:
                            iBT = float(IO)
                        else:
                            iBT = float(IO*np.exp(-kPAR*depth_to_datum))
                        IBT.append(iBT)

                    # LUNAR 
                    if args.lunar:
                        if depth_to_datum==0:
                            LIBT = float(Iol)
                        else:
                            LIBT = float(Iol*np.exp(-kPAR*depth_to_datum))
                        lIBT.append(LIBT)

                    # ALAN
                    if args.ALAN:
                        if depth_to_datum==0:
                            AIBT = float(Aint)
                        else:
                            AIBT = float(Aint*np.exp(-kPAR*depth_to_datum))
                        aIBT.append(AIBT)
        Sto = timeit.default_timer()
        print('Calculations completed for ', date, ' in ', (Sto-Sta)/60, 'minutes')

    sto = timeit.default_timer()
    print('Calculations completed in ', (sto-sta)/60, 'minutes')

    # OUTPUT DATA 
    if args.output:
        Output_fname = f"DATA_{tt_s}-{tt_e}_Bb.csv" # filename of data output
        if args.save:
            Output_fname = f'f"DATA_{tt_s}-{tt_e}_Bb_{annotation}.csv'
        datapath = os.getcwd() + "/Output/" + Output_fname # path of data file output
        sta = timeit.default_timer()
        df0 = pd.DataFrame({'Location_Lat (deg)' : 50, 'Location_Lat_EQ (deg)' : 0,  'Location_Long (deg)' : longitude_deg, 'binary_day' : sol, 'Time increment (Decimal hours)' : args.time, 'decimal_day' : dec_day, 'date' : Dates})
        frames.append(df0)

        if args.solar:
            sta = timeit.default_timer()
            N_Io = [Io[i]/max(Io) for i in range(len(Io))]
            if args.tidal: 
                N_IBT = [IBT[i]/max(IBT) for i in range(len(IBT))]
                df1 = pd.DataFrame({'Solar_Alt(deg)' : altitude_deg, 'Solar_Int_SeaLevel(PPFD)' : Io, 'Solar_Int_SeaLevel_Norm(PPFD)' : N_Io, 'Solar_Alt_EQ(deg)' : altitude_deg_eq, 'Solar_Int_SeaLevel_EQ(PPFD)' : Io_EQ, 'Solar_Int_SeaLevel_Norm_EQ(PPFD)' : N_Io_EQ, 'Int_BelowTide (PFFD)' : IBT, 'Int_BelowTide_Norm(PFFD)' : N_IBT,})   
            else:
                df1 = pd.DataFrame({'Solar_Alt(deg)' : altitude_deg, 'Solar_Int_SeaLevel(PPFD)' : Io, 'Solar_Int_SeaLevel_Norm(PPFD)' : N_Io, 'Solar_Alt_EQ(deg)' : altitude_deg_eq, 'Solar_Int_SeaLevel_EQ(PPFD)' : Io_EQ, 'Solar_Int_SeaLevel_Norm_EQ(PPFD)' : N_Io_EQ})
            frames.append(df1)

        if args.lunar:
            N_I = [I[i]/max(I) for i in range(len(I))]
            N_I_EQ = [I_EQ[i]/max(I_EQ) for i in range(len(I_EQ))]
            if args.tidal:
                N_lIBT = [lIBT[i]/max(lIBT) for i in range(len(lIBT))]
                df2 = pd.DataFrame({'Lunar_Alt(deg)' : alt, 'Lunar_Az(deg)': az, 'phase(normalised)' : phase, 'phase_label' : Phase_lit ,'Lunar_Int_SeaLevel(PPFD)' : I, 'Lunar_Int_SeaLevel_Norm(PPFD)' : N_I, 'Lunar_Alt_EQ(deg)' : alt_EQ, 'Lunar_Az_EQ(deg)': az_EQ, 'phase_EQ(normalised)' : phase_EQ, 'Lunar_Int_SeaLevel_EQ(PPFD)' : I_EQ, 'Lunar_Int_SeaLevel_Norm_EQ(PPFD)' : N_I_EQ, 'Phase_Angle' : PhaseAng, 'Phase_Angle_EQ' : PhaseAng_EQ, 'Lunar_Int_BelowTide(PPFD)' : lIBT, 'Lunar_Int_BelowTide_Norm(PPFD)' : N_lIBT})   
            else:
                df2 = pd.DataFrame({'RightAsc' : RightAscension, 'Decl' : Declination, 'Lunar_Alt(deg)' : alt, 'Lunar_Az(deg)': az, 'phase(normalised)' : phase, 'Lunar_Int_SeaLevel(PPFD)' : I, 'Lunar_Int_SeaLevel_Norm(PPFD)' : N_I, 'Lunar_Alt_EQ(deg)' : alt_EQ, 'Lunar_Az_EQ(deg)': az_EQ, 'phase_EQ(normalised)' : phase_EQ, 'Lunar_Int_SeaLevel_EQ(PPFD)' : I_EQ, 'Lunar_Int_SeaLevel_Norm_EQ(PPFD)' : N_I_EQ, 'Phase_Angle' : PhaseAng, 'Phase_Angle_EQ' : PhaseAng_EQ })   

            frames.append(df2)

        if args.ALAN:
            N_aIBT = [aIBT[i]/max(aIBT) for i in range(len(aIBT))]
            df3 = pd.DataFrame({'ALAN(PFFD)' : A, 'ALAN_BelowTide(PFFD)' : aIBT, 'ALAN_BelowTide_Norm(PFFD)' : N_aIBT})
            frames.append(df3)

        if args.tidal:
            df4 = pd.DataFrame({'water_height_above_datum(m)' : waterdepth, 'tidal_range' : tide_h,  'depth_to_datum' : args.datum})
            frames.append(df4)

        sto = timeit.default_timer()
        print('Data output completed in', (sto-sta)/60, 'minutes')
        result = pd.concat(frames, axis=1)
        result.to_csv(Output_fname)
        
        print(result)

        
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #                                           PLOTS MODULE
        #
        # Plot data based on input from command line
        #
        # 
                        
    if args.plots:
        startdate = int(tt_s)
        enddate = int(tt_e)
        plot_index=args.plots
        figure_name = (f'JDay{startdate}-{enddate}_{year}_Bb.png')
        if args.save:
            figure_name = (f'JDay{startdate}-{enddate}_{year}_Bb_{annotation}.png')
        figurepath = os.getcwd() + "/Output/" + figure_name
        if plot_index==5: 
            
            # This plot is for the publication
            # plot the full year of lunar data for lat = 50 and lat = 0.1 equator
            print('Lat 50 max', max(I))
            print('Lat EQ max', max(I_EQ))
            fig, ax = plt.subplots(2, figsize=(13,7.34))
            ax[0].plot(dec_day, I, color='black')
            ax[0].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[0].set_title('Lunar Irradiance at latitude 50$^\circ$N')
            # ax[0].set_xlim([startdate, enddate])

            ax[1].plot(dec_day, I_EQ, color='black', alpha=0.7)
            ax[1].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[1].set_title('Lunar Irradiance at the equator')
            #ax[1].set_xlim([startdate, enddate])
            
            ax[1].set_xlabel(f'Serial Day {year}')
            months = mdates.MonthLocator()
            months_fmt = mdates.DateFormatter('%B')
            ax1 = ax[1].twiny()
            ax1.plot(Dates, I_EQ, color='white', alpha=0.1)
            # ax1.set_xticks(dec_day)
            # ax1.set_xticklabels(Dates)
            ax1.xaxis.set_major_locator(months)
            ax1.xaxis.set_major_formatter(months_fmt)
            ax1.set_xlabel('Month')
            ax1.xaxis.set_label_position('bottom')
            ax1.xaxis.set_ticks_position('bottom')
            ax1.spines['bottom'].set_position(('outward',36))
            plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.3)
            if args.save:
                fig.savefig(figurepath)
            plt.show()
        

        if plot_index==6: 
            
            # This plot is for the publication
            # plot the full year of lunar data for lat = 50 and lat = 0.1 equator
            print('Lat 50 max', max(I))
            print('Lat EQ max', max(I_EQ))
            fig, ax = plt.subplots(4, figsize=(13,7.34))
            ax[0].plot(dec_day, I, color='black')
            ax[0].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[0].set_title('Lunar Irradiance at latitude 50$^\circ$N')
            # ax[0].set_xlim([startdate, enddate])
            ax[1].plot(dec_day, alt, color='black', linestyle='dashed')
            ax[1].set_ylabel('Altitude ($^\circ$)')
            ax[1].set_title('Lunar Altitude at 50$^\circ$N')

            ax[2].plot(dec_day, I_EQ, color='black', alpha=0.7)
            ax[2].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[2].set_title('Lunar Irradiance at the equator')
            #ax[1].set_xlim([startdate, enddate])
            
            ax[3].plot(dec_day, alt_EQ, color='black', alpha=0.7, linestyle='dashed')
            ax[3].set_xlabel(f'Serial Day {year}')
            ax[3].set_title('Lunar Altitude at the Equator')
            months = mdates.MonthLocator()
            months_fmt = mdates.DateFormatter('%B')
            ax1 = ax[3].twiny()
            ax1.plot(Dates, alt_EQ, color='white', alpha=0.1)
            # ax1.set_xticks(dec_day)
            # ax1.set_xticklabels(Dates)
            ax1.xaxis.set_major_locator(months)
            ax1.xaxis.set_major_formatter(months_fmt)
            ax1.set_xlabel('Month')
            ax1.xaxis.set_label_position('bottom')
            ax1.xaxis.set_ticks_position('bottom')
            ax1.spines['bottom'].set_position(('outward',36))
            plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.3)
            if args.save:
                fig.savefig(figurepath)
            plt.show()
            
        if plot_index==1:
            print('trying to plot')
            fig, ax = plt.subplots(4, figsize=(13,7.34))
            AA = np.ones(len(dec_day), dtype=int)
            aa = 2300*AA
            bb = -300*AA
            # SOLAR INT SEALEVEL
            ax[0].plot(dec_day, Io, color='firebrick')
            ax[0].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[0].set_title('Solar Irradiance at sea level')
            ax[0].set_xlim([startdate, enddate])
            ax0 = ax[0].twinx()
            ax0.set_ylim([(-max(Io)/10), (max(Io)+max(Io)/10)])
            ax0.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
            ax0.axes.get_yaxis().set_visible(False)

            # SOLAR INT SEABED
            ax1 = ax[1].twinx()
            ax1.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
            ax1.axes.get_yaxis().set_visible(False)
            ax[1].plot(dec_day, IBT, color='indianred')
            ax[1].set_title('Solar Irradiance seabed')
            ax[1].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[1].set_xlim([startdate, enddate])
            ax1.set_ylim([(-max(IBT)/10), (max(IBT)+max(IBT)/10)])
            
            # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
            ax2i = ax[2].twinx()
            ax2i.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
            ax2i.axes.get_yaxis().set_visible(False)
            ax2i.set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])
            ax[2].plot(dec_day, waterdepth, color='royalblue')
            ax[2].set_title(f'Height of water column above datum at {args.datum}m')
            ax[2].set_ylabel('Water column\n above datum (m)')
            ax2 = ax[2].twinx()
            ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.5, linestyle='dashed', label='datum')
            ax2.legend(prop={"size":8}, loc='upper right')
            ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.5)
            ax2.set_ylabel('Tidal range (m)', color='cadetblue')
            ax[2].set_xlim([startdate, enddate])
            ax[2].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])

            # DAY/NIGHT
            ax[3].set_title('Daylight hours')
            ax[3].plot(dec_day, sol, color='gold')
            ax[3].set_ylabel('Daylight')
            ax[3].set_xlabel(f"Serial Day {year}")
            ax[3].set_xlim([startdate, enddate])
            loc = mdates.WeekdayLocator(byweekday=mdates.MO)
            loc_fmt = mdates.DateFormatter('%F')
            ax1 = ax[3].twiny()
            ax1.plot(Dates, sol, color='white', alpha=0.1)
            # ax1.set_xticks(dec_day)
            # ax1.set_xticklabels(Dates)
            ax1.xaxis.set_major_locator(loc)
            ax1.xaxis.set_major_formatter(loc_fmt)
            ax1.set_xlabel('Date')
            ax1.xaxis.set_label_position('bottom')
            ax1.xaxis.set_ticks_position('bottom')
            ax1.spines['bottom'].set_position(('outward',36))
            plt.setp(ax1.xaxis.get_majorticklabels(), rotation=20, fontsize=8)
            plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.3)
            if args.save:
                fig.savefig(figurepath)
            print('made it to here')
            plt.show()
            
        if plot_index==2:
            fig, ax = plt.subplots(4, figsize=(13,8))
            loc = mdates.DayLocator(bymonthday=None, interval=3, tz=None)
            loc_minor = mdates.DayLocator(bymonthday=None, interval=1, tz=None)
            loc_fmt = mdates.DateFormatter('%d-%B')
            AA = np.ones(len(dec_day), dtype=int)
            aa = 2300*AA
            bb = -300*AA
            # LUNAR INT SEALEVEL
            ax[0].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[0].set_title('Lunar irradiance at sea level')
            ax[0].plot(dec_day, I, label='Lunar irradiance at sea level', color='sienna')
            ax0 = ax[0].twinx()
            ax0.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
            ax0.axes.get_yaxis().set_visible(False)
            # ax[0].set_xticklabels([])
            # ax[0].axes.get_xaxis().set_visible(False)
            ax[0].set_xlim([startdate, enddate])
            ax[0].set_ylim([(-max(I)/10), (max(I)+max(I)/10)])
            ax0.set_xticklabels([])
            ax0.set_xticks([])
            ax[0].set_xticklabels([])
            ax[0].set_xticks([])
            ax0i = ax[0].twiny()
            ax0i.plot(Dates, I, color='none')
            # ax0i.set_xticks(Dates)
            ax0i.xaxis.set_major_locator(loc)
            ax0i.xaxis.set_major_formatter(loc_fmt)
            ax0i.xaxis.set_minor_locator(loc_minor)
            ax0i.xaxis.set_label_position('bottom')
            ax0i.xaxis.set_ticks_position('bottom')
            ax0i.set_xticklabels([])

            # LUNAR INT SEABED 
            ax1 = ax[1].twinx()
            ax1.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
            ax1.axes.get_yaxis().set_visible(False)
            ax1.set_xticklabels([])
            ax1.set_xticks([])
            ax[1].set_xticklabels([])
            ax[1].set_xticks([])
            # ax[1].axes.get_xaxis().set_visible(False)
            pdb.set_trace()
            ax[1].plot(dec_day, lIBT, label='Lunar irradiance at seabed', color='lightsalmon')
            ax[1].set_title('Lunar irradiance at seabed ')
            ax[1].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[1].set_xlim([startdate, enddate])
            ax[1].set_ylim([(-max(lIBT)/10), (max(lIBT)+max(lIBT)/10)])


            ax1i = ax[1].twiny()
            ax1i.plot(Dates, lIBT, color='none', alpha=0.1)
            # ax1.set_xticks(dec_day)
            # ax1.set_xticklabels(Dates)
            ax1i.xaxis.set_major_locator(loc)
            ax1i.xaxis.set_major_formatter(loc_fmt)
            ax1i.xaxis.set_label_position('bottom')
            ax1i.xaxis.set_ticks_position('bottom')
            ax1i.xaxis.set_minor_locator(loc_minor)
            ax1i.set_xticklabels([])


            # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
            ax2i = ax[2].twinx()
            ax2i.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
            ax2i.axes.get_yaxis().set_visible(False)
            ax[2].plot(dec_day, waterdepth, label='water depth', color='royalblue')
            ax[2].set_ylabel(f'Water column\n above {args.datum}m datum')
            ax2 = ax[2].twinx()
            ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.5, linestyle='dashed', label='datum')
            ax2.legend(prop={"size":8}, loc='upper right')
            ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.5)
            ax2.set_ylabel('Tidal range (m)', color='cadetblue')
            ax2.set_xticklabels([])
            ax2.set_xticks([])
            ax[2].set_xticklabels([])
            ax[2].set_xticks([])
            # ax[2].axes.get_xaxis().set_visible(False)
            ax[2].set_xlim([startdate, enddate])
            ax[2].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])
            ax[2].set_title('Height of water column')              
            
            ax2j = ax[2].twiny()
            ax2j.plot(Dates, sol, color='none')
            ax2j.set_xticks(Dates)
            ax2j.xaxis.set_major_locator(loc)
            ax2j.xaxis.set_minor_locator(loc_minor)
            ax2j.xaxis.set_major_formatter(loc_fmt)
            ax2j.xaxis.set_label_position('bottom')
            ax2j.xaxis.set_ticks_position('bottom')
            ax2j.set_xticklabels([])
            # LUNAR PHASE
            ax[3].plot(dec_day, phase, label='Lunar phase', color='darkgrey')
            ax[3].set_title('Lunar phase cycle')
            ax[3].set_ylabel('Normalised phase')
            # ax[3].set_xlabel("Date")
            ax[3].set_xticklabels([])
            ax[3].set_xticks([])
            
            ax[3].set_xlim([startdate, enddate])
            
            ax1 = ax[3].twiny()
            ax3 = ax[3].twinx()
            ax3.plot(Dates, Phase_lit, color='none', alpha=0.1)
            plt.gca().invert_yaxis()
            # ax3.axes.get_xaxis().set_visible(False)
            ax3.set_xticklabels([])
            ax3.set_xticks([])
            ax1.plot(Dates, phase, color='none', alpha=0.1)
            # ax1.set_xticks(dec_day)
            # ax1.set_xticklabels(Dates)
            ax1.xaxis.set_major_locator(loc)
            ax1.xaxis.set_major_formatter(loc_fmt)
            ax1.xaxis.set_minor_locator(loc_minor)
            ax1.set_xlabel("Date")
            ax1.xaxis.set_label_position('bottom')
            ax1.xaxis.set_ticks_position('bottom')
            # ax1.spines['bottom'].set_position(('outward',36))
            plt.tight_layout() #pad=0.5, w_pad=0.5, h_pad=0.5
            plt.setp(ax1.xaxis.get_majorticklabels(), rotation=25, fontsize=8)
            if args.save:
                filepath = os.getcwd() + figure_name
                fig.savefig(filepath, dpi=450)
            plt.show()
            
        if plot_index==3:
            fig, ax = plt.subplots(3, figsize=(13,7.34))

            # ALAN INT SEALEVEL
            ax[0].set_ylabel('Irradaince\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[0].plot(dec_day, A, label=f'ALAN at sea level ({condition} sky PAR average = {avgPFFD} PFFD)', color='seagreen')
            ax[0].set_title(f'ALAN Irradiance at sea level PAR average {avgPFFD} PFFD - {condition}')
            ax[0].set_xlim([startdate, enddate])

            # ALAN INT SEABED
            ax[1].set_title('ALAN Irradiance at seabed')
            ax[1].plot(dec_day, aIBT, label='ALAN Irradiance at seabed', color='mediumseagreen')
            ax[1].set_ylabel('Irradiance\n ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[1].set_xlim([startdate, enddate])

            # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
            ax[2].set_title(f'Height of water column above datum at {args.datum}m')
            ax[2].plot(dec_day, waterdepth, label='water-column', color='royalblue')
            ax[2].set_ylabel('Water column\n above datum (m)')
            ax2 = ax[2].twinx()
            ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.2, linestyle='dashed', label='datum')
            ax2.legend(prop={"size":8}, loc='upper right')
            ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.2)
            ax2.set_ylabel('Tidal range (m)', color='cadetblue')
            ax[2].set_xlabel("Serial Day")
            ax[2].fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax[2].axes.get_yaxis().set_visible(False)
            ax[2].set_xlim([startdate, enddate])
            ax[2].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])

            plt.tight_layout(pad=0.5, w_pad=0.6, h_pad=0.6)
            if args.save:
                fig.savefig(figurepath)
            plt.show()

        ################## Requires a bit of updating, can plot all plots at once ###################################
            
        if plot_index==4:

            ### PLOT #1 ###
            fig, ax = plt.subplots(4, figsize=(13,7.34))
            # SOLAR INT SEALEVEL
            ax[0].plot(dec_day, Io, color='firebrick')
            ax[0].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[0].set_title('Solar Irradiance at sea level')
            ax[0].set_xlim([startdate, enddate])
            ax0 = ax[0].twinx()
            ax0.set_ylim([(-max(Io)/10), (max(Io)+max(Io)/10)])
            ax0.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax0.axes.get_yaxis().set_visible(False)

            # SOLAR INT SEABED
            ax1 = ax[1].twinx()
            ax1.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax1.axes.get_yaxis().set_visible(False)
            ax[1].plot(dec_day, IBT, color='indianred')
            ax[1].set_title('Solar Irradiance seabed')
            ax[1].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[1].set_xlim([startdate, enddate])
            ax1.set_ylim([(-max(IBT)/10), (max(IBT)+max(IBT)/10)])
            
            # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
            ax2i = ax[2].twinx()
            ax2i.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax2i.axes.get_yaxis().set_visible(False)
            ax2i.set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])
            ax[2].plot(dec_day, waterdepth, color='royalblue')
            ax[2].set_title(f'Height of water column above datum at {args.datum}m')
            ax[2].set_ylabel('water column above datum (m)', color='royalblue')
            ax2 = ax[2].twinx()
            ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.2, linestyle='dashed', label='datum')
            ax2.legend(prop={"size":8}, loc='upper right')
            ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.2)
            ax2.set_ylabel('Tidal range (m)', color='cadetblue')
            ax[2].set_xlim([startdate, enddate])
            ax[2].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])

            # DAY/NIGHT
            ax[3].set_title('Daylight hours')
            ax[3].plot(dec_day, sol, color='gold')
            ax[3].set_ylabel('Daylight')
            ax[3].set_xlabel("Serial Day")
            ax[3].set_xlim([startdate, enddate])
            plt.tight_layout(pad=0.5, w_pad=0.6, h_pad=0.6)
            fig.savefig((figurepath +str('Solar.pdf')))
            
            ### PLOT #2 ###
            fig, ax = plt.subplots(4, figsize=(13,7.34))
            # LUNAR INT SEALEVEL
            ax[0].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$')
            ax[0].set_title('Lunar Irradiance at sea level')
            ax[0].plot(dec_day, I, label='Lunar Irradiance at sea level', color='sienna')
            ax0 = ax[0].twinx()
            ax0.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax0.axes.get_yaxis().set_visible(False)
            ax[0].set_xlim([startdate, enddate])
            ax[0].set_ylim([(-max(I)/10), (max(I)+max(I)/10)])   

            # LUNAR INT SEABED 
            ax1 = ax[1].twinx()
            ax1.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax1.axes.get_yaxis().set_visible(False)
            ax[1].plot(dec_day, lIBT, label='Lunar Irradiance at seabed', color='lightsalmon')
            ax[1].set_title('Lunar Irradiance at seabed ')
            ax[1].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[1].set_xlim([startdate, enddate])
            ax[1].set_ylim([(-max(lIBT)/10), (max(lIBT)+max(lIBT)/10)])

            # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
            ax2i = ax[2].twinx()
            ax2i.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax2i.axes.get_yaxis().set_visible(False)
            ax[2].plot(dec_day, waterdepth, label='water depth', color='royalblue')
            ax[2].set_ylabel('water column above datum (m)', color='royalblue')
            ax2 = ax[2].twinx()
            ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.2, linestyle='dashed', label='datum')
            ax2.legend(prop={"size":8}, loc='upper right')
            ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.2)
            ax2.set_ylabel('Tidal range (m)', color='cadetblue')
            ax[2].set_xlim([startdate, enddate])
            ax[2].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])
            ax[2].set_title(f'Height of water column above datum at {args.datum}m')              

            # LUNAR PHASE
            ax[3].plot(dec_day, phase, label='Lunar phase', color='darkgrey')
            ax[3].set_title('Lunar Phase cycle')
            ax[3].set_ylabel('Normalised Phase')
            ax[3].set_xlabel("Serial Day")
            ax[3].set_xlim([startdate, enddate])
            plt.tight_layout(pad=0.5, w_pad=0.6, h_pad=0.6)
            fig.savefig((figurepath +str('Lunar.pdf')))
            
            ### PLOT #3 ###
            fig, ax = plt.subplots(3, figsize=(13,7.34))
            # ALAN INT SEALEVEL
            ax[0].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[0].plot(dec_day, A, label=f'ALAN at sea level ({condition} sky PAR average = {avgPFFD} PFFD)', color='seagreen')
            ax[0].set_title(f'ALAN Irradiance at sea level PAR average {avgPFFD} PFFD - {condition}')
            ax[0].set_xlim([startdate, enddate])

            # ALAN INT SEABED
            ax[1].set_title('ALAN Irradiance at seabed')
            ax[1].plot(dec_day, aIBT, label='ALAN Irradiance at seabed', color='mediumseagreen')
            ax[1].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')
            ax[1].set_xlim([startdate, enddate])

            # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
            ax[2].set_title(f'Height of water column above datum at {args.datum}m')
            ax[2].plot(dec_day, waterdepth, label='water-column', color='royalblue')
            ax[2].set_ylabel('water column above datum (m)', color='royalblue')
            ax2 = ax[2].twinx()
            ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.2, linestyle='dashed', label='datum')
            ax2.legend(prop={"size":8}, loc='upper right')
            ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.2)
            ax2.set_ylabel('Tidal range (m)', color='cadetblue')
            ax[2].set_xlabel("Serial Day")
            ax2i = ax[2].twinx()
            ax2i.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax2i.axes.get_yaxis().set_visible(False)
            ax[2].set_xlim([startdate, enddate])
            ax[2].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])
            plt.tight_layout(pad=0.5, w_pad=0.6, h_pad=0.6)
            fig.savefig((figurepath +str('ALAN.pdf')))

            ### PLOT #4 ###
            # INTENSITY COMPARISON
            fig, ax = plt.subplots(4, figsize=(13,7.34))
            ax[0].plot(dec_day, Io, label='SolarSS', color='firebrick')
            ax[0].plot(dec_day, I, label='LunarSS', color='sienna', alpha=0.7)
            ax[0].plot(dec_day, A, label=f'ALANSS {condition}', color='seagreen', alpha=0.4)
            ax[0].set_title('Intensity at Sea Surface')
            ax[0].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')

            ax[0].legend(prop={"size":8}, loc='upper right')
            ax0 = ax[0].twinx()
            ax0.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax0.axes.get_yaxis().set_visible(False)
            
            # INTENSITY AT SEABED 
            ax[1].plot(dec_day, IBT, label='SolarSB', color='indianred')
            ax[1].plot(dec_day, lIBT, label='LunarSB', color='lightsalmon', alpha=0.7)
            ax[1].plot(dec_day, aIBT, label=f'ALANSB {condition}', color='mediumseagreen', alpha=0.4)
            ax[1].set_title('Intensity at Sea Bed')
            ax[1].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')

            ax[1].legend(prop={"size":8}, loc='upper right')
            ax1 = ax[1].twinx()
            ax1.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax1.axes.get_yaxis().set_visible(False)
            
            # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
            ax[2].set_title(f'Height of water column above datum at {args.datum}m')
            ax[2].plot(dec_day, waterdepth, label='water-column', color='royalblue')
            ax[2].set_ylabel('water column above datum (m)', color='royalblue')
            ax2 = ax[2].twinx()
            ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.2, linestyle='dashed', label='datum')
            ax2.legend(prop={"size":8}, loc='upper right')
            ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.2)
            ax2.set_ylabel('Tidal range (m)', color='cadetblue')
            ax2i = ax[2].twinx()
            ax2i.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax2i.axes.get_yaxis().set_visible(False)
            ax[2].axes.get_yaxis().set_visible(False)
            ax[2].set_xlim([startdate, enddate])
            ax[2].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])

            # LUNAR PHASE
            ax[3].plot(dec_day, phase, label='Lunar phase', color='darkgrey')
            ax[3].set_title('Lunar Phase cycle')
            ax[3].set_ylabel('Normalised Phase')
            ax[3].set_xlabel("Serial Day")
            ax[3].set_xlim([startdate, enddate])
            plt.tight_layout(pad=0.5, w_pad=0.6, h_pad=0.6)
            fig.savefig((figurepath +str('INT_COMPARISON.pdf')))
            

            ### PLOT #5 ###
            fig, ax = plt.subplots(4, figsize=(13,7.34))
            # JUST LUNAR AND ALAN
            ax[0].plot(dec_day, I, label='LunarSS', color='sienna')
            ax[0].plot(dec_day, A, label=f'ALANSS {condition}', color='seagreen', alpha=0.7)
            ax[0].set_title('Evening Intensity at Sea Surface')
            ax[0].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')

            ax[0].legend(prop={"size":8}, loc='upper right')
            ax0 = ax[0].twinx()
            ax0.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax0.axes.get_yaxis().set_visible(False)

            # EVENING INTENSITY AT SEABED 
            ax[1].plot(dec_day, lIBT, label='LunarSB', color='lightsalmon')
            ax[1].plot(dec_day, aIBT, label=f'ALANSB {condition}', color='mediumseagreen', alpha=0.7)
            ax[1].set_title('Evening Intensity at Sea Bed')
            ax[1].set_ylabel('Irr ($\mu$mol m$^{-2}$ s$^{-1}$)')

            ax[1].legend(prop={"size":8}, loc='upper right')
            ax1 = ax[1].twinx()
            ax1.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax1.axes.get_yaxis().set_visible(False)

            # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
            ax[2].set_title(f'Height of water column above datum at {args.datum}m')
            ax[2].plot(dec_day, waterdepth, label='water-column', color='royalblue')
            ax[2].set_ylabel('water column above datum (m)', color='royalblue')
            ax2 = ax[2].twinx()
            ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.2, linestyle='dashed', label='datum')
            ax2.legend(prop={"size":8}, loc='upper right')
            ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.2)
            ax2.set_ylabel('Tidal range (m)', color='cadetblue')
            ax[2].set_xlabel("Serial Day")
            ax2i = ax[2].twinx()
            ax2i.fill_between(dec_day, night, z, where=[night[i] > z[i] for i in night], facecolor='lightgrey', alpha=0.3)
            ax2i.axes.get_yaxis().set_visible(False)
            ax[2].set_xlim([startdate, enddate])
            ax[2].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])

            # LUNAR PHASE
            ax[3].plot(dec_day, phase, label='Lunar phase', color='darkgrey')
            ax[3].set_title('Lunar Phase cycle')
            ax[3].set_ylabel('Normalised Phase')
            ax[3].set_xlabel("Serial Day")
            plt.tight_layout(pad=0.5, w_pad=0.6, h_pad=0.6)
            fig.savefig((figurepath +str('EVENING_INT_COMPARISON.pdf')))
            plt.show()
        plt.show()
        

# Run script if called from command line.   
if __name__=='__main__':
    
    main()
    
        
        
       
                     





