"""
Make prediction in _predict_tide

        lats : np.ndarray
            Latitudes of the positions to predict.
        times : np.ndarray
            Array of matplotlib datenums (see `matplotlib.dates.num2date').

        coef : utide.utilities.Bunch
            Configuration options for utide.
        amplitudes : np.ndarray
            Amplitude of the relevant constituents shaped [nconst].
        phases : np.ndarray
            Array of the phase of the relevant constituents shaped [nconst].
"""

from utide import reconstruct, ut_constants
from utide.utilities import Bunch
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import pdb
import os
from datetime import datetime, timedelta
pd.options.mode.chained_assignment = None

def get_TidalCoef(geo_location, reftime):

   constituents = ['M2', 'K1', 'S2', 'O1']
   const_idx = np.asarray([ut_constants['const']['name'].tolist().index(i) for i in constituents])
   frq = ut_constants['const']['freq'][const_idx]

   coef = Bunch(name=constituents, mean=0, slope=0)
   coef['aux'] = Bunch(reftime=reftime, lind=const_idx, frq=frq)
   coef['aux']['opt'] = Bunch(twodim=False, nodsatlint=False,
                              nodsatnone=True, gwchlint=False, gwchnone=False, nodiagn=True, 
                              notrend=True, prefilt=[])

   # read in the TPXO file
   # locate the coefficients belonging to the geo_location name (e.g. 'Tokyo')
   TPXOpath = os.getcwd() + "/Required/TPXO/TPXO_ext_tides.csv" # path of data file output
   df_tpxo = pd.read_csv(TPXOpath, delimiter=',', engine='python')  
   print("Extracting data from specific location")
   print(geo_location)

   df_tpxo_ext = df_tpxo.loc[df_tpxo['Name'] == geo_location]
   lats = df_tpxo_ext['Lat'].to_numpy().flatten()
   lons = df_tpxo_ext['Lon'].to_numpy().flatten()

   amplitude = df_tpxo_ext[['M2_amp','K1_amp','S2_amp','O1_amp']].to_numpy().flatten()
   phase = df_tpxo_ext[['M2_phs','K1_phs','S2_phs','O1_phs']].to_numpy().flatten()   
   
   coef['aux']['lat'] = float(lats)  # float
   
   coef['A'] = amplitude
   coef['g'] = phase
   coef['A_ci'] = np.zeros(amplitude.shape)
   coef['g_ci'] = np.zeros(phase.shape)

   return coef
