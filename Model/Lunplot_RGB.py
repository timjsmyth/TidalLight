# Lunar plots
import matplotlib.dates as mdates
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import datetime
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

################# Residuals #################
def LunRes(dec_day, LunSpec, LunI_LS, LunI_LSb, LunI_LSRes, tide_h, waterdepth, sol, datum):
    AA = np.ones(len(dec_day), dtype=int)
    aa = 2300*AA
    bb = -300*AA
    fig, ax = plt.subplots(2, figsize=(13,7.34))
    for i in range(len(LunSpec.columns)):
        if i==0:
            colour = 'firebrick'
            a = 0.7
            label='Red' # '620nm=<Red '
        elif i==1:
            a = 0.7
            colour = 'seagreen'
            label='Green' # '490nm<Green<560nm'
        elif i==2:
            colour = 'steelblue'
            a = 0.7
            label='Blue' # '400nm<Blue<490nm'
        else:
            colour = 'grey'
            a = 0.3
            label=None
               
        ax[0].plot(dec_day, LunI_LSRes.iloc[:,i], label=label + 'Res', color=colour)
        ax[0].plot(dec_day, LunI_LSb.iloc[:,i], label=label + 'Sb', color=colour, linestyle='dashed', alpha=0.5)
        ax[0].plot(dec_day, LunI_LS.iloc[:,i], label=label + 'S', color=colour, linestyle='dotted', alpha=0.5)
        ax[0].set_title('ALAN \n Difference between Surface and Seabed (Datum) (dot - dashed)')
        ax[0].set_ylabel('Irradiance difference\n (\u03bcW m$^{-2}$)')
        ax[0].legend(prop={"size":7}, loc='upper right')

    ax2i = ax[1].twinx()
    ax2i.fill_between(dec_day, aa, bb, where=sol < AA, facecolor='lightgrey', alpha=0.3)
    ax2i.axes.get_yaxis().set_visible(False)
    ax2i.set_ylim([-(max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])
    ax[1].plot(dec_day, waterdepth, color='royalblue')
    ax[1].set_title(f'Height of water column above datum at {datum}m')
    ax[1].set_ylabel('water column above datum (m)', color='royalblue')
    ax2 = ax[1].twinx()
    ax2.axhline(y=2.5, xmin=0, xmax=366, color='cadetblue', alpha = 0.5, linestyle='dashed', label='datum')
    ax2.legend(prop={"size":8}, loc='upper right')
    ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.5)
    ax2.set_ylabel('Tidal range (m)', color='cadetblue')

    ax[1].set_ylim([-(max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])


############################## 3D Normalized Residuals/Spectral Irradiance ##########################
def Lun3d(dec_day, LunSpec, LunI_LSb, LunI_LSRes, LunI_LS, datum):
    fig = plt.figure(figsize=(13,7.34))
    el=6; az=-99; distance=7
    ylab='Wavelength (nm)'; xlab='Serial Day'; zlab='Normalized Irradiance\n (\u03bcW m$^{-2}$)'
    fig.suptitle(f'Normalized Spectral Irradiance at {datum}m below Mean Sea-level')    
    ax = fig.add_subplot(1,2,1, projection='3d')
    el=6; az=-99; distance=7
    ax.view_init(elev=el, azim=az)
    ax.dist=distance
    ax.set_zlabel(zlab)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    DATA = LunI_LSb # plot 1,1,1
    DATA1 = LunI_LSRes # plot 1,1,2
    ax.set_title('Lunar Seabed Spectra: 380-700nm') # DATA title 
    for i in range(len(LunSpec.columns)):
        if i==0:
            colour = 'firebrick'
            a = 0.7
            label='620nm=<Red '
        elif i==1:
            a = 0.7
            colour = 'seagreen'
            label='490nm<Green<560nm'
        elif i==2:
            colour = 'steelblue'
            a = 0.7
            label='400nm<Blue<490nm'
        else:
            colour = 'grey'
            a = 0.3
            label=None
        Axes3D.plot(ax, xs=dec_day, ys=DATA.iloc[:,i]/max(DATA.iloc[:,i]), zs=i, zdir='y', color=colour)
        ax.set_yticks([])

    ax = fig.add_subplot(1,2,2, projection='3d')
    ylab='Wavelength (nm)'; xlab='Serial Day'; zlab='Normalized Irradiance\n (\u03bcW/m$^{-2}$)'
    fig.suptitle(f'Normalized Spectral Irradiance at {datum}m below Mean Sea-level')

    ax.view_init(elev=el, azim=az)
    ax.dist=distance
    ax.set_zlabel(zlab)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    ax.set_title('Lunar through water Residuals Spectra: 380-700nm') # DATA1 title 
    for i in range(len(LunSpec.columns)):
        if i==0:
            colour = 'firebrick'
            a = 0.7
            label='620nm=<Red '
        elif i==1:
            a = 0.7
            colour = 'seagreen'
            label='490nm<Green<560nm'
        elif i==2:
            colour = 'steelblue'
            a = 0.7
            label='400nm<Blue<490nm'
        else:
            colour = 'grey'
            a = 0.3
            label=None
        
        Axes3D.plot(ax, xs=dec_day, ys=DATA1.iloc[:,i]/max(DATA1.iloc[:,i]), zs=i, zdir='y', color=colour)
        ax.set_yticks([])
###################################### Overlay plots ################################
def LunOverlay(dec_day, LunSpec, LunI_LS, LunI_LSb, I, tide_h, waterdepth, sol, lIBT,datum, datum_percentage, phase, location, figurepath):
    fig, ax = plt.subplots(4, figsize=(13.5,7.34))
    fig.suptitle(f'Lunar (Spectral) - {location}\n')
    AA = np.ones(len(dec_day), dtype=int)
    aa = 2300*AA
    bb = -300*AA
    foundmaxLS = []

    # LUNAR BROADBAND SEALEVEL
    ax[0].plot(dec_day, I, color='black')
    #ax0 = ax[0].twinx()

    for i in range(len(LunSpec.columns)):
        if i==0:
            colour = 'firebrick'
            a = 0.7
            label='620nm=<Red '
        elif i==1:
            a = 0.7
            colour = 'seagreen'
            label='490nm<Green<560nm'
        elif i==2:
            colour = 'steelblue'
            a = 0.7
            label='400nm<Blue<490nm'
        else:
            colour = 'grey'
            a = 0.3
            label=None
        maximumLS = max(LunI_LS.iloc[:,i])
        foundmaxLS.append(maximumLS)
        ax[0].plot(dec_day, LunI_LS.iloc[:,i], label=label, color=colour, alpha=a)
    # LUNAR INT SEALEVEL
    #ax[0].set_ylabel('Irr (umol m^-2 s^-1')
    ax[0].set_ylabel('Irradiance\n (\u03bcW m$^{-2}$)')
    # ax[0].set_yscale('symlog')
    ax[0].yaxis.set_major_formatter(ScalarFormatter())
    ax[0].set_title('Surface Lunar Irradiance')

    ax0 = ax[0].twinx()
    ax0.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
    ax0.axes.get_yaxis().set_visible(False)
    ax[0].set_ylim([-(max(foundmaxLS)/10), (max(foundmaxLS)+(max(foundmaxLS)/10))]) 
    ax[0].set_ylim([(-max(I)/10), (max(I)+max(I)/10)])

    # LUNAR INT SEABED 
    foundmaxLSb = []
    ax1 = ax[1].twinx()
    ax1.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
    ax1.axes.get_yaxis().set_visible(False)
    
    ax[1].plot(dec_day, lIBT, color='black')
    for i in range(len(LunSpec.columns)):
        if i==0:
            colour = 'firebrick'
            a = 0.7
            label='620nm=<Red '
        elif i==1:
            a = 0.7
            colour = 'seagreen'
            label='490nm<Green<560nm'
        elif i==2:
            colour = 'steelblue'
            a = 0.7
            label='400nm<Blue<490nm'
        else:
            colour = 'grey'
            a = 0.3
            label=None
        maximumLSb = max(LunI_LSb.iloc[:,i])
        foundmaxLSb.append(maximumLSb)
        ax[1].plot(dec_day, LunI_LSb.iloc[:,i], label=label, color=colour, alpha=a)
    # ax[1].plot(dec_day, lIBT, label='Lunar Irradiance at seabed', color='lightsalmon')
    ax[1].set_title('Intertial Lunar Irradiance')
    #ax[1].set_ylabel('Irr (umol m^-2 s^-1)')
    ax[1].set_ylabel('Irradiance\n (\u03bcW m$^{-2}$)')
    # ax[1].set_yscale('symlog')
    ax[1].yaxis.set_major_formatter(ScalarFormatter())
    # ax[1].set_xlim([182.1, 182.3])
    #ax[1].set_ylim([-(max(foundmaxLSb)/10), (max(foundmaxLSb)+(max(foundmaxLSb)/10))])   
    ax[1].set_ylim([(-max(I)/10), (max(I)+max(I)/10)])

    # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
    #########################
    ax2i = ax[2].twinx()
    ax2i.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
    ax2i.axes.get_yaxis().set_visible(False)
    ax[2].plot(dec_day, waterdepth, label='water depth', color='royalblue')
    ax[2].set_title('Height of water column')
    ax[2].set_ylabel(f'Water column\n above datum at {int(datum_percentage*100)}%', color='black')
    ax2 = ax[2].twinx()
    ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.5, linestyle='dashed', label='datum')
    ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.5)
    ax2.legend(prop={"size":8}, loc='upper right')
    ax2.set_ylabel('Tidal range (m)', color='cadetblue', labelpad=15)
    ax[2].set_ylim([-(max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])

    # LUNAR PHASE
    ax[3].plot(dec_day, phase, label='Lunar phase', color='darkgrey')
    ax[3].set_title('Lunar Phase cycle')
    ax[3].set_ylabel('Normalised Phase')
    ax[3].set_xlabel("Julian Day")
    # ax[3].set_xlim([182.1, 182.3])

    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.5)
    fig.savefig(figurepath +str('Lunar.png'))
