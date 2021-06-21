# Solar plots
import matplotlib.dates as mdates
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import datetime
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import pandas 

######################### Residuals #######################
def ARes(dec_day, ASpec, AI_AS, AI_ASb, AI_ASRes, tide_h, waterdepth, col_names_SS, night, sol, datum):

    fig, ax = plt.subplots(2, figsize=(13,7.34))
    AA = np.ones(len(dec_day), dtype=int)
    aa = 2300*AA
    bb = -300*AA
    for i in range(len(Aspec.columns)):
            if 400<SS_wavenm[i]<500:
                colour = 'steelblue'
            elif 490<SS_wavenm[i]<560:
                colour = 'seagreen'
            elif SS_wavenm[i]>=620:
                colour = 'firebrick'
            else:
                colour = 'black'
            ax[0].plot(dec_day, AI_ASRes.iloc[:,i], color=colour)
            ax[0].plot(dec_day, AI_ASb.iloc[:,i], color=colour, linestyle='dashed', alpha=0.5)
            ax[0].plot(dec_day, AI_AS.iloc[:,i], color=colour, linestyle='dotted', alpha=0.5)


    ax2i = ax[1].twinx()
    ax2i.fill_between(dec_day, aa, bb, where=sol < AA, facecolor='lightgrey', alpha=0.3)
    ax2i.axes.get_yaxis().set_visible(False)
    ax2i.set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])
    ax[1].plot(dec_day, waterdepth, color='royalblue')
    ax[1].set_title(f'Height of water column above datum at {datum}m')
    ax[1].set_ylabel('water column above datum (m)', color='royalblue')
    ax2 = ax[1].twinx()
    ax2.axhline(y=2.5, xmin=0, xmax=366, color='cadetblue', alpha = 0.5, linestyle='dashed', label='datum')
    ax2.legend(prop={"size":8}, loc='upper right')
    ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.5)
    ax2.set_ylabel('Tidal range (m)', color='cadetblue')

    ax[1].set_ylim([(-max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])



############################## 3D Normalised Residuals/ Irradiance ########################### 
def A3d(dec_day, SS_wavenm, AI_ASb, AI_ASRes, AI_AS, col_names_SS, datum, location, skycondition):
  
   

    fig = plt.figure(figsize=(13,7.34))
    fig.suptitle(f'ALAN ({skycondition})- {location}\n')
    ax = fig.add_subplot(1,2,1, projection='3d')
    el=6; az=-99; distance=7
    
    ylab='Wavelength (nm)'; xlab='Julian Day'; zlab='Normalized Irradiance\n (W m$^{-2}$)'
    fig.suptitle(f'Normalized Spectral Irradiance at {datum}m below Mean Sea-level')
    #DATA = SolI_SSb # plot 1,1,1
    #DATA1 = SolI_SSRes # plot 1,1,2

    ax.view_init(elev=el, azim=az)
    ax.dist=distance
    ax.set_zlabel(zlab)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title('Solar Seabed Spectra: 380-700nm') # DATA title
    for i in range(len(col_names_SS)):
        if 400<SS_wavenm[i]<500:
            colour = 'steelblue'
        elif 490<SS_wavenm[i]<560:
            colour = 'seagreen'
        elif SS_wavenm[i]>=620:
            colour = 'firebrick'
        else:
            colour = 'black'
        

        Axes3D.plot(ax, xs=dec_day, ys=AI_ASb.iloc[:,i]/(max(AI_ASb.iloc[:,i])), zs=SS_wavenm[i], zdir='y', color=colour)
    ax = fig.add_subplot(1,2,2, projection='3d')
    el=6; az=-99; distance=7
    
    ylab='Wavelength (nm)'; xlab='Julian Day'; zlab='Normalized Irradiance (W/m$^{-2}$)'
    fig.suptitle(f'Normalized Spectral Irradiance at {datum}m below Mean Sea-level')

    ax.view_init(elev=el, azim=az)
    ax.dist=distance
    ax.set_zlabel(zlab)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title('Solar through water Residuals Spectra: 380-700nm') # DATA1 title
    for i in range(len(col_names_SS)):
        if 400<SS_wavenm[i]<500:
            colour = 'steelblue'
        elif 490<SS_wavenm[i]<560:
            colour = 'seagreen'
        elif SS_wavenm[i]>=620:
            colour = 'firebrick'
        else:
            colour = 'black'
        
        Axes3D.plot(ax, xs=dec_day, ys=AI_ASRes.iloc[:,i]/max(AI_ASRes.iloc[:,i]), zs=SS_wavenm[i], zdir='y', color=colour)

###################### Overlay Plots ####################################   
def AOverlay(dec_day, ASpec, AI_AS, AI_ASb, A, tide_h, waterdepth, night, sol, aIBT, col_names_SS, datum, datum_percentage, location, skycondition, ALAN_TYPE, figurepath):
    fig, ax = plt.subplots(4, figsize=(13.5,7.34))
    fig.suptitle(f'{ALAN_TYPE} - ALAN (Spectral + Broadband) ({skycondition})- {location}\n')
    #fig.suptitle(f'Chlorophyll = {args.chlorophyll}, fCDOM = {args.CDOM}, backscatter = {args.backscatter}')
    # ALAN INT SEALEVEL
    ax[0].plot(dec_day, A, color='black')
    ax0 = ax[0].twinx()
    #ax0.set_ylim([(-max(Io)/10), (max(Io)+max(Io)/10)])

    AA = np.ones(len(dec_day), dtype=int)
    aa = 2300*AA
    bb = -300*AA

    ax0.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
    ax0.axes.get_yaxis().set_visible(False)
    for x in range(len(ASpec.columns)):
        if x==0:
            colour = 'orangered'
            a = 0.7
        elif x==1:
            a = 0.7
            colour = 'limegreen'
        elif x==2:
            colour = 'deepskyblue'
            a = 0.7
        else:
            colour = 'grey'
            a = 0.3
        ax[0].plot(dec_day, AI_AS.iloc[:,x], color=colour, alpha=a)
    # ax[0].set_yscale('symlog')
    ax[0].yaxis.set_major_formatter(ScalarFormatter())
    # ax[0].set_ylabel('Irr (umol m^-2 s^-1)')
    ax[0].set_ylabel('Irradiance\n (\u03bcW m$^{-2}$)')
    ax[0].set_title('ALAN Irradiance at sea level')
    ax0.set_ylim([(-max(A)/10), (max(A)+max(A)/10)])
    

    # INT SEABED
    ax1 = ax[1].twinx()
    ax[1].plot(dec_day, aIBT, color='black')
    # ax[1].plot(dec_day, A, color='black')
    for i in range(len(col_names_SS)):
        if i==0:
            colour = 'orangered'
            a = 0.7
            label='620nm=<Red '
        elif i==1:
            a = 0.7
            colour = 'limegreen'
            label='490nm<Green<560nm'
        elif i==2:
            colour = 'deepskyblue'
            a = 0.7
            label='400nm<Blue<490nm'
        else:
            colour = 'grey'
            a = 0.3
        ax[1].plot(dec_day, AI_ASb.iloc[:,i], label=label, color=colour, alpha=a)
    #ax[1].legend()
    ax[1].set_title('ALAN Irradiance at datum')
    #ax[1].set_ylabel('Irr (umol m^-2 s^-1)')
    ax[1].set_ylabel('Irradiance\n (\u03bcW m$^{-2}$)')
    # ax[1].set_yscale('symlog')
    ax[1].yaxis.set_major_formatter(ScalarFormatter())
    ax1.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_ylim([(-max(aIBT)/10), (max(aIBT)+max(aIBT)/10)])
    
    # TIDAL MODEL WITH REFERENCE DATUM AND WATERCOLUMN HEIGHT
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

    # DAY/NIGHT
    ax[3].set_title('Daylight hours')
    ax[3].plot(dec_day, sol, color='gold')
    ax[3].set_ylabel('Daylight')
    ax[3].set_xlabel("Julian Day")
##            ax[3].set_xlim([startdate, enddate])
##            plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.4)
    fig.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.5)
    fig.savefig(figurepath +str('ALAN.png'))




if __name__ == '__main__':
    print('This script should be run in conjunction with IntensityModel_Modules{variants}.py')
