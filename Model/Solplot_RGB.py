# Solar plots
import matplotlib.dates as mdates
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import datetime
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas 
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

def do_clipping(patches, special_y, keep_below=True, ax=None):
    ax = ax or plt.gca()
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    height = ymax - ymin
    if keep_below:
        height = -height
    clip_rect = Rectangle((x, special_y), xmax - xmin, height,
                          transform=ax.transData)
    for p in patches:
        p.set_clip_path(clip_rect)
    
######################### Residuals #######################
def SolRes(dec_day, SolSpec, SolI_SS, SolI_SSb, SolI_SSRes, tide_h, waterdepth, sol, datum, datum_percentage):
    AA = np.ones(len(dec_day), dtype=int)
    aa = 2300*AA
    bb = -300*AA
    fig, ax = plt.subplots(2, figsize=(13,7.34))
    for i in range(len(SolSpec.columns)):
        if i==0:
            colour = 'firebrick'
            a = 0.7
            label='Red '
        elif i==1:
            a = 0.7
            colour = 'seagreen'
            label='Green'
        elif i==2:
            colour = 'steelblue'
            a = 0.7
            label='Blue'
        else:
            colour = 'grey'
            a = 0.3
            label=None

        ax[0].plot(dec_day, SolI_SSRes.iloc[:,i], label=label + 'Res', color=colour)
        ax[0].plot(dec_day, SolI_SSb.iloc[:,i], label=label + 'Sb', color=colour, linestyle='dashed', alpha=0.5)
        ax[0].plot(dec_day, SolI_SS.iloc[:,i], label=label + 'S', color=colour, linestyle='dotted', alpha=0.5)
        ax[0].set_title('Solar \n Difference between Surface and Seabed (Datum) (dot - dashed)')
        ax[0].set_ylabel('Irradiance difference\n (W m$^{-2}$)')
        ax[0].legend(prop={"size":7}, loc='upper right')
          
    
    ax2i = ax[1].twinx()
    ax2i.fill_between(dec_day, aa, bb, where=sol < AA, facecolor='lightgrey', alpha=0.3)
    ax2i.axes.get_yaxis().set_visible(False)
    ax[1].plot(dec_day, waterdepth, color='royalblue')
    ax[1].set_title(f'Height of water column above {datum_percentage*100}% max depth')
    ax[1].set_ylabel('water column above datum (m)', color='royalblue')
    ax2 = ax[1].twinx()
    ax2.axhline(y=datum, xmin=0, xmax=366, color='cadetblue', alpha = 0.2, linestyle='dashed', label='datum')
    ax2.legend(prop={"size":8}, loc='upper right')
    ax2.plot(dec_day, tide_h, label='tide', color='cadetblue', alpha=0.2)
    ax2.set_ylabel('Tidal range (m)', color='cadetblue')
    ax[1].set_ylim([-(max(tide_h)/10), (max(tide_h)+max(tide_h)/10)])



############################## 3D Normalised Residuals/ Irradiance ########################### 
def Sol3d(dec_day, SolSpec, SolI_SSb, SolI_SSRes, SolI_SS, datum, datum_percentage):

    fig = plt.figure(figsize=(13,7.34))
    ax = fig.add_subplot(1,2,1, projection='3d')
    el=6; az=-99; distance=7
    
    ylab='Wavelength (nm)'; xlab='Julian Day'; zlab='Normalized Irradiance W/m$^{-2}$'
    fig.suptitle(f'Normalized Spectral Irradiance at {datum}m below Mean Sea-level and Residuals')


    ax.view_init(elev=el, azim=az)
    ax.dist=distance
    ax.set_zlabel(zlab)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    
    ax.set_title('Solar Seabed Spectra: 380-700nm') # DATA title
    for i in range(len(SolSpec.columns)):
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
        
                #### ISSUE IS 'zs' requires a numeric value, could just have 0, 1, 2 but plot 1 makes this obsolet with only 3 lines
        Axes3D.plot(ax, xs=dec_day, ys=SolI_SSb.iloc[:,i]/(max(SolI_SSb.iloc[:,i])), zs=i, zdir='y', color=colour)
        
        

    ax = fig.add_subplot(1,2,2, projection='3d')
    el=6; az=-99; distance=7
    
    ylab='Wavelength (nm)'; xlab='Julian Day'; zlab='Normalized Irradiance (umol/m^2/s)'
    fig.suptitle(f'Normalized Spectral Irradiance at {datum}m below Mean Sea-level')

    ax.view_init(elev=el, azim=az)
    ax.dist=distance
    ax.set_zlabel(zlab)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    
    ax.set_title('Solar through water Residuals Spectra: 380-700nm') # DATA1 title
    for i in range(len(SolSpec.columns)):
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
               
        Axes3D.plot(ax, xs=dec_day, ys=SolI_SSRes.iloc[:,i]/max(SolI_SSRes.iloc[:,i]), zs=i, zdir='y', color=colour)
        ax.set_yticks([])

###################### Overlay Plots ####################################   
def SolOverlay(dec_day, SolSpec, SolI_SS, SolI_SSb, Io, tide_h, waterdepth, sol, IBT, datum, datum_percentage,location, figurepath):
    fig, ax = plt.subplots(4, figsize=(13.5,7.34))
    # fig.suptitle(f'Chlorophyll = {Chlorophyll}, fCDOM = {fCDOM}, backscatter = {backscatter}')
    fig.suptitle(f'{location}\n')

    # SOLAR INT SEALEVEL
    ax[0].plot(dec_day, Io, color='black')
    ax0 = ax[0].twinx()
    #ax0.set_ylim([(-max(Io)/10), (max(Io)+max(Io)/10)])

    AA = np.ones(len(dec_day), dtype=int)
    aa = 2300*AA
    bb = -300*AA

    ax0.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
    ax0.axes.get_yaxis().set_visible(False)
    for x in range(len(SolSpec.columns)):
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
        ax[0].plot(dec_day, SolI_SS.iloc[:,x], color=colour, alpha=a)
    ax[0].set_yscale('symlog')
    ax[0].yaxis.set_major_formatter(ScalarFormatter())
    # ax[0].set_ylabel('Irr (umol m^-2 s^-1)')
    ax[0].set_ylabel('Irradiance\n (W m$^{-2}$)')
    ax[0].set_title('Surface Solar Irradiance')
##            ax[0].set_xlim([startdate, enddate])
    ax0.set_ylim([(-max(Io)/10), (max(Io)+max(Io)/10)])
    

    # SOLAR INT SEABED
    ax1 = ax[1].twinx()

    ax[1].plot(dec_day, IBT, color='black')
    for i in range(len(SolSpec.columns)):
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
            label=None
        ax[1].plot(dec_day, SolI_SSb.iloc[:,i], label=label, color=colour, alpha=a)
    #ax[1].legend()
    ax[1].set_title('Intertidal Solar Irradiance')
    #ax[1].set_ylabel('Irr (umol m^-2 s^-1)')
    ax[1].set_ylabel('Irradiance\n (W m$^{-2}$)')
    ax[1].set_yscale('symlog')
    ax[1].yaxis.set_major_formatter(ScalarFormatter())
    ax1.fill_between(dec_day, aa, bb, where= sol < AA, facecolor='grey', alpha=0.2)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_ylim([(-max(IBT)/10), (max(IBT)+max(IBT)/10)])
    
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


    # # DAY/NIGHT
    # for i in range(len(sol)):
    #     if sol[i] < 1:
    #         solar_vars = 1
    #     else: 
    #         solar_vars = 0
         

    ax[3].set_title('Daylight hours')
    ax[3].plot(dec_day, sol, color='gold')

    ax[3].set_ylabel('Daylight')
    ax[3].set_xlabel("Julian Day")
#            ax[3].set_xlim([startdate, enddate])
#            plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.4)
    fig.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.5)
    fig.savefig(figurepath +str('Solar.png'))



if __name__ == '__main__':
    print('This script should be run in conjunction with TidalLight_Model.py')
