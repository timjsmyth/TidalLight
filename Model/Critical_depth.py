import pandas as pd
import numpy as np 
import os
import pdb
import scipy.integrate
import argparse

aparser = argparse.ArgumentParser()
aparser.add_argument("-irr", "--irradiance", action="store", type=float, help="Input ALAN Irradiance in units uW/m^2")
args = aparser.parse_args()
# ALAN_total = args.irradiance
Kds = pd.read_csv("Required/Kd_Eilat_RGB_12Months.csv")
A_zBb = []
zBb = []
KD_Bb = []
critical_depth_Bb = []
crit_depth = np.arange(0,200,0.25)
threshold_intensity = 0.102

columns = ["#1", "#2","#3","#4", "#5","#6","#7", "#8","#9","#10", "#11","#12","#13","#14", "#15","#16","#17", "#18","#19"]
Zc = pd.DataFrame(columns=columns)

loc_name = "/home/awr/OneDrive_PML/Seabed_project-PML/Tamir/Model/Output/OriginalFalchi"
for j in sorted(os.listdir(loc_name)):
    A_zBb = []
    zBb = []
    KD_Bb = []
    critical_depth_Bb = []
    if j[0:6]=="ALAN_i":
        # pdb.set_trace()

        ID = j[6:8]
        if ID[-1] == "_":
            ID = ID[0]
        else:
            pass
        ID = int(ID)
        path = os.path.join(loc_name, j)
        tmp = pd.read_csv(path)
        ALAN_total = float(tmp.iloc[2,15])
        for i in range(len(Kds)):
            # pdb.set_trace()
            KD = Kds.iloc[i,:]
            Kd_Bb = scipy.integrate.simps(KD, dx=1)
            KD_Bb.append(Kd_Bb)

            for zz in range(len(crit_depth)):
                AIBT_z = float(ALAN_total*np.exp(-Kd_Bb*crit_depth[zz]))
                A_zBb.append(AIBT_z)
                zBb.append(crit_depth[zz])
                if AIBT_z < threshold_intensity: # this value is in uW/m^2, 
                    # print("Broadband Intesity", AIBT_z, " Broadband Depth", crit_depth[zz])
                    critical_depth_Bb.append(crit_depth[zz])
                    break
                else:
                    continue
        Zc[f"#{ID}"] = critical_depth_Bb
    else:
        pass
Zc.to_csv("Zc.csv")
