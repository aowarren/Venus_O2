# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:18:27 2021

@author: sasha
"""

import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt
import glob2 as glob
import seaborn

#%%

#Find all .csv files for model runs from directory where they are saved:
file_list = glob.glob("C:/sasha/Box/Clean Code/results/sensitivity/*.csv")
title = 'data_plots_b_sensitivity.csv'


results = np.zeros((1,28))

#Set arrays for model parameters:
t_arr = np.array([0.5,1.5,3])
fext_arr = np.array([0.1,0.3,0.5,1])
fvolc_arr = np.array([0.1,0.5,0.9,1])
gwat_arr = np.array([100,300,500,700,1000])
CO2_arr = np.array([300,500,1000,2000])
H2O_arr = np.array([0.0001,0.1,0.2,0.5,0.7,1.])

i1 = 0
i2 = 0

for x in file_list:
    i1 = i1+1
    print(i1)

    f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28),header=None)
    last_line = f1.iloc[-1:]
    last_line = last_line.to_numpy()

    b_multi = last_line[0,28]
    tstart = last_line[0,1]
    f_ext = last_line[0,3]
    f_volc = last_line[0,4]
    gwater = last_line[0,8]
    CO2 = last_line[0,9]*1e6
    H2O = last_line[0,10]*1e2
    FMQ = last_line[0,11]

    
    if np.size(last_line) == 0:
        continue
    else:
        results = np.append(results,last_line,axis=0)
    
    if i1 == 10:
        break

results2 = results[1:-1]
data_rem = pd.DataFrame(results2,index=None,columns=None)
data_rem.set_axis(['i1','tstart','time','f_ext','f_volc','b_err','z_erupt','q','gwater',
                    'CO2','H2O','FMQ','Pressure','moles_atm','moles_CO2','moles_CO',
                    'moles_H2','moles_H2O','moles_O2','Rmm_atm','T_surf','ash','lava',
                    'magma','nonthermal','redox','z_melt','moles_CH4'],
                    axis=1,inplace=True)
data_sav = data_rem.drop_duplicates(subset=['tstart','f_ext','f_volc','CO2','H2O','FMQ','gwater'],keep="first",ignore_index=True)

#Save combined data to .csv
data_sav.to_csv(title,header=True,index=False)

#%%

#Import .csv file with all 40Ar model run results (note: MUST have column headers)
# Ar_runs = pd.read_csv('C:/Users/sasha/Box/third_year/venus_earlyhab/volccalc/VolcGases-master/52522_all_Ar.csv')

# #Calculate mean # of successful 40Ar runs for unique combinations of parameters from O2 loss model:
# Ar_summary = Ar_runs.groupby(['tstart', 'f_ext', 'f_volc', 'CO2'])[['success']].agg('mean')

# #Import .csv file with all O2 loss model run results (note: MUST have column headers)
# OG_runs = pd.read_csv('results/'+title)

# output2 = pd.merge(OG_runs,Ar_summary, how = 'left', left_on = ['tstart', 'f_ext', 'f_volc', 'CO2'], right_on = ['tstart','f_ext', 'f_volc', 'CO2'])
# outpu2 = output2.drop_duplicates()

# output2.to_csv('results/'+title+'_Ar.csv',header = True,index=False) 

