# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 14:36:07 2022

@author: sasha
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 09:51:00 2021

@author: sasha
"""


import numpy as np
import math
import pandas as pd
import scipy.stats 
import random
from matplotlib import pyplot as plt
import glob2 as glob
import seaborn

file = 'results/extra/data_plots_wdH2O20_Ar_FULL.csv'
# file = 'results/data_plots_melt_wd.csv_Ar.csv'
# mantle_redox = -2.
#Import .csv file with last line results (generated using gen_data_forplots.py) for all model runs with runaway greenhouse surface melting with 40Ar results added:
results_m_pd = pd.read_csv(file)

file = 'results/extra/data_plots_nomelt_0_Ar_extra.csv'
#Import .csv file with last line results (generated using gen_data_forplots.py) for all baseline model runs with 40Ar results added:
results_nm_pd = pd.read_csv(file)


results_m = results_m_pd.to_numpy()
results_nm = results_nm_pd.to_numpy()

print(np.size(results_m[:,0]))

#%%
# 

measuredH2O = 100.e-6
measuredO2 = 50.e-6
measuredCO = 50.e-6
Pmax = 93.e5
crust_Ar = 100.e3
g = 8.87
Rmm_CO2 = 0.044
rho = 2800 #density of basalt in kg m^{-3}

for i in [1,2]:
    
    if i == 1:
        results = results_m
        melt = 1
        
    if i == 2:
        results = results_nm
        melt = 0

    final_tstep = results[:,0]
    
    FMQ = results[:,11]
    FMQ_arr = np.unique(FMQ)
    
    print(np.size(final_tstep))
    
    # results = results[FMQ==0,:]
    
    tstart = results[:,1]
    t_arr = np.unique(tstart)
    
    
    
    f_ext = results[:,3]
    fext_arr = np.unique(f_ext)
    
    
    f_volc = results[:,4]
    fvolc_arr = np.unique(f_volc)
    
    
    # b_err = results[:,5]
    # b_arr = np.unique(f_ext)
    
    z_erupt = results[:,6]
    
    # q = results[:,7]
    
    gwater = results[:,8]
    gwat_arr = np.unique(gwater)
    
    
    concCO2 = results[:,9]
    CO2_arr = np.unique(concCO2)
    CO2_plotarr = np.array([300,500,1000,2000])
    
    
    concH2O = results[:,10]
    concH2O[concH2O>0.6e-2] = np.round(concH2O[concH2O>0.6e-2],3)
    H2O_arr = np.unique(concH2O)
    
    
    Pres = results[:,12]
    moles_atm = results[:,13]
    moles_CO2 = results[:,14]
    moles_CO = results[:,15]
    moles_H2 = results[:,16]
    moles_H2O = results[:,17]
    moles_O2 = results[:,18]
    Rmm_atm = results[:,19]
    T_surf = results[:,20]
    ash = results[:,21]
    lava = results[:,22]
    magma = results[:,23]
    nonthermal = results[:,24]
    redox = results[:,25]
    z_erupt = ((moles_CO2*Rmm_CO2)-((1-f_volc)*(90e5/g)))/(1000*g*rho*concCO2)
    moles_CH4 = results[:,27]
    success_Ar = results[:,28]
    # results = np.load('dataset_forplots_magma.npy')
    
    t_arr = np.array([0.5,1.5,3])
    
    fext_arr = np.array([0.1,0.3,0.5,1])
    
    fvolc_arr = np.array([0.1,0.5,0.9,1])
    
    gwat_arr = np.array([10,50,100,300,500,700,1000])
    
    CO2_arr = np.array([0.0003,0.0005,0.001,0.002])
    
    H2O_arr = 1e-2*np.array([0.001,0.1,0.2,0.5,0.7,1.,2.])
    
    # H2O_arr = 1e-2*np.array([0.001,0.1,0.5,0.7,1.])
    
    

    
    #success conditions: Pres = 93 bar, O2 < measured, H2O < measured
    
    # if mantle_redox<0:
    #     print('modifying O totals for redox')
    #     if mantle_redox == -1.:
    #         moles_O2 = moles_O2 - 1134869.15085443
    #     if mantle_redox == -2.:
    #         moles_O2 = moles_O2 - 2958902.2110705
    #     if mantle_redox == -3.:
    #         moles_O2 = moles_O2 - 6017222.78276291
    #     if mantle_redox == -4.:
    #         moles_O2 = moles_O2 - 8938942.35908869
    
    end = np.size(FMQ)
    
    success = np.zeros(end)
    
    for i1 in range(0,end):
        frac_o2 = moles_O2[i1]/moles_atm[i1]
        frac_h2o = moles_H2O[i1]/moles_atm[i1]
        frac_co = moles_CO[i1]/moles_atm[i1]
        
        # and frac_co <= measuredCO
        if frac_o2<= measuredO2 and frac_h2o <= measuredH2O and frac_co <= measuredCO and final_tstep[i1] >= 99900.:
            success[i1] = 1.
        else:
            success[i1] = 0.
        if concH2O[i1]==0:
            concH2O[i1] = 0.001e-2
        # if concH2O[i1] >= 0.69e-2 and concH2O[i1] <= 0.71e-2:
        #     concH2O[i1] = 0.7e-2
    print(success)
    
    
    t_totals = 100*np.array([np.sum(success[tstart == 0.5]),np.sum(success[tstart == 1.5]),np.sum(success[tstart == 3])])/(15680)
    g_totals = 100*np.array([np.sum(success[gwater == 10]),np.sum(success[gwater == 50]),np.sum(success[gwater == 100]),np.sum(success[gwater == 300]),np.sum(success[gwater == 500]),np.sum(success[gwater == 700]),np.sum(success[gwater == 1000])])/(6720)
    c_totals = 100*np.array([np.sum(success[concCO2 == 300e-6]),np.sum(success[concCO2 == 500e-6]),np.sum(success[concCO2 == 1000e-6]),np.sum(success[concCO2 == 2000e-6])])/(11760)
    w_totals = 100*np.array([np.sum(success[concH2O == 0.001e-2]),np.sum(success[concH2O == 0.1e-2]),np.sum(success[concH2O == 0.2e-2]),np.sum(success[concH2O == 0.5e-2]),np.sum(success[concH2O == 0.7e-2]),np.sum(success[concH2O == 1.0e-2]),np.sum(success[concH2O == 2.0e-2])])/(6720)
    e_totals = 100*np.array([np.sum(success[f_ext == 0.1]),np.sum(success[f_ext == 0.3]),np.sum(success[f_ext == 0.5]),np.sum(success[f_ext == 1.0])])/(11760)
    v_totals = 100*np.array([np.sum(success[f_volc == 0.1]),np.sum(success[f_volc == 0.5]),np.sum(success[f_volc == 0.9]),np.sum(success[f_volc == 1.0])])/(11760)
    f_totals = 100*np.array([np.sum(success[FMQ == 0]),np.sum(success[FMQ == -1]),np.sum(success[FMQ ==-2]),np.sum(success[FMQ ==-3]),np.sum(success[FMQ ==-4])])/(9408)
    
    # t_totals_Ar = np.array([np.sum(success_Ar[tstart == 0.5]),np.sum(success_Ar[tstart == 1.5]),np.sum(success_Ar[tstart == 3])])/(np.sum(tstart==1.5))
    # g_totals_Ar = np.array([np.sum(success_Ar[gwater == 100]),np.sum(success_Ar[gwater == 300]),np.sum(success_Ar[gwater == 500]),np.sum(success_Ar[gwater == 700]),np.sum(success_Ar[gwater == 1000])])/(np.sum(gwater==100))
    # c_totals_Ar = np.array([np.sum(success_Ar[concCO2 == 300e-6]),np.sum(success_Ar[concCO2 == 500e-6]),np.sum(success_Ar[concCO2 == 1000e-6]),np.sum(success_Ar[concCO2 == 2000e-6])])/(np.sum(concCO2==1000e-6))
    # w_totals_Ar = np.array([np.sum(success_Ar[concH2O == 0.001e-2]),np.sum(success_Ar[concH2O == 0.1e-2]),np.sum(success_Ar[concH2O == 0.2e-2]),np.sum(success_Ar[concH2O == 0.5e-2]),np.sum(success_Ar[concH2O == 0.7e-2]),np.sum(success_Ar[concH2O == 1.0e-2])])/(np.sum(concH2O==1e-2))
    # e_totals_Ar = np.array([np.sum(success_Ar[f_ext == 0.1]),np.sum(success_Ar[f_ext == 0.3]),np.sum(success_Ar[f_ext == 0.5]),np.sum(success_Ar[f_ext == 1.0])])/(np.sum(f_ext==1.0))
    # v_totals_Ar = np.array([np.sum(success_Ar[f_volc == 0.1]),np.sum(success_Ar[f_volc == 0.5]),np.sum(success_Ar[f_volc == 0.9]),np.sum(success_Ar[f_volc == 1.0])])/(np.sum(f_volc==1.0))
    
    t_totals_Ar = t_totals*np.array([np.sum(success_Ar[tstart == 0.5]),np.sum(success_Ar[tstart == 1.5]),np.sum(success_Ar[tstart == 3])])/(15680)
    g_totals_Ar = g_totals*np.array([np.sum(success[gwater == 10]),np.sum(success[gwater == 50]),np.sum(success_Ar[gwater == 100]),np.sum(success_Ar[gwater == 300]),np.sum(success_Ar[gwater == 500]),np.sum(success_Ar[gwater == 700]),np.sum(success_Ar[gwater == 1000])])/(6720)
    c_totals_Ar = c_totals*np.array([np.sum(success_Ar[concCO2 == 300e-6]),np.sum(success_Ar[concCO2 == 500e-6]),np.sum(success_Ar[concCO2 == 1000e-6]),np.sum(success_Ar[concCO2 == 2000e-6])])/(11760)
    w_totals_Ar = w_totals*np.array([np.sum(success_Ar[concH2O == 0.001e-2]),np.sum(success_Ar[concH2O == 0.1e-2]),np.sum(success_Ar[concH2O == 0.2e-2]),np.sum(success_Ar[concH2O == 0.5e-2]),np.sum(success_Ar[concH2O == 0.7e-2]),np.sum(success_Ar[concH2O == 1.0e-2]),np.sum(success_Ar[concH2O == 2.0e-2])])/(6720)
    e_totals_Ar = e_totals*np.array([np.sum(success_Ar[f_ext == 0.1]),np.sum(success_Ar[f_ext == 0.3]),np.sum(success_Ar[f_ext == 0.5]),np.sum(success_Ar[f_ext == 1.0])])/(11760)
    v_totals_Ar = v_totals*np.array([np.sum(success_Ar[f_volc == 0.1]),np.sum(success_Ar[f_volc == 0.5]),np.sum(success_Ar[f_volc == 0.9]),np.sum(success_Ar[f_volc == 1.0])])/(11760)
    f_totals_Ar = f_totals*np.array([np.sum(success_Ar[FMQ == 0]),np.sum(success_Ar[FMQ == -1]),np.sum(success_Ar[FMQ == -2]),np.sum(success_Ar[FMQ == -3]),np.sum(success_Ar[FMQ == -4])])/(9408)
    #%% calc stats
    
    sample_mean_all = []
    
    sample_mean_t05 = []
    sample_mean_t15 = []
    sample_mean_t3 = []
    
    sample_mean_g100 = []
    sample_mean_g300 = []
    sample_mean_g500 = []
    sample_mean_g700 = []
    sample_mean_g1000 = []
    
    sample_mean_c300 = []
    sample_mean_c500 = []
    sample_mean_c1000 = []
    sample_mean_c2000 = []
    
    sample_mean_e01 = []
    sample_mean_e03 = []
    sample_mean_e05 = []
    sample_mean_e10 = []
    
    sample_mean_v01 = []
    sample_mean_v05 = []
    sample_mean_v09 = []
    sample_mean_v10 = []
    
    sample_mean_h0001 = []
    sample_mean_h01 = []
    sample_mean_h02 = []
    sample_mean_h05 = []
    sample_mean_h07 = []
    sample_mean_h1 = []
    sample_mean_h2 = []
    
    sample_mean_f0 = []
    sample_mean_f1 = []
    sample_mean_f2 = []
    sample_mean_f3 = []
    sample_mean_f4 = []
    
    for i in range(100000):
      y = random.choices(success.tolist())
      avg = np.mean(y)
      sample_mean_all.append(avg)
      
      y = random.choices(success[tstart == 0.5].tolist())
      avg = np.mean(y)
      sample_mean_t05.append(avg)
      
      y = random.choices(success[tstart == 1.5].tolist())
      avg = np.mean(y)
      sample_mean_t15.append(avg)
      
      y = random.choices(success[tstart == 3.0].tolist())
      avg = np.mean(y)
      sample_mean_t3.append(avg)
      
      y = random.choices(success[gwater == 100].tolist())
      avg = np.mean(y)
      sample_mean_g100.append(avg)
      
      y = random.choices(success[gwater == 300].tolist())
      avg = np.mean(y)
      sample_mean_g300.append(avg)
      
      y = random.choices(success[gwater == 500].tolist())
      avg = np.mean(y)
      sample_mean_g500.append(avg)
      
      y = random.choices(success[gwater == 700].tolist())
      avg = np.mean(y)
      sample_mean_g700.append(avg)
      
      y = random.choices(success[gwater == 1000].tolist())
      avg = np.mean(y)
      sample_mean_g1000.append(avg)
      
      y = random.choices(success[concCO2 == 300e-6].tolist())
      avg = np.mean(y)
      sample_mean_c300.append(avg)
      
      y = random.choices(success[concCO2 == 500e-6].tolist())
      avg = np.mean(y)
      sample_mean_c500.append(avg)
      
      y = random.choices(success[concCO2 == 1000e-6].tolist())
      avg = np.mean(y)
      sample_mean_c1000.append(avg)
      
      y = random.choices(success[concCO2 == 2000e-6].tolist())
      avg = np.mean(y)
      sample_mean_c2000.append(avg)
      
      y = random.choices(success[f_ext == 0.1].tolist())
      avg = np.mean(y)
      sample_mean_e01.append(avg)
      
      y = random.choices(success[f_ext == 0.3].tolist())
      avg = np.mean(y)
      sample_mean_e03.append(avg)
      
      y = random.choices(success[f_ext == 0.5].tolist())
      avg = np.mean(y)
      sample_mean_e05.append(avg)
      
      y = random.choices(success[f_ext == 1.0].tolist())
      avg = np.mean(y)
      sample_mean_e10.append(avg)
      
      y = random.choices(success[f_volc == 0.1].tolist())
      avg = np.mean(y)
      sample_mean_v01.append(avg)
      
      y = random.choices(success[f_volc == 0.5].tolist())
      avg = np.mean(y)
      sample_mean_v05.append(avg)
      
      y = random.choices(success[f_volc == 0.9].tolist())
      avg = np.mean(y)
      sample_mean_v09.append(avg)
      
      y = random.choices(success[f_volc == 1.0].tolist())
      avg = np.mean(y)
      sample_mean_v10.append(avg)
      
      y = random.choices(success[concH2O == 0.001e-2].tolist())
      avg = np.mean(y)
      sample_mean_h0001.append(avg)
      
      y = random.choices(success[concH2O == 0.1e-2].tolist())
      avg = np.mean(y)
      sample_mean_h01.append(avg)
      
      y = random.choices(success[concH2O == 0.2e-2].tolist())
      avg = np.mean(y)
      sample_mean_h02.append(avg)
      
      y = random.choices(success[concH2O == 0.5e-2].tolist())
      avg = np.mean(y)
      sample_mean_h05.append(avg)
      
      y = random.choices(success[concH2O == 0.7e-2].tolist())
      avg = np.mean(y)
      sample_mean_h07.append(avg)
      
      y = random.choices(success[concH2O == 1e-2].tolist())
      avg = np.mean(y)
      sample_mean_h1.append(avg)
      
      y = random.choices(success[concH2O == 2e-2].tolist())
      avg = np.mean(y)
      sample_mean_h1.append(avg)
      
      y = random.choices(success[FMQ == 0].tolist())
      avg = np.mean(y)
      sample_mean_f0.append(avg)
      
      y = random.choices(success[FMQ == -1].tolist())
      avg = np.mean(y)
      sample_mean_f1.append(avg)
      
      y = random.choices(success[FMQ == -2].tolist())
      avg = np.mean(y)
      sample_mean_f2.append(avg)
      
      y = random.choices(success[FMQ == -3].tolist())
      avg = np.mean(y)
      sample_mean_f3.append(avg)
      
      y = random.choices(success[FMQ == -4].tolist())
      avg = np.mean(y)
      sample_mean_f4.append(avg)
     
    mean_succ = (np.mean(sample_mean_all))
    
    conf = scipy.stats.norm.interval(alpha=0.95, loc=mean_succ, scale=scipy.stats.sem(sample_mean_all))
    print('mean',mean_succ,'95% confidence intervals',conf)
    
    mean_t05  = (np.mean(sample_mean_t05))
    mean_t15  = (np.mean(sample_mean_t15))
    mean_t3  = (np.mean(sample_mean_t3))
    mean_g100  = (np.mean(sample_mean_g100))
    mean_g300  = (np.mean(sample_mean_g300))
    mean_g500  = (np.mean(sample_mean_g500))
    mean_g700  = (np.mean(sample_mean_g700))
    mean_g1000  = (np.mean(sample_mean_g1000))
    mean_c300  = (np.mean(sample_mean_c300))
    mean_c500  = (np.mean(sample_mean_c500))
    mean_c1000  = (np.mean(sample_mean_c1000))
    mean_c2000  = (np.mean(sample_mean_c2000))
    mean_e01  = (np.mean(sample_mean_e01))
    mean_e03  = (np.mean(sample_mean_e03))
    mean_e05  = (np.mean(sample_mean_e05))
    mean_e10  = (np.mean(sample_mean_e10))
    mean_v01  = (np.mean(sample_mean_v01))
    mean_v05  = (np.mean(sample_mean_v05))
    mean_v09  = (np.mean(sample_mean_v09))
    mean_v10  = (np.mean(sample_mean_v10))
    mean_h0001  = (np.mean(sample_mean_h0001))
    mean_h01  = (np.mean(sample_mean_h01))
    mean_h02  = (np.mean(sample_mean_h02))
    mean_h05  = (np.mean(sample_mean_h05))
    mean_h07  = (np.mean(sample_mean_h07))
    mean_h1  = (np.mean(sample_mean_h1))
    mean_h2  = (np.mean(sample_mean_h2))
    mean_f0  = (np.mean(sample_mean_f0))
    mean_f1  = (np.mean(sample_mean_f1))
    mean_f2  = (np.mean(sample_mean_f2))
    mean_f3  = (np.mean(sample_mean_f3))
    mean_f4  = (np.mean(sample_mean_f4))
    
    mean_t = np.array([ mean_t05, mean_t15, mean_t3])
    mean_g = np.array([ mean_g100, mean_g300, mean_g500, mean_g700, mean_g1000])
    mean_c = np.array([ mean_c300,mean_c500,mean_c1000,mean_c2000])
    mean_h = np.array([ mean_h0001,mean_h01,mean_h02,mean_h05,mean_h07,mean_h1,mean_h2])
    mean_e = np.array([ mean_e01,mean_e03,mean_e05,mean_e10])
    mean_v = np.array([ mean_v01,mean_v05,mean_v09,mean_v10])
    mean_f = np.array([ mean_f0,mean_f1,mean_f2,mean_f3,mean_f4])
    
    
    conf1 = scipy.stats.norm.interval(alpha=0.95, loc=mean_t05, scale=scipy.stats.sem(sample_mean_t05))
    conf2 = scipy.stats.norm.interval(alpha=0.95, loc=mean_t15, scale=scipy.stats.sem(sample_mean_t15))
    conf3 = scipy.stats.norm.interval(alpha=0.95, loc=mean_t3, scale=scipy.stats.sem(sample_mean_t3))
    
    conf_t = np.array([(mean_t05 - conf1[0]) + (conf1[1] - mean_t05),
                       (mean_t15 - conf2[0]) + (conf2[1] - mean_t15),
                       (mean_t3 - conf3[0]) + (conf3[1] - mean_t3) ])

    conf1 = scipy.stats.norm.interval(alpha=0.95, loc=mean_g100, scale=scipy.stats.sem(sample_mean_g100))
    conf2 = scipy.stats.norm.interval(alpha=0.95, loc=mean_g300, scale=scipy.stats.sem(sample_mean_g300))
    conf3 = scipy.stats.norm.interval(alpha=0.95, loc=mean_g500, scale=scipy.stats.sem(sample_mean_g500))
    conf4 = scipy.stats.norm.interval(alpha=0.95, loc=mean_g700, scale=scipy.stats.sem(sample_mean_g700))
    conf5 = scipy.stats.norm.interval(alpha=0.95, loc=mean_g1000, scale=scipy.stats.sem(sample_mean_g1000))
    
    conf_g = np.array([(mean_g100 - conf1[0]) + (conf1[1] - mean_g100),
                       (mean_g300 - conf2[0]) + (conf2[1] - mean_g300),
                       (mean_g500 - conf3[0]) + (conf3[1] - mean_g500),
                       (mean_g700 - conf4[0]) + (conf4[1] - mean_g700),
                       (mean_g1000 - conf5[0]) + (conf5[1] - mean_g1000)])

    conf1 = scipy.stats.norm.interval(alpha=0.95, loc=mean_h0001, scale=scipy.stats.sem(sample_mean_h0001))
    conf2 = scipy.stats.norm.interval(alpha=0.95, loc=mean_h01, scale=scipy.stats.sem(sample_mean_h01))
    conf3 = scipy.stats.norm.interval(alpha=0.95, loc=mean_h02, scale=scipy.stats.sem(sample_mean_h02))
    conf4 = scipy.stats.norm.interval(alpha=0.95, loc=mean_h05, scale=scipy.stats.sem(sample_mean_h05))
    conf5 = scipy.stats.norm.interval(alpha=0.95, loc=mean_h07, scale=scipy.stats.sem(sample_mean_h07))
    conf6 = scipy.stats.norm.interval(alpha=0.95, loc=mean_h1, scale=scipy.stats.sem(sample_mean_h1))
    conf7 = scipy.stats.norm.interval(alpha=0.95, loc=mean_h2, scale=scipy.stats.sem(sample_mean_h2))
    
    
    conf_h = np.array([(mean_h0001 - conf1[0]) + (conf1[1] - mean_h0001),
                       (mean_h01 - conf2[0]) + (conf2[1] - mean_h01),
                       (mean_h02 - conf3[0]) + (conf3[1] - mean_h02),
                       (mean_h05 - conf4[0]) + (conf4[1] - mean_h05),
                       (mean_h07 - conf5[0]) + (conf5[1] - mean_h07),
                       (mean_h1 - conf6[0]) + (conf6[1] - mean_h1),
                       (mean_h2 - conf7[0]) + (conf7[1] - mean_h2)])
    
    conf1 = scipy.stats.norm.interval(alpha=0.95, loc=mean_c300, scale=scipy.stats.sem(sample_mean_c300))
    conf2 = scipy.stats.norm.interval(alpha=0.95, loc=mean_c500, scale=scipy.stats.sem(sample_mean_c500))
    conf3 = scipy.stats.norm.interval(alpha=0.95, loc=mean_c1000, scale=scipy.stats.sem(sample_mean_c1000))
    conf4 = scipy.stats.norm.interval(alpha=0.95, loc=mean_c2000, scale=scipy.stats.sem(sample_mean_c2000))
    
    conf_c = np.array([(mean_c300 - conf1[0]) + (conf1[1] - mean_c300),
                       (mean_c500 - conf2[0]) + (conf2[1] - mean_c500),
                       (mean_c1000 - conf3[0]) + (conf3[1] - mean_c1000),
                       (mean_c2000 - conf4[0]) + (conf4[1] - mean_c2000)])
    
    conf1 = scipy.stats.norm.interval(alpha=0.95, loc=mean_e01, scale=scipy.stats.sem(sample_mean_e01))
    conf2 = scipy.stats.norm.interval(alpha=0.95, loc=mean_e03, scale=scipy.stats.sem(sample_mean_e03))
    conf3 = scipy.stats.norm.interval(alpha=0.95, loc=mean_e05, scale=scipy.stats.sem(sample_mean_e05))
    conf4 = scipy.stats.norm.interval(alpha=0.95, loc=mean_e10, scale=scipy.stats.sem(sample_mean_e10))
    
    conf_e = np.array([(mean_e01 - conf1[0]) + (conf1[1] - mean_e01),
                       (mean_e03 - conf2[0]) + (conf2[1] - mean_e03),
                       (mean_e05 - conf3[0]) + (conf3[1] - mean_e05),
                       (mean_e10 - conf4[0]) + (conf4[1] - mean_e10)])
    
    conf1 = scipy.stats.norm.interval(alpha=0.95, loc=mean_v01, scale=scipy.stats.sem(sample_mean_v01))
    conf2 = scipy.stats.norm.interval(alpha=0.95, loc=mean_v05, scale=scipy.stats.sem(sample_mean_v05))
    conf3 = scipy.stats.norm.interval(alpha=0.95, loc=mean_v09, scale=scipy.stats.sem(sample_mean_v09))
    conf4 = scipy.stats.norm.interval(alpha=0.95, loc=mean_v10, scale=scipy.stats.sem(sample_mean_v10))
    
    conf_v = np.array([(mean_v01 - conf1[0]) + (conf1[1] - mean_v01),
                       (mean_v05 - conf2[0]) + (conf2[1] - mean_v05),
                       (mean_v09 - conf3[0]) + (conf3[1] - mean_v09),
                       (mean_v10 - conf4[0]) + (conf4[1] - mean_v10)])
    
    
    conf1 = scipy.stats.norm.interval(alpha=0.95, loc=mean_f0, scale=scipy.stats.sem(sample_mean_f0))
    conf2 = scipy.stats.norm.interval(alpha=0.95, loc=mean_f1, scale=scipy.stats.sem(sample_mean_f1))
    conf3 = scipy.stats.norm.interval(alpha=0.95, loc=mean_f2, scale=scipy.stats.sem(sample_mean_f2))
    conf4 = scipy.stats.norm.interval(alpha=0.95, loc=mean_f3, scale=scipy.stats.sem(sample_mean_f3))
    conf5 = scipy.stats.norm.interval(alpha=0.95, loc=mean_f4, scale=scipy.stats.sem(sample_mean_f4))
    
    
    conf_f = np.array([(mean_f0 - conf1[0]) + (conf1[1] - mean_f0),
                       (mean_f1 - conf2[0]) + (conf2[1] - mean_f1),
                       (mean_f2 - conf3[0]) + (conf3[1] - mean_f2),
                       (mean_f3 - conf4[0]) + (conf4[1] - mean_f3),
                       (mean_f4 - conf5[0]) + (conf5[1] - mean_f4)])
    #%%
    if melt==1:
        title = 'e_wd_statistics_melt_CO_FMQ0.npz'

        
        # np.savez('melt_totals.npz',g_totals_m=g_totals,e_totals_m=e_totals,t_totals_m=t_totals,c_totals_m=c_totals,w_totals_m=w_totals,v_totals_m=v_totals,mean_succ_m = mean_succ,conf_m = conf)
        np.savez(title,t_m = t_totals,
                     g_m = g_totals,
                     c_m = c_totals,
                     h_m = w_totals,
                     e_m = e_totals,
                     v_m = v_totals,
                     f_m = f_totals,
                     all_m = mean_succ,
                     t_m_Ar =t_totals_Ar,
                     g_m_Ar =g_totals_Ar,
                     c_m_Ar =c_totals_Ar,
                     h_m_Ar =w_totals_Ar,
                     e_m_Ar =e_totals_Ar,
                     v_m_Ar =v_totals_Ar,
                     f_m_Ar =f_totals_Ar,
                     ct_m = conf_t,
                     cg_m = conf_g,
                     cc_m = conf_c,
                     ch_m = conf_h,
                     ce_m = conf_e,
                     cv_m = conf_v,
                     cf_m = conf_f,
                     call_m = conf)
        # melt_stats = df
        # df.to_csv('statistics_melt.csv',header = False) 
    elif melt==0:
        
        title = 'e_statistics_nomelt_CO_FMQ0.npz'
        
        # np.savez('nomelt_totals.npz',g_totals_nm=g_totals,e_totals_nm=e_totals,t_totals_nm=t_totals,c_totals_nm=c_totals,w_totals_nm=w_totals,v_totals_nm=v_totals,mean_succ_nm = mean_succ,conf_nm = conf)
        # df.to_csv('statistics_nomelt.csv',header = False) 
    
        np.savez(title,t_nm = t_totals,
                     g_nm = g_totals,
                     c_nm = c_totals,
                     h_nm = w_totals,
                     e_nm = e_totals,
                     v_nm = v_totals,
                     f_nm = f_totals,
                     all_nm = mean_succ,
                     t_nm_Ar =t_totals_Ar,
                     g_nm_Ar =g_totals_Ar,
                     c_nm_Ar =c_totals_Ar,
                     h_nm_Ar =w_totals_Ar,
                     e_nm_Ar =e_totals_Ar,
                     v_nm_Ar =v_totals_Ar,
                     f_nm_Ar =f_totals_Ar,
                     ct_nm = conf_t,
                     cg_nm = conf_g,
                     cc_nm = conf_c,
                     ch_nm = conf_h,
                     ce_nm = conf_e,
                     cv_nm = conf_v,
                     cf_nm = conf_f,
                     call_nm = conf)
    


