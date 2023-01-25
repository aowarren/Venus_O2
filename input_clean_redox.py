# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 09:51:22 2020

@author: aowarren
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 08:38:09 2020

@author: aowarren
"""

import numpy as np
import random
import pandas as pd
from datetime import datetime
from modular_functions_clean_redox import time_evol
import time
from matplotlib import pyplot as plt

now = datetime.now()
dt_string = now.strftime("%d-%m-%Y_%H%M")
nsteps = 100000 #number of timesteps

#Atmospheric species
Rmm_CO2 = 0.04401 #kg mol^-1
Rmm_CO = 0.02801 #kg mol^-1
Rmm_CH4 = 0.016 #kg mol^-1
Rmm_H2O = 0.01801528 #kg mol^-1
Rmm_H = 0.001 #kg mol^-1 
Rmm_H2 = 2*Rmm_H
Rmm_O2 = 0.032 #kg mol^-1 
Rmm_O = Rmm_O2/2
Avo = 6.023e23 #Avogadro's number
R = 8.31 #universal gas constant
secingyr = 1e9*365*24*60*60 #seconds in 1 Gyr

#Venus parameters
a_venus = 108.21e6 * 1.e3 #semimajor axis m
g = 8.87 #m s^-2
G = 6.67e-11 #gravitational constant
A_venus = 4.6e14 #m2
R_venus = 6.0518e6 #m
M_venus = 4.867e24 #kg

#Set to 1 for pure metamorphic degassing (no volcanism)
metamorphic_degas = 0


tend = 4.5 #age of Sun Gyr after CAIs
#crustal heatflow, Wm^-2
q = 50e-3
phi0 = random.uniform(0.01,0.2) #uncompressed porosity of Venus crust
#Final atmospheric pressure
P_final = 93.e5
g = 8.87
rho = 2800.
b_err = 1.
T = 1000. + 273. #melt temperature, K
epsilon = 0.125
#thermosphere temperature, K
T_atm = 200.
dzdt_pres = 1e9


rand_start = random.randint(0,28800) # Use if running a large number of runs to select random runs to start from to minimize chance of identical runs occuring simultaneously
# for i0 in range(rand_start,28799):
# rand_start = 10000
    
# for i0 in range(rand_start,rand_start+1): # For individual runs with manually set parameters
for i0 in range(0,28800): # For individual runs with manually set parameters
    print(i0)
    
    # # Load parameters for runs from 'gridsearch_all.csv'
    # param = np.genfromtxt('gridsearch_all_redox.csv',delimiter=',',invalid_raise = False)
    # p_count = np.size(param)
    
    # while p_count < 7: # if p_count file being accessed by another iteration of script, this will tell code to wait 1 second then try to read again to prevent an error
    #     time.sleep(1)
    #     param = np.genfromtxt('gridsearch_all_redox.csv',delimiter=',',invalid_raise = False)
    #     p_count = np.size(param)

    
    # param = np.genfromtxt('gridsearch_all_redox.csv',delimiter=',',invalid_raise = False,usecols=(0,1,2,3,4,5,6))

    # # Extract 
    # FMQ = param[i0,0]
    # tstart = param[i0,1]
    # mH2Otot = param[i0,2] * 1e-2
    # gwater_GEL = param[i0,3]
    # mCO2tot = param[i0,4] * 1e-6
    # frac_volc = param[i0,5]
    # frac_ext = param[i0,6]
    
    FMQ = -1
    tstart = 1.5
    mH2Otot = 0.01 * 1e-2
    gwater_GEL = 300
    mCO2tot = 300 * 1e-6
    frac_volc = 9.0
    frac_ext = 1.0
    
    print('FMQ',FMQ,'tstart',tstart,'q',q,'epsilon',epsilon,'gwater_GEL',gwater_GEL,'CO2',mCO2tot,'H2O',mH2Otot,'FMQ',FMQ,'fext',frac_ext,
          'fvolc',frac_volc,'b_er=',b_err,'Tmelt',T,'T_atm',T_atm)
    #%%
  
    # Run tmospheric evolution model
    t,i1,Pres,moles_atm,Rmm_atm,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,moles_O2,dz,T_surf,ash,lava,magma,nt,redox,z_meltlayer_tot,z_max = time_evol(tstart,nsteps,metamorphic_degas,
                                                                                                 frac_volc,P_final,b_err,dzdt_pres,
                                                                                                 gwater_GEL,T,epsilon,T_atm,mCO2tot,
                                                                                                 mH2Otot,FMQ,frac_ext,q)
    
    
    total_erupt = np.sum(dz)*(t[1]-t[0])/secingyr
    
    if np.sum(t)==0:
        continue

    #%%
    # Save results

    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%H%M")
    gwater_string = str(int(gwater_GEL))
    h2o_string = str(int(1000*mH2Otot))
    f_ext_string = str(int(10*frac_ext))
    f_volc_string = str(int(10*frac_volc))
    CO2_string = str(int(1e6*mCO2tot))
    t_string = str(int(10*tstart))
    FMQ_string = str(int(FMQ))
    
    i1 = int(np.max(i1))
    
    if total_erupt > z_max: # Saves FAILED runs where erupted volume needed to produce all of Venus' modern CO2 exceeds mantle processing upper limit
        title2 = ("C:\\Users\\sasha\\Documents\\results\\nomelt_gwat300_t5\\NMtstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + "_FMQ" + FMQ_string + ".csv")
        print('saving fail')
    else:
        if i1 >= int(nsteps-1000): # Saves runs as complete if they have reached end of nsteps
            title2 = ("C:\\Users\\sasha\\Documents\\results\\nomelt_gwat300_t5\\NMtstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + "_FMQ" + FMQ_string + ".csv")
            print('saving success')
        else: # Saves partial runs to be picked up later
            title2 = ("C:\\Users\\sasha\\Documents\\results\\nomelt_gwat300_t5\\NMtstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + "_FMQ" + FMQ_string + ".csv")
            print('saving partial')    

    # Save data from run:
    data = {'tstart':tstart,
            'time':t[0:i1],
            'f_ext':frac_ext,
            'f_volc':frac_volc,
            'b_err':b_err,
            'z_erupt':dz[0:i1],
            'q':q,
            'gwater':gwater_GEL,
            'CO2':mCO2tot,
            'H2O':mH2Otot,
            'FMQ':FMQ,
            'Pressure':Pres[0:i1],
            'moles_atm':moles_atm[0:i1],
            'moles_CO2':moles_CO2[0:i1],
            'moles_CO':moles_CO[0:i1],
            'moles_H2':moles_H2[0:i1],
            'moles_H2O':moles_H2O[0:i1],
            'moles_O2':moles_O2[0:i1],
            'Rmm_atm':Rmm_atm[0:i1],
            'T_surf':T_surf[0:i1],
            'ash':ash[0:i1],
            'lava':lava[0:i1],
            'magma':magma[0:i1],
            'nonthermal':nt[0:i1],
            'redox':redox[0:i1],
            'z_melt':z_meltlayer_tot[0:i1],
            'moles_CH4':moles_CH4[0:i1]}
        
    df = pd.DataFrame(data)
    # print(df)
    df.to_csv(title2,header = False) 
    
    del df
    del data
    
    #%% ASH VS. LAVA VS. MAGMA VS. SPACE
    
    # i2 = int(i1/2)
    # i2 = i1
    
    # fig1d, ax1d = plt.subplots(figsize=(20,15))
    # plot_ash, = plt.plot(t[0:i2]/(1.e9*365.*60.*60.*24.),(ash[0:i2]),label='Ash',linewidth=4)
    # plot_lava, = plt.plot(t[0:i2]/(1.e9*365.*60.*60.*24.),(lava[0:i2]),label = 'Lava',linewidth=4)
    # plot_magma, = plt.plot(t[0:i2]/(1.e9*365.*60.*60.*24.),(magma[0:i2]),label='Magma',linewidth=4)
    # plot_nt, = plt.plot(t[0:i2]/(1.e9*365.*60.*60.*24.),((nt[0:i2])),label='Nonthermal',linewidth=4)
    # plot_red, = plt.plot(t[0:i2]/(1.e9*365.*60.*60.*24.),((redox[0:i2])),label='Redox',linewidth=4)
    # plot_atm, = plt.plot(t[0:i2]/(1.e9*365.*60.*60.*24.),(moles_O2[0:i2]),label='Atmosphere',linewidth=4)
    # plt.yscale("log")
    # # plt.xlim(0.5,0.501)
    # # plt.xscale("log")
    # plt.legend(handles=[plot_ash,plot_lava,plot_magma,plot_nt,plot_red,plot_atm], fontsize=40)
    # # plt.axhline(y=(1100 + 273), color='k', linestyle='--')
    # ax1d.set_ylabel('Moles O2', fontsize=40)
    # ax1d.set_xlabel('Time (Gyr)', fontsize=40)
    # plt.gca().set_ylim(bottom=1e0,top=1e8)
    # ax1d.tick_params(axis='both', which='major', labelsize=30)
    
    # #%%

    # fig2, ax2 = plt.subplots(figsize=(20,15))
    # plot_colmass, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=4,color='black')
    # plot_colCO2, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=4,color ='indianred')
    # plot_colH2O, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=4,color ='slateblue')
    # plot_colO2, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=4,color ='lightseagreen')
    # plot_colCO, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_CO[0:i1]*Rmm_CO,label='CO',linewidth=4,color ='orange')
    # plot_colH2, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_H2[0:i1]*Rmm_H2,label='H\N{SUBSCRIPT TWO}',linewidth=4,color ='dodgerblue')
    # plot_colCH4, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_CH4[0:i1]*Rmm_CH4,label='CH\N{SUBSCRIPT FOUR}',linewidth=4,color ='yellowgreen')
    # plt.yscale("log")
    # plt.axhline(y=(100.e-6 * moles_atm[i1] * Rmm_H2O),color ='slateblue', linestyle='--',linewidth=4)
    # plt.axhline(y=(50.e-6 * moles_atm[i1] * Rmm_O2), color ='lightseagreen', linestyle='--',linewidth=4)
    # # plt.axhline(y=1, color ='black', linestyle='--',linewidth=2)
    # plt.ylim((1.e-6,1.e8))
    # # plt.xlim(2,2.1)
    # plt.legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2,plot_colCO,plot_colH2,plot_colCH4], fontsize=40)
    # plt.title(str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot)))+ ' wt% H$\mathregular{_{2}}$O' , fontsize=40) 
    # ax2.set_ylabel('Column Mass (kg)', fontsize=40)
    # ax2.set_xlabel('Time (Gyr)', fontsize=40)
    # plt.gca().set_ylim(bottom=1e0)
    # ax2.tick_params(axis='both', which='major', labelsize=30)
    

    # #%% ERUPTED 
    # fig1c, ax1c = plt.subplots(figsize=(20,15))
    # plot_dz = plt.plot(t[np.where((lava+ash)>0)]/(1.e9*365.*60.*60.*24.),dz[np.where((lava+ash)>0)]/1e3,color='indigo',linewidth=4)
    # ax1c.set_ylabel('Eruption rate (km/Gyr)', fontsize=40)
    # ax1c.set_xlabel('Time (Gyr after CAIs)', fontsize=40)
    # ax1c.tick_params(axis='both', which='major', labelsize=30)
    
