# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 11:52:58 2022

@author: sasha
"""

import numpy as np
import pandas as pd
from datetime import datetime
from scipy.interpolate import Rbf
from scipy.interpolate import interpn
import random
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import rc
import glob2 as glob
import seaborn

#%%
#Atmospheric species
Rmm_CO2 = 0.04401 #kg mol^-1
Rmm_CO = 0.02801 #kg mol^-1
Rmm_CH4 = 0.016
Rmm_H2O = 0.01801528 #kg mol^-1
Rmm_H = 0.001 #kg mol^-1 
Rmm_H2 = 2*Rmm_H
Rmm_O2 = 0.032 #kg mol^-1 
Rmm_O = Rmm_O2/2
Avo = 6.023e23 #Avogadro's number
R = 8.31 #universal gas constant 
rho = 2800 #density of basalt kg m^-3

#Venus parameters
a_venus = 108.21e6 * 1.e3 #semimajor axis m
g = 8.87 #m s^-2
G = 6.67e-11 #gravitational constant
A_venus = 4.6e14 #m2
R_venus = 6.0518e6 #m
M_venus = 4.867e24 #kg

secingyr = 1.e9*365.*60.*60.*24.

t_arr = np.array([0.5,1.5,3])
fext_arr = np.array([0.1,0.3,0.5,1])
fvolc_arr = np.array([0.1,0.5,0.9,1])
gwat_arr = np.array([10,50,100,300,500,700,1000])
CO2_arr = np.array([0.0003,0.0005,0.001,0.002])
H2O_arr = 1e-2*np.array([0.001,0.1,0.2,0.5,0.7,1.,2.])
H2O_arr[4] = 7e-3
FMQ_arr = np.array([0,-1,-2,-3,-4])

moles_atm_old = 0
Rmm_atm_old = 0
moles_CO2_old = 0
moles_H2O_old = 0
moles_O2_old = 0
#%%




plot_number = 4 #For Fig 2: 1, Fig 3: 4. For SI Appendix Fig S1: 2, Fig S2: 3, Fig S8: 5, Fig S9: 6.

melt = 0 #can set to 1 if want to plot runs with melting

counter = 0
i2 = 0

#%%
if plot_number == 1:
    
    if melt == 0:
        file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\nomelt_gwat100_t5\\*.csv")
        # file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\nomelt_gwat300_t5\\*.csv")
    elif melt == 1:
        file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\melt_wd_gwat500_5\\*.csv")
    
    fig1, ax1 = plt.subplots(2,figsize=(45,60))
    fig1.tight_layout(pad=25.0)
    
    ax1[0].set_ylim((1.e1,3.e6))
    ax1[1].set_ylim((1.e1,3.e6))
    
    ax1[0].set_yscale("log")
    ax1[1].set_yscale("log")
    ax1[0].invert_xaxis()
    ax1[1].invert_xaxis()
    
    
    ax1[0].set_ylabel('Column Mass (kg)', fontsize=100,weight='bold')
    ax1[0].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    ax1[1].set_ylabel('Column Mass (kg)', fontsize=100,weight='bold')
    ax1[1].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    
    # plt.gca().set_ylim(bottom=1e0)
    ax1[0].tick_params(axis='both', which='major', labelsize=100)
    ax1[1].tick_params(axis='both', which='major', labelsize=100)
    # plt.ylim((1.e-6,1.e8))
    # plt.yscale("log")
    
    # ax1[0].set_title('H\N{SUBSCRIPT TWO}O', fontsize=120)
    # ax1[1].set_title('O\N{SUBSCRIPT TWO}', fontsize=120)
    
    #%%
    
    for x in file_list:
        counter = counter + 1
        # print('counter',counter)
        
        # if counter == 3:
        #     break
        # if i2 == 3:
        #     break
        
        f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27),header=None)
        
        results = f1.to_numpy()
        
        i = results[:,0]
        i = i.astype(int)
        i1 = i[-1]
        
        final_tstep = results[:,0]
        
        tstart = results[0,1]
        # t_arr = np.unique(tstart)
        
        t = results[:,2]
        
        f_ext = results[0,3]
        # fext_arr = np.unique(f_ext)
        
        f_volc = results[0,4]
        # fvolc_arr = np.unique(f_volc)
        
        z_erupt = results[:,6]
        
        gwater = results[0,8]
        gwat_arr = np.unique(gwater)
        
        mCO2tot = results[0,9]
        CO2_arr = np.unique(mCO2tot)
        CO2_plotarr = np.array([300,500,1000,2000])
        
        mH2Otot = np.round(results[0,10],3)
        # mH2Otot[mH2Otot>0.6e-2] = np.round(mH2Otot[mH2Otot>0.6e-2],3)
        
        FMQ = results[0,11]
        
        
        
        
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
        z_erupt = ((moles_CO2*Rmm_CO2)-((1-f_volc)*(90e5/g)))/(1000*g*rho*mCO2tot)
        moles_CH4 = results[:,27]
        # success_Ar = results[:,28]
        
        
        capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
        
        
        # z_erupt = (((moles_CO2 - np.append(0,moles_CO2[0:-1]))*Rmm_CO2)/(g*rho*mCO2tot))/((t - np.append(0,t[0:-1]))/secingyr)
    
        # z_erupt[T_surf>0]==0
        
        if all(t_arr != tstart):
            continue
        if all(CO2_arr != mCO2tot):
            continue
        if all(H2O_arr != mH2Otot):
            continue
        if all(gwat_arr != gwater):
            continue
        if all(fext_arr != f_ext):
            continue
        if all(fvolc_arr != f_volc):
            continue 
        
        if melt == 1:
            gwater_string = str(int(gwater))
            h2o_string = str(int(1000*mH2Otot))
            f_ext_string = str(int(10*f_ext))
            f_volc_string = str(int(10*f_volc))
            CO2_string = str(int(1e6*mCO2tot))
            t_string = str(int(10*tstart))
            FMQ_string = str(int(FMQ))
            
            title = "C:\\Users\\sasha\\Box\\Clean Code\\results\\reduced_melt\\withdiss_t" + t_string + "_gwater" + gwater_string + "_f_volc" + f_volc_string + ".csv"
            
            melt_res = pd.read_csv(title)
        
            diss_dat = melt_res.to_numpy()
            
            time_d = diss_dat[:,4]
            Pres_d = diss_dat[:,5]
            moles_atm_d = diss_dat[:,6]
            moles_CO2_d = diss_dat[:,7]
            moles_H2O_d = diss_dat[:,8]
            moles_O2_d = diss_dat[:,9]
            nonthermal_d = diss_dat[:,10]
            magma_d = diss_dat[:,11]
            Rmm_atm_d = diss_dat[:,12]
            
            blanks = np.zeros(np.size(time_d))
            
            Pres = np.append(Pres_d,Pres)
            moles_atm = np.append(moles_atm_d,moles_atm)
            moles_CO2 = np.append(moles_CO2_d,moles_CO2)
            moles_CO = np.append(blanks,moles_CO)
            moles_H2 = np.append(blanks,moles_H2)
            moles_H2O = np.append(moles_H2O_d,moles_H2O)
            moles_O2 = np.append(moles_O2_d,moles_O2)
            Rmm_atm = np.append(Rmm_atm_d,Rmm_atm)
            T_surf = np.append(blanks,T_surf)
            ash = np.append(blanks,ash)
            lava = np.append(blanks,lava)
            magma = np.append(magma_d,magma)
            nonthermal = np.append(nonthermal_d,nonthermal)
            redox = np.append(blanks,redox)
            z_erupt = np.append(blanks,z_erupt)
            moles_CH4 = np.append(blanks,moles_CH4)
        # success_Ar = results[:,28]
        
        
            capN = capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
        
        
        
        i2 = i2 + 1
        print('valid, i2', i2,' gwat', gwater, ' tstart', tstart)
        
        # plot_colmass, = plt.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=10,color='lightgrey',linestyle=':')
        
        # plot_colCO2, = plt.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=10,color ='bisque')
        
        ax1[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=15,color = 'slateblue',alpha=0.01)
        
        ax1[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=15,color ='lightseagreen',alpha =0.01)
        
        
        # plt.axhline(y=(100.e-6 * moles_atm[i1] * Rmm_H2O),color ='lightsteelblue', linestyle='--',linewidth=4)
        # plt.axhline(y=(50.e-6 * moles_atm[i1] * Rmm_O2), color ='paleturquoise', linestyle='--',linewidth=4)
        # plt.axhline(y=1, color ='lightgrey', linestyle='--',linewidth=2)
        
        # plt.xlim(2,2.1)
        # plt.legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=40)
        # plt.legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=50)
        # plt.title(str(float('%g' % (gwater))) + ' m GEL, '+ str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) + ' wt% H$\mathregular{_{2}}$O, f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))), fontsize=40) 
        # ax1.set_ylabel('Column Mass (kg)', fontsize=50,weight='bold')
        # ax1.set_xlabel('Time before present (Gyr)', fontsize=50,weight='bold')
        # plt.gca().set_ylim(bottom=1e0)
        # ax1.tick_params(axis='both', which='major', labelsize=40)
        
        # if counter > 1:
            
        #     ax1.fill_between(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),
        #                                      moles_atm[0:i1]*Rmm_atm[0:i1],
        #                                      moles_atm_old[0:i1]*Rmm_atm_old[0:i1],
        #                                      color='lightgrey')
        #     ax1.fill_between(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),
        #                                     moles_CO2[0:i1]*Rmm_CO2,
        #                                     moles_CO2_old[0:i1]*Rmm_CO2,
        #                                     color ='bisque')
        #     ax1.fill_between(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),
        #                                     moles_H2O[0:i1]*Rmm_H2O,
        #                                     moles_H2O_old[0:i1]*Rmm_H2O,
        #                                     color ='lightsteelblue')
        #     ax1.fill_between(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),
        #                                    moles_O2[0:i1]*Rmm_O2,
        #                                    moles_O2_old[0:i1]*Rmm_O2,
        #                                    color ='paleturquoise')
        
        # moles_atm_old = moles_atm
        # Rmm_atm_old = Rmm_atm
        # moles_CO2_old = moles_CO2
        # moles_H2O_old = moles_H2O
        # moles_O2_old = moles_O2
        
    #%%
    # fig1, ax1 = plt.subplots(2,figsize=(45,60))
    # fig1.tight_layout(pad=25.0)
    
    # ax1[0].set_ylim((1.e1,3.e6))
    # ax1[1].set_ylim((1.e1,3.e6))
    
    # ax1[0].set_yscale("log")
    # ax1[1].set_yscale("log")
    # ax1[0].invert_xaxis()
    # ax1[1].invert_xaxis()
    
    
    # ax1[0].set_ylabel('Column Mass (kg)', fontsize=100,weight='bold')
    # ax1[0].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    # ax1[1].set_ylabel('Column Mass (kg)', fontsize=100,weight='bold')
    # ax1[1].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    
    # # plt.gca().set_ylim(bottom=1e0)
    # ax1[0].tick_params(axis='both', which='major', labelsize=100)
    # ax1[1].tick_params(axis='both', which='major', labelsize=100)
    #Highlight specific plot
    # f1 = pd.read_csv("C:\\Users\\sasha\\Box\\third_year\\venus_earlyhab\\volccalc\\VolcGases-master\\results\\repaired_all\\old\\tstart5_gwater500_H2O0_CO2300_f_ext10_f_volc9.csv",
    #                   usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),header=None)
    # f1 = pd.read_csv("C:\\Users\\sasha\\Box\\third_year\\venus_earlyhab\\volccalc\\VolcGases-master\\results\\no_melting\\old\\NM_tstart5_gwater100_H2O0_CO2300_f_ext10_f_volc9.csv",
    #                     usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),header=None)
    if melt == 0:
        f1 = pd.read_csv("C:\\Users\\sasha\\Documents\\results\\nomelt_gwat100_t5\\NMtstart5_gwater100_H2O0_CO2500_f_ext10_f_volc9_FMQ-1.csv",
                        usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27),header=None)
    elif melt == 1:
        f1 = pd.read_csv("C:\\Users\\sasha\\Documents\\results\\melt_wd_gwat500_5\\WDtstart5_gwater500_H2O0_CO2500_f_ext10_f_volc9_FMQ-1.csv",
                        usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27),header=None)
        
        
    
        
    results = f1.to_numpy()
        
    i = results[:,0]
    i = i.astype(int)
    i1 = i[-1]
    
    final_tstep = results[:,0]
    
    tstart = results[0,1]
    # t_arr = np.unique(tstart)
    
    t = results[:,2]
    
    f_ext = results[0,3]
    # fext_arr = np.unique(f_ext)
    
    f_volc = results[0,4]
    # fvolc_arr = np.unique(f_volc)
    
    z_erupt = results[:,6]
    
    gwater = results[0,8]
    # gwat_arr = np.unique(gwater)
    
    mCO2tot = results[0,9]
    # CO2_arr = np.unique(mCO2tot)
    CO2_plotarr = np.array([300,500,1000,2000])
    
    mH2Otot = results[0,10]
    if mH2Otot>0.6e-2:
        mH2Otot = np.round(mH2Otot,3)
    
    FMQ = results[0,11]
    
    
    
    
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
    z_erupt = ((moles_CO2*Rmm_CO2)-((1-f_volc)*(90e5/g)))/(1000*g*rho*mCO2tot)
    moles_CH4 = results[:,27]
    # success_Ar = results[:,28]
    
    
    capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
    
    if melt == 1:
        gwater_string = str(int(gwater))
        h2o_string = str(int(1000*mH2Otot))
        f_ext_string = str(int(10*f_ext))
        f_volc_string = str(int(10*f_volc))
        CO2_string = str(int(1e6*mCO2tot))
        t_string = str(int(10*tstart))
        FMQ_string = str(int(FMQ))
        
        title = "C:\\Users\\sasha\\Box\\Clean Code\\results\\reduced_melt\\withdiss_t" + t_string + "_gwater" + gwater_string + "_f_volc" + f_volc_string + ".csv"
        
        melt_res = pd.read_csv(title)
    
        diss_dat = melt_res.to_numpy()
        
        time_d = diss_dat[:,4]
        Pres_d = diss_dat[:,5]
        moles_atm_d = diss_dat[:,6]
        moles_CO2_d = diss_dat[:,7]
        moles_H2O_d = diss_dat[:,8]
        moles_O2_d = diss_dat[:,9]
        nonthermal_d = diss_dat[:,10]
        magma_d = diss_dat[:,11]
        Rmm_atm_d = diss_dat[:,12]
        
        blanks = np.zeros(np.size(time_d))
        
        Pres = np.append(Pres_d,Pres)
        moles_atm = np.append(moles_atm_d,moles_atm)
        moles_CO2 = np.append(moles_CO2_d,moles_CO2)
        moles_CO = np.append(blanks,moles_CO)
        moles_H2 = np.append(blanks,moles_H2)
        moles_H2O = np.append(moles_H2O_d,moles_H2O)
        moles_O2 = np.append(moles_O2_d,moles_O2)
        Rmm_atm = np.append(Rmm_atm_d,Rmm_atm)
        T_surf = np.append(blanks,T_surf)
        ash = np.append(blanks,ash)
        lava = np.append(blanks,lava)
        magma = np.append(magma_d,magma)
        nonthermal = np.append(nonthermal_d,nonthermal)
        redox = np.append(blanks,redox)
        z_erupt = np.append(blanks,z_erupt)
        moles_CH4 = np.append(blanks,moles_CH4)
    # success_Ar = results[:,28]
    
    
        capN = capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)

    ax1[0].axhline(y=(100.e-6 * moles_atm[i1] * Rmm_H2O),color ='midnightblue', linestyle='--',linewidth=15)
    ax1[1].axhline(y=(50.e-6 * moles_atm[i1] * Rmm_O2), color ='lightseagreen', linestyle='--',linewidth=15)
    
    plot_colmass, = ax1[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=15,color='black',linestyle=':')
    plot_colCO2, = ax1[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=15,color ='red')
    plot_colH2O, = ax1[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=15,color ='midnightblue')
    # plot_colO2, = ax1[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=15,color ='lightseagreen')
    
    plot_colmass, = ax1[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=15,color='black',linestyle=':')
    plot_colCO2, = ax1[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=15,color ='red')
    # plot_colH2O, = ax1[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=15,color ='midnightblue')
    plot_colO2, = ax1[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=15,color ='teal')
    
    
    
    ax1[0].axhline(y=(100.e-6 * moles_atm[i1] * Rmm_H2O),color ='midnightblue', linestyle='--',linewidth=15)
    ax1[1].axhline(y=(50.e-6 * moles_atm[i1] * Rmm_O2), color ='teal', linestyle='--',linewidth=15)
    
    ax1[0].text(1.4,(100.e-6 * moles_atm[i1] * Rmm_H2O)+15,'Modern Venus H$\mathregular{_{2}}$O', fontsize=90,weight='bold',color='midnightblue')
    ax1[1].text(1.3,(50.e-6 * moles_atm[i1] * Rmm_O2)+15,'Modern Venus O$\mathregular{_{2}}$', fontsize=90,weight='bold',color='teal')
    
    ax1[0].text(4.1,1e6,'i. H$\mathregular{_{2}}$O', fontsize=120,weight='bold')
    ax1[1].text(4.1,1e6,'ii. O$\mathregular{_{2}}$', fontsize=120,weight='bold')
    # plt.axhline(y=1, color ='lightgrey', linestyle='--',linewidth=2)
    
    # plt.xlim(2,2.1)
    # plt.legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=40)
    ax1[0].legend(handles=[plot_colmass,plot_colCO2,plot_colH2O], fontsize=100,loc='center right')
    ax1[1].legend(handles=[plot_colmass,plot_colCO2,plot_colO2], fontsize=100,loc='center right')
    # ax1[1].legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=100)
    if melt == 0:
        
        ax1[0].set_title('$\mathbf{a. No \; runaway \; greenhouse \; surface \; melting}$' '\n' + str(float( '%g' % (gwater))) + ' m GEL, '+ str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))) + ', $\mathregular{\Delta}$FMQ = ' + str(float('%.1g' % (FMQ))), fontsize=110) 
    # ax1[0].set_title('$\mathbf{b. With\;runaway\;greenhouse\;surface\;melting}$' '\n' + str(float( '%g' % (gwater))) + ' m GEL, '+ str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))), fontsize=110) 
    elif melt == 1 :
        ax1[0].set_title('$\mathbf{b. Runaway \; greenhouse \; surface \; melting}$' '\n' + str(float( '%g' % (gwater))) + ' m GEL, '+ str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))) + ', $\mathregular{\Delta}$FMQ = ' + str(float('%.1g' % (FMQ))), fontsize=110) 
    
    
    ax1[0].annotate('', xy=(3.5,10),  xycoords='data',
            xytext=(3.5, (100.e-6 * moles_atm[i1] * Rmm_H2O)), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[0].annotate('', xy=(2.5,10),  xycoords='data',
            xytext=(2.5, (100.e-6 * moles_atm[i1] * Rmm_H2O)), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[0].annotate('', xy=(1.5,10),  xycoords='data',
            xytext=(1.5, (100.e-6 * moles_atm[i1] * Rmm_H2O)), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[0].annotate('', xy=(0.5,10),  xycoords='data',
            xytext=(0.5, (100.e-6 * moles_atm[i1] * Rmm_H2O)), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    
    ax1[1].annotate('', xy=(3.5,10),  xycoords='data',
            xytext=(3.5, (50.e-6 * moles_atm[i1] * Rmm_O2)), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[1].annotate('', xy=(2.5,10),  xycoords='data',
            xytext=(2.5, (50.e-6 * moles_atm[i1] * Rmm_O2)), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[1].annotate('', xy=(1.5,10),  xycoords='data',
            xytext=(1.5, (50.e-6 * moles_atm[i1] * Rmm_O2)), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[1].annotate('', xy=(0.5,10),  xycoords='data',
            xytext=(0.5, (50.e-6 * moles_atm[i1] * Rmm_O2)), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    
#%%

if plot_number == 2:
    
    file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\nomelt_gwat300_t5\\NMtstart*_gwater300_H2O0_CO2500_f_ext10_f_volc9_FMQ-4.csv")
    # file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\melt_wd_gwat500_5\\WDtstart*_gwater500_H2O0_CO2500_f_ext10_f_volc9_FMQ-4.csv")
    
    # print(file_list)
    #%%
    fig2, ax2 = plt.subplots(figsize=(50,40))
    
    ax2.set_ylim((1.e1,3.e6))

    
    ax2.set_yscale("log")
    ax2.invert_xaxis()

    ax2.set_ylabel('Column Mass (kg)', fontsize=100,weight='bold')
    ax2.set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    ax2.tick_params(axis='both', which='major', labelsize=100)

    
    for x in file_list:
        counter = counter + 1
        # print('counter',counter)
        
        # if counter == 3:
        #     break
        # if i2 == 3:
        #     break
        
        f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27),header=None)
        
        results = f1.to_numpy()
        
        i = results[:,0]
        i = i.astype(int)
        i1 = i[-1]
        
        final_tstep = results[:,0]
        
        tstart = results[0,1]
        # t_arr = np.unique(tstart)
        
        t = results[:,2]
        
        f_ext = results[0,3]
        # fext_arr = np.unique(f_ext)
        
        f_volc = results[0,4]
        # fvolc_arr = np.unique(f_volc)
        
        z_erupt = results[:,6]
        
        gwater = results[0,8]
        gwat_arr = np.unique(gwater)
        
        mCO2tot = results[0,9]
        CO2_arr = np.unique(mCO2tot)
        CO2_plotarr = np.array([300,500,1000,2000])
        
        mH2Otot = results[0,10]
        if mH2Otot > 0.6e-2:
            mH2Otot = np.round(results[0,10],3)
        
        FMQ = results[0,11]
        
        
        
        
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
        z_erupt = ((moles_CO2*Rmm_CO2)-((1-f_volc)*(90e5/g)))/(1000*g*rho*mCO2tot)
        moles_CH4 = results[:,27]
        # success_Ar = results[:,28]
        
        
        capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
        
        
        z_erupt = (((moles_CO2 - np.append(0,moles_CO2[0:-1]))*Rmm_CO2)/(g*rho*mCO2tot))/((t - np.append(0,t[0:-1]))/secingyr)
        CO2_mass_er = np.cumsum(z_erupt)*(t[1]-t[0])*rho*mCO2tot*g/secingyr
        z_erupt[T_surf>0]==0
        
        if melt == 1:
            gwater_string = str(int(gwater))
            h2o_string = str(int(1000*mH2Otot))
            f_ext_string = str(int(10*f_ext))
            f_volc_string = str(int(10*f_volc))
            CO2_string = str(int(1e6*mCO2tot))
            t_string = str(int(10*tstart))
            FMQ_string = str(int(FMQ))
            
            title = "C:\\Users\\sasha\\Box\\Clean Code\\results\\reduced_melt\\withdiss_t" + t_string + "_gwater" + gwater_string + "_f_volc" + f_volc_string + ".csv"
            
            melt_res = pd.read_csv(title)
        
            diss_dat = melt_res.to_numpy()
            
            time_d = diss_dat[:,4]
            Pres_d = diss_dat[:,5]
            moles_atm_d = diss_dat[:,6]
            moles_CO2_d = diss_dat[:,7]
            moles_H2O_d = diss_dat[:,8]
            moles_O2_d = diss_dat[:,9]
            nonthermal_d = diss_dat[:,10]
            magma_d = diss_dat[:,11]
            Rmm_atm_d = diss_dat[:,12]
            
            blanks = np.zeros(np.size(time_d))
            
            Pres = np.append(Pres_d,Pres)
            moles_atm = np.append(moles_atm_d,moles_atm)
            moles_CO2 = np.append(moles_CO2_d,moles_CO2)
            moles_CO = np.append(blanks,moles_CO)
            moles_H2 = np.append(blanks,moles_H2)
            moles_H2O = np.append(moles_H2O_d,moles_H2O)
            moles_O2 = np.append(moles_O2_d,moles_O2)
            Rmm_atm = np.append(Rmm_atm_d,Rmm_atm)
            T_surf = np.append(blanks,T_surf)
            ash = np.append(blanks,ash)
            lava = np.append(blanks,lava)
            magma = np.append(magma_d,magma)
            nonthermal = np.append(nonthermal_d,nonthermal)
            redox = np.append(blanks,redox)
            z_erupt = np.append(blanks,z_erupt)
            moles_CH4 = np.append(blanks,moles_CH4)
        # success_Ar = results[:,28]
        
        
            capN = capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)

        # if all(t_arr != tstart):
        #     continue
        # if all(CO2_arr != mCO2tot):
        #     continue
        # if all(H2O_arr != mH2Otot):
        #     continue
        # if all(gwat_arr != gwater):
        #     continue
        # if all(fext_arr != f_ext):
        #     continue
        # if all(fvolc_arr != f_volc):
        #     continue 
        
        ax2.axvline(x=(4.5-tstart), color ='grey', linestyle=':',linewidth=15)
        
        if tstart == 0.5:
            # plot_colmass, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=15,color='black')
            plot_colCO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=15,color ='red')
            # plot_colCO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),CO2_mass_er[0:i1],label='CO\N{SUBSCRIPT TWO}',linewidth=15,color ='red',linestyle='-')
            plot_colH2O, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=15,color ='midnightblue')
            plot_colO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=15,color ='lightseagreen')
            plot_colCO, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO[0:i1]*Rmm_CO,label='CO',linewidth=15,color ='darkorange')
            plot_colCH4, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CH4[0:i1]*Rmm_CH4,label='CH\N{SUBSCRIPT FOUR}',linewidth=15,color ='hotpink')
            plot_colH2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2[0:i1]*Rmm_H2,label='H\N{SUBSCRIPT TWO}',linewidth=15,color ='mediumorchid')
            
            
        if tstart == 1.5:
            # plot_colmass, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=15,color='black',linestyle='--')
            plot_colCO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=15,color ='red',linestyle='--')
            # plot_colCO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),CO2_mass_er[0:i1],label='CO\N{SUBSCRIPT TWO}',linewidth=15,color ='red',linestyle='--')
            plot_colH2O, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=15,color ='midnightblue',linestyle='--')
            plot_colO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=15,color ='lightseagreen',linestyle='--')
            plot_colCO, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO[0:i1]*Rmm_CO,label='CO',linewidth=15,color ='darkorange',linestyle='--')
            plot_colCH4, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CH4[0:i1]*Rmm_CH4,label='CH\N{SUBSCRIPT FOUR}',linewidth=15,color ='hotpink',linestyle='--')
            plot_colH2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2[0:i1]*Rmm_H2,label='H\N{SUBSCRIPT TWO}',linewidth=15,color ='mediumorchid',linestyle='--')
            
        ax2.annotate('', xy=(3.5,10),  xycoords='data',
            xytext=(3.5, (50.e-6 * moles_atm[i1] * Rmm_O2)), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
        ax2.annotate('', xy=(2.5,10),  xycoords='data',
            xytext=(2.5, (100.e-6 * moles_atm[i1] * Rmm_H2O)), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
        ax2.annotate('', xy=(1.5,10),  xycoords='data',
            xytext=(1.5, (50.e-6 * moles_atm[i1] * Rmm_O2)), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
        ax2.annotate('', xy=(-0.15,10),  xycoords='data',
            xytext=(-0.15, (100.e-6 * moles_atm[i1] * Rmm_H2O)), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
        
    #%%
    
    ax2.axhline(y=(100.e-6 * moles_atm[i1] * Rmm_H2O),color ='midnightblue', linestyle='--',linewidth=15)
    ax2.axhline(y=(50.e-6 * moles_atm[i1] * Rmm_O2), color ='teal', linestyle='--',linewidth=15)
    
    ax2.text(1.4,(100.e-6 * moles_atm[i1] * Rmm_H2O)+15,'Modern Venus H$\mathregular{_{2}}$O', fontsize=90,weight='bold',color='midnightblue')
    ax2.text(1.3,(50.e-6 * moles_atm[i1] * Rmm_O2)-15,'Modern Venus O$\mathregular{_{2}}$', fontsize=90,weight='bold',color='teal')
    
    ax2.text(3.95,1e6,'$\mathregular{t_{0}}=4.0$', fontsize=90,color='grey')
    ax2.text(2.95,1e6,'$\mathregular{t_{0}}=3.0$', fontsize=90,color='grey')
    # plt.axhline(y=1, color ='lightgrey', linestyle='--',linewidth=2)
    
    # plt.xlim(2,2.1)
    # plt.legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=40)
    ax2.legend(handles=[plot_colO2,plot_colCO2,plot_colCO,plot_colCH4,plot_colH2O,plot_colH2], fontsize=100,loc='best',bbox_to_anchor=(1.22, 0.8))
    ax2.set_title('$\mathbf{No \; runaway \; greenhouse \; surface \; melting}$' '\n' + str(float( '%g' % (gwater))) + ' m GEL, '+ str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))) + ', $\mathregular{\Delta}$FMQ = ' + str(float('%.1g' % (FMQ))), fontsize=100) 

    
if plot_number == 3:
    file_list = glob.glob("C:\\Users\\sasha\\Box\\third_year\\venus_earlyhab\\volccalc\\VolcGases-master\\results\\no_melting\\old\\NM_tstart*_gwater100_H2O0_CO2*f_ext10_f_volc10.csv")
    counter = 0
    
    
    fig3, ax3 = plt.subplots(figsize=(50,40))
    
    # ax3.set_ylim((1.e1,3.e6))

    
    # ax3.set_yscale("log")
    ax3.invert_xaxis()

    ax3.set_ylabel('Eruption rate (km Gyr$\mathregular{^{-1}}$)', fontsize=100,weight='bold')
    ax3.set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    ax3.tick_params(axis='both', which='major', labelsize=100)

    for x in file_list:
        counter = counter + 1
        # print('counter',counter)
        
        if counter == 100:
            break
        # if i2 == 3:
        #     break
        
        f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27),header=None)
        
        results = f1.to_numpy()
        
        i = results[:,0]
        i = i.astype(int)
        i1 = i[-1]
        
        final_tstep = results[:,0]
        
        tstart = results[0,1]
        # t_arr = np.unique(tstart)
        
        t = results[:,2]
        
        f_ext = results[0,3]
        # fext_arr = np.unique(f_ext)
        
        f_volc = results[0,4]
        # fvolc_arr = np.unique(f_volc)
        
        z_erupt = results[:,6]
        
        gwater = results[0,8]
        gwat_arr = np.unique(gwater)
        
        mCO2tot = results[0,9]
        CO2_arr = np.unique(mCO2tot)
        CO2_plotarr = np.array([300,500,1000,2000])
        
        mH2Otot = results[0,10]
        if mH2Otot > 0.6e-2:
            mH2Otot = np.round(results[0,10],3)
        
        FMQ = results[0,11]
        
        
        
        
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
        z_erupt = ((moles_CO2*Rmm_CO2)-((1-f_volc)*(90e5/g)))/(1000*g*rho*mCO2tot)
        moles_CH4 = results[:,27]
        # success_Ar = results[:,28]
        
        
        capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
        
        
        
        if all(t_arr != tstart):
            continue
        if all(CO2_arr != mCO2tot):
            continue
        if all(H2O_arr != mH2Otot):
            continue
        if all(gwat_arr != gwater):
            continue
        if all(fext_arr != f_ext):
            continue
        if all(fvolc_arr != f_volc):
            continue 
        
        if mCO2tot == 300e-6:
            lst = str(':')
            
        if mCO2tot == 500e-6:
            lst = str('--')
            
        if mCO2tot == 1000e-6:
            lst = str('-.')
            
        if mCO2tot == 2000e-6:
            lst = str('-')
            
            
        if tstart == 3:
            clr = str('limegreen')
            
        if tstart == 1.5:
            clr = str('lightseagreen')
        if tstart == 0.5:
            clr = str('midnightblue')
        
        ax3.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),z_erupt[0:i1]/1000,linewidth=15,color=clr,linestyle=lst)
        
        
        
        print(moles_CO2[-1])
        
        
if plot_number == 4:
    if melt == 0:
        file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\nomelt_gwat100_t5\\NMtstart5_gwater100_H2O6_CO21000_f_ext3_f_volc9_FMQ*.csv")
    elif melt == 1:
        # file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\melt_wd_gwat500_5\\WDtstart5_gwater500_H2O2_CO2500_f_ext10_f_volc9_FMQ*.csv")
        file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\melt_wd_gwat1000\\WDtstart5_gwater1000_H2O6_CO21000_f_ext3_f_volc9_FMQ*.csv")
    
    
    fig4, ax4 = plt.subplots(2,figsize=(50,70))
    
    ax4[0].set_ylim((1.e1,3.e6))
    ax4[0].set_yscale("log")
    ax4[0].invert_xaxis()

    ax4[0].set_ylabel('Column Mass O$\mathregular{_{2}}$ (kg)', fontsize=100,weight='bold')
    ax4[0].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    ax4[0].tick_params(axis='both', which='major', labelsize=100)
    
    ax4[1].set_ylim((1.e1,3.e6))

    
    ax4[1].set_yscale("log")
    ax4[1].invert_xaxis()

    ax4[1].set_ylabel('Column Mass C species (kg)', fontsize=100,weight='bold')
    ax4[1].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    ax4[1].tick_params(axis='both', which='major', labelsize=100)
    
    # ax4[2].set_ylim((1.e1,3.e6))

    # ax4[2].set_yscale("log")
    # ax4[2].invert_xaxis()

    # ax4[2].set_ylabel('Column Mass H species (kg)', fontsize=100,weight='bold')
    # ax4[2].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    # ax4[2].tick_params(axis='both', which='major', labelsize=100)
    
    
    
    for x in file_list:
        counter = counter + 1
        # print('counter',counter)
        
        # if counter == 3:
        #     break
        # if i2 == 3:
        #     break
        
        f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27),header=None)
        
    
        results = f1.to_numpy()
            
        i = results[:,0]
        i = i.astype(int)
        i1 = i[-1]
        
        final_tstep = results[:,0]
        
        tstart = results[0,1]
        # t_arr = np.unique(tstart)
        
        t = results[:,2]
        
        f_ext = results[0,3]
        # fext_arr = np.unique(f_ext)
        
        f_volc = results[0,4]
        # fvolc_arr = np.unique(f_volc)
        
        z_erupt = results[:,6]
        
        gwater = results[0,8]
        gwat_arr = np.unique(gwater)
        
        mCO2tot = results[0,9]
        CO2_arr = np.unique(mCO2tot)
        CO2_plotarr = np.array([300,500,1000,2000])
        
        mH2Otot = results[0,10]
        if mH2Otot > 0.6e-2:
            mH2Otot = np.round(results[0,10],3)
        
        FMQ = results[0,11]
        
        
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
        z_erupt = ((moles_CO2*Rmm_CO2)-((1-f_volc)*(90e5/g)))/(1000*g*rho*mCO2tot)
        moles_CH4 = results[:,27]
        # success_Ar = results[:,28]
        
        
        capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
        
        
        z_erupt = (((moles_CO2 - np.append(0,moles_CO2[0:-1]))*Rmm_CO2)/(g*rho*mCO2tot))/((t - np.append(0,t[0:-1]))/secingyr)
        CO2_mass_er = np.cumsum(z_erupt)*(t[1]-t[0])*rho*mCO2tot*g/secingyr
        z_erupt[T_surf>0]==0
        
        if melt == 1:
            gwater_string = str(int(gwater))
            h2o_string = str(int(1000*mH2Otot))
            f_ext_string = str(int(10*f_ext))
            f_volc_string = str(int(10*f_volc))
            CO2_string = str(int(1e6*mCO2tot))
            t_string = str(int(10*tstart))
            FMQ_string = str(int(FMQ))
            
            title = "C:\\Users\\sasha\\Box\\Clean Code\\results\\reduced_melt\\withdiss_t" + t_string + "_gwater" + gwater_string + "_f_volc" + f_volc_string + ".csv"
            
            melt_res = pd.read_csv(title)
        
            diss_dat = melt_res.to_numpy()
            
            time_d = diss_dat[:,4]
            Pres_d = diss_dat[:,5]
            moles_atm_d = diss_dat[:,6]
            moles_CO2_d = diss_dat[:,7]
            moles_H2O_d = diss_dat[:,8]
            moles_O2_d = diss_dat[:,9]
            nonthermal_d = diss_dat[:,10]
            magma_d = diss_dat[:,11]
            Rmm_atm_d = diss_dat[:,12]
            
            blanks = np.zeros(np.size(time_d))
            
            Pres = np.append(Pres_d,Pres)
            moles_atm = np.append(moles_atm_d,moles_atm)
            moles_CO2 = np.append(moles_CO2_d,moles_CO2)
            moles_CO = np.append(blanks,moles_CO)
            moles_H2 = np.append(blanks,moles_H2)
            moles_H2O = np.append(moles_H2O_d,moles_H2O)
            moles_O2 = np.append(moles_O2_d,moles_O2)
            Rmm_atm = np.append(Rmm_atm_d,Rmm_atm)
            T_surf = np.append(blanks,T_surf)
            ash = np.append(blanks,ash)
            lava = np.append(blanks,lava)
            magma = np.append(magma_d,magma)
            nonthermal = np.append(nonthermal_d,nonthermal)
            redox = np.append(blanks,redox)
            z_erupt = np.append(blanks,z_erupt)
            moles_CH4 = np.append(blanks,moles_CH4)
        # success_Ar = results[:,28]
        
        
            capN = capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
        
        ax4[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_O2[0:i1]*Rmm_O2,label=('$\mathregular{\Delta}$FMQ = ' + str(FMQ)),linewidth=15,color = cm.gnuplot((np.abs(FMQ))/(np.max(1.2*np.abs(FMQ_arr)))))
        
        ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO2[0:i1]*Rmm_CO2,linewidth=15,color = cm.gnuplot(np.abs(FMQ)/np.max(1.2*np.abs(FMQ_arr))))
        ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO[0:i1]*Rmm_CO,linewidth=15,color = cm.gnuplot(np.abs(FMQ)/np.max(1.2*np.abs(FMQ_arr))),linestyle='--')
        ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CH4[0:i1]*Rmm_CH4,linewidth=15,color = cm.gnuplot(np.abs(FMQ)/np.max(1.2*np.abs(FMQ_arr))),linestyle=':')
        
        if counter == 1:
            ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),1e-10*moles_CO2[0:i1]*Rmm_CO2,label = 'CO$\mathregular{_{2}}$', linewidth=15,color = 'black')
            ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),1e-10*moles_CO[0:i1]*Rmm_CO,label = 'CO',linewidth=15,color ='black',linestyle='--')
            ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),1e-10*moles_CH4[0:i1]*Rmm_CH4,label = 'CH$\mathregular{_{4}}$',linewidth=15,color ='black',linestyle=':')
        
        # ax4[2].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=15,color = cm.gnuplot(np.abs(FMQ)/np.max(1.2*np.abs(FMQ_arr))))
        # ax4[2].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2[0:i1]*Rmm_H2,label='H\N{SUBSCRIPT TWO}',linewidth=15,color = cm.gnuplot(np.abs(FMQ)/np.max(1.2*np.abs(FMQ_arr))),linestyle='--')

    
    ax4[1].axhline(y=(50.e-6 * (93e5/(g*Rmm_CO2)) * Rmm_CO),color ='black', linestyle='--',linewidth=15)
    ax4[0].axhline(y=(50.e-6 * (93e5/(g*Rmm_CO2)) * Rmm_O2), color ='black', linestyle='--',linewidth=15)
    
    if melt == 0:
        
        ax4[0].text(4.1,1e6,'i.', fontsize=120,weight='bold')
        ax4[1].text(4.1,1e6,'ii.', fontsize=120,weight='bold')
    
    if melt == 1:
        ax4[0].text(3.9,1e6,'i.', fontsize=120,weight='bold')
        ax4[1].text(3.9,1e6,'ii.', fontsize=120,weight='bold')
        ax4[0].text(3.9,10**5.5,'Surface molten', fontsize=90,weight='bold',color='grey')
        ax4[1].text(3.9,10**5.5,'Surface molten', fontsize=90,weight='bold',color='grey')
        ax4[0].legend(fontsize=100,loc='best',bbox_to_anchor=(0.9, 0.8))
        ax4[1].legend(fontsize=100,loc='best',bbox_to_anchor=(1, 0.8))
    
    ax4[1].text(1.4,(50.e-6 * (93e5/(g*Rmm_CO2)) * Rmm_CO)+15,'Modern Venus CO', fontsize=90,weight='bold',color='black')
    ax4[0].text(1.3,(50.e-6 * (93e5/(g*Rmm_CO2)) * Rmm_O2)+15,'Modern Venus O$\mathregular{_{2}}$', fontsize=90,weight='bold',color='black')
    
    
    
    
    ax4[0].annotate('', xy=(3.5,10),  xycoords='data',
        xytext=(3.5, (50.e-6 * (93e5/(g*Rmm_CO2)) * Rmm_O2)), textcoords='data',
        arrowprops=dict(facecolor='black', shrink=0,width=20,headwidth=50,headlength=60),
        horizontalalignment='left', verticalalignment='bottom',
        )
    ax4[0].annotate('', xy=(1.5,10),  xycoords='data',
        xytext=(1.5, (50.e-6 *(93e5/(g*Rmm_CO2)) * Rmm_O2)), textcoords='data',
        arrowprops=dict(facecolor='black', shrink=0,width=20,headwidth=50,headlength=60),
        horizontalalignment='left', verticalalignment='bottom',
        )
    
    ax4[1].annotate('', xy=(3.5,10),  xycoords='data',
        xytext=(3.5, (50.e-6 * (93e5/(g*Rmm_CO2)) * Rmm_CO)), textcoords='data',
        arrowprops=dict(facecolor='black', shrink=0,width=20,headwidth=50,headlength=60),
        horizontalalignment='left', verticalalignment='bottom',
        )
    ax4[1].annotate('', xy=(1.5,10),  xycoords='data',
        xytext=(1.5, (50.e-6 * (93e5/(g*Rmm_CO2)) * Rmm_CO)), textcoords='data',
        arrowprops=dict(facecolor='black', shrink=0,width=20,headwidth=50,headlength=60),
        horizontalalignment='left', verticalalignment='bottom',
        )
    
    if melt == 0:
        ax4[0].set_title('$\mathbf{a. No \; runaway \; greenhouse \; surface \; melting - }$' + r'$\mathbf{' + str(float( '%g' % (gwater))) + '}$' + ' $\mathbf{m \;GEL,}$ ' + '\n'  + str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))), fontsize=100) 
    elif melt == 1:
        ax4[0].set_title('$\mathbf{b. Runaway \; greenhouse \; surface \; melting - }$' + r'$\mathbf{' + str(float( '%g' % (gwater))) + '}$' + ' $\mathbf{m \;GEL,}$ ' + '\n'  + str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))), fontsize=100) 
    
if plot_number == 5:
    
    option = 1
    
    max_O2_modern = 50.e-6*93.e5/(g*Rmm_CO2) 
    max_H2O_modern = 100.e-6*93.e5/(g*Rmm_CO2) 
    max_CO_modern = 50.e-6*93.e5/(g*Rmm_CO2) 
    
    f_volc_var = np.array([1,5,9,10])
    
    if melt == 0:
        file_list = glob.glob("C:\\Users\\sasha\\Box\\Clean Code\\results\\sensitivity\\varyb*NMtstart15_gwater500_H2O5_CO2300_f_ext3_f_volc10_FMQ-2.csv")
    elif melt == 1:
        # file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\melt_wd_gwat500_5\\WDtstart5_gwater500_H2O2_CO2500_f_ext10_f_volc9_FMQ*.csv")
        file_list = glob.glob("C:\\Users\\sasha\\Documents\\results\\melt_wd_gwat1000\\WDtstart5_gwater1000_H2O6_CO21000_f_ext3_f_volc9_FMQ*.csv")
    
    
    fig4, ax4 = plt.subplots(3,figsize=(50,105))
    
    
    ax4[2].set_xlim((-0.1,3.1))
    ax4[0].invert_xaxis()

    
    
    # ax4[1].set_ylim((0,1e4))
    # ax4[1].set_ylim((1.e1,3.e6))

    
    # ax4[1].set_yscale("log")
    ax4[1].invert_xaxis()
    ax4[2].invert_xaxis()

    if option == 1:
        
        ax4[0].set_ylabel('% Difference in  [O$\mathregular{_{2}}$]$\mathregular{_{atm}}$', fontsize=100,weight='bold')
        ax4[0].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[0].tick_params(axis='both', which='major', labelsize=100) 
        
        ax4[1].set_ylabel('% Difference in  [H$\mathregular{_{2}}$O]$\mathregular{_{atm}}$', fontsize=100,weight='bold')
        ax4[1].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[1].tick_params(axis='both', which='major', labelsize=100)
        
        ax4[2].set_ylabel('% Difference in [CO]$\mathregular{_{atm}}$', fontsize=100,weight='bold')
        ax4[2].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[2].tick_params(axis='both', which='major', labelsize=100)
    
    if option == 2:
        
        ax4[0].set_ylabel('[O$\mathregular{_{2}}$]$\mathregular{_{model}}$ / [O$\mathregular{_{2}}$]$\mathregular{_{modern}}$', fontsize=100,weight='bold')
        ax4[0].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[0].tick_params(axis='both', which='major', labelsize=100)
        
        ax4[1].set_ylabel('[H$\mathregular{_{2}}$O]$\mathregular{_{model}}$ / [H$\mathregular{_{2}}$O]$\mathregular{_{modern}}$ ', fontsize=100,weight='bold')
        ax4[1].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[1].tick_params(axis='both', which='major', labelsize=100)
        
        ax4[2].set_ylabel('[CO]$\mathregular{_{model}}$ / [CO]$\mathregular{_{modern}}$', fontsize=100,weight='bold')
        ax4[2].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[2].tick_params(axis='both', which='major', labelsize=100)
        
        ax4[0].axhline(y=1,color ='black', linestyle='--',linewidth=15)
        ax4[1].axhline(y=1,color ='black', linestyle='--',linewidth=15)
        ax4[2].axhline(y=1,color ='black', linestyle='--',linewidth=15)
        
        ax4[0].set_yscale("log")
        ax4[1].set_yscale("log")
        ax4[2].set_yscale("log")
        ax4[2].set_ylim((1e-5,10))
    
        
        
    
    # if melt == 0:
    if option == 1:
        
        ax4[0].text(3.1,0.03,'a.', fontsize=120,weight='bold')
        ax4[1].text(3.1,45,'b.', fontsize=120,weight='bold')
        ax4[2].text(3.1,0.045,'c.', fontsize=120,weight='bold')
       
    if option == 2:
        
        ax4[0].text(3.1,10**4,'d.', fontsize=120,weight='bold')
        ax4[1].text(3.1,10**4.1,'e.', fontsize=120,weight='bold')
        ax4[2].text(3.1,5,'f.', fontsize=120,weight='bold')
    
    
    # ax4[2].invert_xaxis()

    # ax4[2].set_ylabel('Column Mass H species (kg)', fontsize=100,weight='bold')
    # ax4[2].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
    # ax4[2].tick_params(axis='both', which='major', labelsize=100)
    
    
    for i_n in range(0,np.size(f_volc_var)):
        
        f_volc_pick = f_volc_var[i_n]
        
        file_list = glob.glob("C:\\Users\\sasha\\Box\\Clean Code\\results\\sensitivity\\varyb*NMtstart15_gwater500_H2O5_CO2500_f_ext3*f_volc" + str(f_volc_pick) + "_FMQ-3.csv")
        
        if f_volc_pick == 1:
            ln_str = '-'
        if f_volc_pick == 5:
            ln_str = '--'
        if f_volc_pick == 9:
            ln_str = ':'
        if f_volc_pick == 10:
            ln_str = '-.'
        
        counter = 0
        
        if np.size(file_list) == 3:
            for x in file_list:
                counter = counter + 1
                print('counter',counter)
                
                # if counter == 3:
                #     break
                # if i2 == 3:
                #     break
                
                f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28),header=None)
                
            
                results = f1.to_numpy()
                    
                i = results[:,0]
                i = i.astype(int)
                i1 = i[-1]
                
                final_tstep = results[:,0]
                
                tstart = results[0,1]
                # t_arr = np.unique(tstart)
                
                t = results[:,2]
                
                f_ext = results[0,3]
                # fext_arr = np.unique(f_ext)
                
                f_volc = results[0,4]
                # fvolc_arr = np.unique(f_volc)
                
                z_erupt = results[:,6]
                
                gwater = results[0,8]
                gwat_arr = np.unique(gwater)
                
                mCO2tot = results[0,9]
                CO2_arr = np.unique(mCO2tot)
                CO2_plotarr = np.array([300,500,1000,2000])
                
                mH2Otot = results[0,10]
                if mH2Otot > 0.6e-2:
                    mH2Otot = np.round(results[0,10],3)
                
                FMQ = results[0,11]
                
                
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
                z_erupt = ((moles_CO2*Rmm_CO2)-((1-f_volc)*(90e5/g)))/(1000*g*rho*mCO2tot)
                moles_CH4 = results[:,27]
                bmulti = results[0,28]
                # success_Ar = results[:,28]
                
                print('final O',moles_O2[-1])
                
                if counter == 1:
                    if option == 2:
                        
                        O2_1 = moles_O2/max_O2_modern
                        H2O_1 = moles_H2O/max_H2O_modern
                        CO_1 = moles_CO/max_CO_modern
                        
                    if option == 1:
                        O2_1 = moles_O2/moles_atm
                        H2O_1 = moles_H2O/moles_atm
                        CO_1 = np.nan_to_num(moles_CO/moles_atm)
                elif counter == 2:
                    if option == 2:
                        O2_2 = moles_O2/max_O2_modern
                        H2O_2 = moles_H2O/max_H2O_modern
                        CO_2 = moles_CO/max_CO_modern
                    if option == 1:
                        O2_2 = moles_O2/moles_atm
                        H2O_2 = moles_H2O/moles_atm
                        CO_2 = np.nan_to_num(moles_CO/moles_atm)
                elif counter == 3:
                    if option == 2:
                        O2_3 = moles_O2/max_O2_modern
                        H2O_3 = moles_H2O/max_H2O_modern
                        CO_3 = moles_CO/max_CO_modern
                    if option == 1:
                        O2_3 = moles_O2/moles_atm
                        H2O_3 = moles_H2O/moles_atm
                        CO_3 = np.nan_to_num(moles_CO/moles_atm)
                
                capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
                
                
                z_erupt = (((moles_CO2 - np.append(0,moles_CO2[0:-1]))*Rmm_CO2)/(g*rho*mCO2tot))/((t - np.append(0,t[0:-1]))/secingyr)
                CO2_mass_er = np.cumsum(z_erupt)*(t[1]-t[0])*rho*mCO2tot*g/secingyr
                z_erupt[T_surf>0]==0
                
                if option == 2:
                    if counter == 1:
                        ax4[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),O2_1[0:99998],label = 'b$\mathregular{_{H}}$, f_$\mathregular{_{volc}}$ = ' + str(f_volc_pick/10), linewidth=15,color = 'black',linestyle = ln_str,alpha=0.5)
                        ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),H2O_1[0:99998],linewidth=15,color = 'black',linestyle = ln_str,alpha=0.5)
                        ax4[2].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),np.nan_to_num(CO_1[0:99998]),linewidth=15,color = 'black',linestyle = ln_str,alpha=0.5)
                
                if counter == 2:
                    if option == 1:
                        ax4[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),100*(O2_2[0:99998] - O2_1[0:99998])/O2_1[0:99998],label = '1.3b$\mathregular{_{H}}$, f_$\mathregular{_{volc}}$ = ' + str(f_volc_pick/10), linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                        ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),100*(H2O_2[0:99998] - H2O_1[0:99998])/H2O_1[0:99998],linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                        ax4[2].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),100*(CO_2[0:99998] - CO_1[0:99998])/CO_1[0:99998],linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                    if option == 2:
                        ax4[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),O2_2[0:99998],linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str) # 
                        ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),H2O_2[0:99998],linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                        ax4[2].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),np.nan_to_num(CO_2[0:99998]),linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                
                elif counter == 3:
                    if option == 1:
                        ax4[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),100*(O2_3[0:99998] - O2_1[0:99998])/O2_1[0:99998],label = '0.7b$\mathregular{_{H}}$, f_$\mathregular{_{volc}}$ = ' + str(f_volc_pick/10),linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                        ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),100*(H2O_3[0:99998] - H2O_1[0:99998])/H2O_1[0:99998],linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                        ax4[2].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),100*(CO_3[0:99998] - CO_1[0:99998])/CO_1[0:99998],linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                    if option == 2:
                        ax4[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),O2_3[0:99998],linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str) #
                        ax4[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),H2O_3[0:99998],linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                        ax4[2].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),np.nan_to_num(CO_3[0:i1]),linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                    
                ax4[0].legend(fontsize=100,loc='best',bbox_to_anchor=(-0.5, 1.0))
                
                
                if melt == 0:
                    ax4[0].set_title('$\mathbf{No \; runaway \; greenhouse \; surface \; melting - }$' + r'$\mathbf{' + str(float( '%g' % (gwater))) + '}$' + ' $\mathbf{m \;GEL,}$ ' + '\n'  + str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O, $\mathregular{f_{O_{2}}}$=FMQ' + str(float('%.1g' % (FMQ))) + ',\n'' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = 0.1-9', fontsize=100)
                    
if plot_number == 6:
    
    option = 2
    
    search_str = str("*tstart15_gwater300_H2O5_CO2500_f_ext3_*_FMQ-3.csv")
    
    max_O2_modern = 50.e-6*93.e5/(g*Rmm_CO2) 
    max_H2O_modern = 100.e-6*93.e5/(g*Rmm_CO2) 
    max_CO_modern = 50.e-6*93.e5/(g*Rmm_CO2) 
    
    f_volc_var = np.array([1,5,9,10])
    
    fig4, ax4 = plt.subplots(3,figsize=(50,105))

    ax4[0].invert_xaxis()
    ax4[2].set_xlim((-0.2,4.2))

    ax4[1].invert_xaxis()
    ax4[2].invert_xaxis()

    if option == 1:
        
        ax4[0].set_ylabel('% Difference in  [O$\mathregular{_{2}}$]$\mathregular{_{atm}}$', fontsize=100,weight='bold')
        ax4[0].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[0].tick_params(axis='both', which='major', labelsize=100) 
        
        ax4[1].set_ylabel('% Difference in  [H$\mathregular{_{2}}$O]$\mathregular{_{atm}}$', fontsize=100,weight='bold')
        ax4[1].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[1].tick_params(axis='both', which='major', labelsize=100)
        
        ax4[2].set_ylabel('% Difference in [CO]$\mathregular{_{atm}}$', fontsize=100,weight='bold')
        ax4[2].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[2].tick_params(axis='both', which='major', labelsize=100)
    
    if option == 2:
        
        ax4[0].set_ylabel('[O$\mathregular{_{2}}$]$\mathregular{_{model}}$ / [O$\mathregular{_{2}}$]$\mathregular{_{modern}}$', fontsize=100,weight='bold')
        ax4[0].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[0].tick_params(axis='both', which='major', labelsize=100)
        
        ax4[1].set_ylabel('[H$\mathregular{_{2}}$O]$\mathregular{_{model}}$ / [H$\mathregular{_{2}}$O]$\mathregular{_{modern}}$ ', fontsize=100,weight='bold')
        ax4[1].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[1].tick_params(axis='both', which='major', labelsize=100)
        
        ax4[2].set_ylabel('[CO]$\mathregular{_{model}}$ / [CO]$\mathregular{_{modern}}$', fontsize=100,weight='bold')
        ax4[2].set_xlabel('Time before present (Gyr)', fontsize=100,weight='bold')
        ax4[2].tick_params(axis='both', which='major', labelsize=100)
        
        ax4[0].axhline(y=1,color ='black', linestyle='--',linewidth=15)
        ax4[1].axhline(y=1,color ='black', linestyle='--',linewidth=15)
        ax4[2].axhline(y=1,color ='black', linestyle='--',linewidth=15)
        
        ax4[0].set_yscale("log")
        ax4[1].set_yscale("log")
        ax4[2].set_yscale("log")
        ax4[2].set_ylim((1e-5,1e2))

    if option == 1:
        
        ax4[0].text(4.1,0.05,'a.', fontsize=120,weight='bold')
        ax4[1].text(4.1,8,'b.', fontsize=120,weight='bold')
        ax4[2].text(4.1,0.0025,'c.', fontsize=120,weight='bold')
       
    if option == 2:
        
        ax4[0].text(4.1,10**4,'d.', fontsize=120,weight='bold')
        ax4[1].text(4.1,10**4.1,'e.', fontsize=120,weight='bold')
        ax4[2].text(4.1,10**1,'f.', fontsize=120,weight='bold')
    
    
    
    for i_n in range(0,np.size(f_volc_var)):
        
        f_volc_pick = f_volc_var[i_n]
        
        search_str = str("*tstart5_gwater500_H2O5_CO2500_f_ext3_f_volc" + str(f_volc_pick) + "_FMQ-3.csv")
        
        file_list1 = glob.glob("C:\\Users\\sasha\\Box\\Clean Code\\results\\sensitivity\\t_atm\\200_res\\" + search_str)
        file_list2 = glob.glob("C:\\Users\\sasha\\Box\\Clean Code\\results\\sensitivity\\t_atm\\215_res\\" + search_str)
        file_list3 = glob.glob("C:\\Users\\sasha\\Box\\Clean Code\\results\\sensitivity\\t_atm\\170_res\\" + search_str)
        # file_list3 = glob.glob("C:\\Users\\sasha\\Box\\Clean Code\\results\\final_revised\\nomelt_gwat300_t30\\" + search_str)
        file_list = file_list1 + file_list2 + file_list3
        
        if f_volc_pick == 1:
            ln_str = '-'
        if f_volc_pick == 5:
            ln_str = '--'
        if f_volc_pick == 9:
            ln_str = ':'
        if f_volc_pick == 10:
            ln_str = '-.'
        
        counter = 0
        
        print(np.size(file_list),file_list)
        
        if np.size(file_list) == 3:
            for x in file_list:
                counter = counter + 1
                print('counter',counter)

                f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27),header=None)
                
            
                results = np.nan_to_num(f1.to_numpy())
                    
                i = results[:,0]
                i = i.astype(int)
                i1 = i[-1]
                print(x,i1)
                
                final_tstep = results[:,0]
                
                tstart = results[0,1]
                # t_arr = np.unique(tstart)
                
                t = results[:,2]
                
                f_ext = results[0,3]
                # fext_arr = np.unique(f_ext)
                
                f_volc = results[0,4]
                # fvolc_arr = np.unique(f_volc)
                
                z_erupt = results[:,6]
                
                gwater = results[0,8]
                gwat_arr = np.unique(gwater)
                
                mCO2tot = results[0,9]
                CO2_arr = np.unique(mCO2tot)
                CO2_plotarr = np.array([300,500,1000,2000])
                
                mH2Otot = results[0,10]
                if mH2Otot > 0.6e-2:
                    mH2Otot = np.round(results[0,10],3)
                
                FMQ = results[0,11]
                
                
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
                z_erupt = ((moles_CO2*Rmm_CO2)-((1-f_volc)*(90e5/g)))/(1000*g*rho*mCO2tot)
                moles_CH4 = results[:,27]
                # bmulti = results[0,28]
                # success_Ar = results[:,28]
                
                print('final O',moles_O2[-1])
                
                if counter == 1:
                    if option == 2:
                        
                        O2_1 = moles_O2/max_O2_modern
                        H2O_1 = moles_H2O/max_H2O_modern
                        CO_1 = moles_CO/max_CO_modern
                        
                    if option == 1:
                        O2_1 = moles_O2/moles_atm
                        H2O_1 = moles_H2O/moles_atm
                        CO_1 = np.nan_to_num(moles_CO/moles_atm)
                elif counter == 2:
                    if option == 2:
                        O2_2 = moles_O2/max_O2_modern
                        H2O_2 = moles_H2O/max_H2O_modern
                        CO_2 = moles_CO/max_CO_modern
                    if option == 1:
                        O2_2 = moles_O2/moles_atm
                        H2O_2 = moles_H2O/moles_atm
                        CO_2 = np.nan_to_num(moles_CO/moles_atm)
                elif counter == 3:
                    if option == 2:
                        O2_3 = moles_O2/max_O2_modern
                        H2O_3 = moles_H2O/max_H2O_modern
                        CO_3 = moles_CO/max_CO_modern
                    if option == 1:
                        O2_3 = moles_O2/moles_atm
                        H2O_3 = moles_H2O/moles_atm
                        CO_3 = np.nan_to_num(moles_CO/moles_atm)
                
                capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)
                
                
                z_erupt = (((moles_CO2 - np.append(0,moles_CO2[0:-1]))*Rmm_CO2)/(g*rho*mCO2tot))/((t - np.append(0,t[0:-1]))/secingyr)
                CO2_mass_er = np.cumsum(z_erupt)*(t[1]-t[0])*rho*mCO2tot*g/secingyr
                z_erupt[T_surf>0]==0
                
                if option == 2:
                    if counter == 1:
                        ax4[0].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),O2_1[1:i1],label = '200K, f_$\mathregular{_{volc}}$ = ' + str(f_volc_pick/10), linewidth=15,color = 'black',linestyle = ln_str,alpha=0.5)
                        ax4[1].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),H2O_1[1:i1],linewidth=15,color = 'black',linestyle = ln_str,alpha=0.5)
                        ax4[2].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),np.nan_to_num(CO_1[1:i1]),linewidth=15,color = 'black',linestyle = ln_str,alpha=0.5)
                
                if counter == 2:
                    if option == 1:
                        ax4[0].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),100*(O2_2[1:i1] - O2_1[1:i1])/O2_1[1:i1],label = '210K, f_$\mathregular{_{volc}}$ = ' + str(f_volc_pick/10), linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                        ax4[1].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),100*(H2O_2[1:i1] - H2O_1[1:i1])/H2O_1[1:i1],linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                        ax4[2].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),100*(CO_2[1:i1] - CO_1[1:i1])/CO_1[1:i1],linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                    if option == 2:
                        ax4[0].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),O2_2[1:i1],linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str) # 
                        ax4[1].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),H2O_2[1:i1],linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                        ax4[2].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),np.nan_to_num(CO_2[1:i1]),linewidth=15,color = 'red',alpha=0.5,linestyle = ln_str)
                
                elif counter == 3:
                    if option == 1:
                        ax4[0].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),100*(O2_3[1:i1] - O2_1[1:i1])/O2_1[1:i1],label = '170K, f_$\mathregular{_{volc}}$ = ' + str(f_volc_pick/10),linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                        ax4[1].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),100*(H2O_3[1:i1] - H2O_1[1:i1])/H2O_1[1:i1],linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                        ax4[2].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),100*(CO_3[1:i1] - CO_1[1:i1])/CO_1[1:i1],linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                    if option == 2:
                        ax4[0].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),O2_3[1:i1],linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str) #
                        ax4[1].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),H2O_3[1:i1],linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                        ax4[2].plot(4.5-(t[1:i1]/(1.e9*365.*60.*60.*24.)),np.nan_to_num(CO_3[1:i1]),linewidth=15,color = 'blue',alpha=0.5,linestyle = ln_str)
                    
                ax4[0].legend(fontsize=100,loc='best',bbox_to_anchor=(-0.5, 1.0))
                
                
                if melt == 0:
                    ax4[0].set_title('$\mathbf{No \; runaway \; greenhouse \; surface \; melting - }$' + r'$\mathbf{' + str(float( '%g' % (gwater))) + '}$' + ' $\mathbf{m \;GEL,}$ ' + '\n'  + str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O, $\mathregular{f_{O_{2}}}$=FMQ' + str(float('%.1g' % (FMQ))) + ',\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$=0.1-1', fontsize=100) 