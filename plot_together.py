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
import glob2 as glob
import seaborn

#%%
#Atmospheric species
Rmm_CO2 = 0.04401 #kg mol^-1
Rmm_H2O = 0.01801528 #kg mol^-1
Rmm_H = 0.001 #kg mol^-1 
Rmm_O2 = 0.032 #kg mol^-1 
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

gwat_arr = np.array([100,300,500,700,1000])

CO2_arr = np.array([0.0003,0.0005,0.001,0.002])

H2O_arr = 1e-2*np.array([0.001,0.1,0.2,0.5,0.7,1.])

moles_atm_old = 0
Rmm_atm_old = 0
moles_CO2_old = 0
moles_H2O_old = 0
moles_O2_old = 0
#%%
plot_number = 1
# file_list = glob.glob("C:\\Users\\sasha\\Box\\third_year\\venus_earlyhab\\volccalc\\VolcGases-master\\results\\repaired_all\\old\\tstart5*gwater500*")
file_list = glob.glob("C:\\Users\\sasha\\Box\\Clean Code\\test*")

# file_list = glob.glob("C:\\Users\\sasha\\Box\\third_year\\venus_earlyhab\\volccalc\\VolcGases-master\\results\\no_melting\\old\\NM_tstart5*gwater300_*")

counter = 0
i2 = 0

#%%
if plot_number == 1:
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
        
        f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),header=None)
        
        results = f1.to_numpy()
        
        i = results[:,0]
        i = i.astype(int)
        i1 = i[-1]
        
        tstart = results[0,1]
        t = secingyr*np.linspace(tstart,4.5,num=(np.max(i1)+1),endpoint=True)
        f_ext = results[0,3]
        f_volc = results[0,4]
        
        # print('ext,volc',f_ext,f_volc)
        b_err = results[0,5]
        z_erupt = results[:,6]
        q = results[0,7]
        gwater = results[0,8]
        mCO2tot = results[0,9]
        mH2Otot = results[0,10]
        FMQ = results[0,11]
        Pres = results[:,12]
        moles_atm = results[:,13]
        moles_CO2 = results[:,14]
        moles_H2O = results[:,15]
        molesH2O_loss = results[:,16]
        moles_O2 = results[:,17]
        Rmm_atm = results[:,18]
        T_surf = results[:,19]
        ash = results[:,20]
        lava = results[:,21]
        magma = results[:,22]
        nt = results[:,23]
        
        
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
        
        i2 = i2 + 1
        print('valid, i2', i2,' gwat', gwater, ' tstart', tstart)
        
        # plot_colmass, = plt.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=10,color='lightgrey',linestyle=':')
        
        # plot_colCO2, = plt.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=10,color ='bisque')
        
        ax1[0].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=15,color ='slateblue',alpha=0.05)
        
        ax1[1].plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=15,color ='lightseagreen',alpha =0.05)
        
        
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
    
    #Highlight specific plot
    # f1 = pd.read_csv("C:\\Users\\sasha\\Box\\third_year\\venus_earlyhab\\volccalc\\VolcGases-master\\results\\repaired_all\\old\\tstart5_gwater500_H2O0_CO2300_f_ext10_f_volc9.csv",
    #                   usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),header=None)
    
    # f1 = pd.read_csv("C:\\Users\\sasha\\Box\\third_year\\venus_earlyhab\\volccalc\\VolcGases-master\\results\\no_melting\\old\\NM_tstart5_gwater100_H2O0_CO2300_f_ext10_f_volc9.csv",
    #                     usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),header=None)
    f1 = pd.read_csv("C:\\Users\\sasha\\Box\\Clean Code\\test_plottstart5_gwater100_H2O0_CO2300_f_ext10_f_volc9.csv",
                        usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),header=None)
    
    
        
    results = f1.to_numpy()
    
    i = results[:,0]
    i = i.astype(int)
    i1 = i[-1]
    
    tstart = results[0,1]
    t = secingyr*np.linspace(tstart,4.5,num=(np.max(i1)+1),endpoint=True)
    f_ext = results[0,3]
    f_volc = results[0,4]
    
    # print('ext,volc',f_ext,f_volc)
    b_err = results[0,5]
    z_erupt = results[:,6]
    q = results[0,7]
    gwater = results[0,8]
    mCO2tot = results[0,9]
    mH2Otot = results[0,10]
    if mH2Otot == 0:
        mH2Otot = 0.001e-2
    FMQ = results[0,11]
    Pres = results[:,12]
    moles_CO2 = results[:,14]
    moles_H2O = results[:,15]
    moles_O2 = results[:,17]
    moles_atm = moles_CO2 + moles_H2O + moles_O2
    Rmm_atm = results[:,18]
    T_surf = results[:,19]
    ash = results[:,20]
    lava = results[:,21]
    magma = results[:,22]
    nt = results[:,23]
    
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
    
    ax1[0].text(1.4,35,'Modern Venus H$\mathregular{_{2}}$O', fontsize=90,weight='bold',color='midnightblue')
    ax1[1].text(1.3,32,'Modern Venus O$\mathregular{_{2}}$', fontsize=90,weight='bold',color='teal')
    
    ax1[0].text(4.1,1e6,'i. H$\mathregular{_{2}}$O', fontsize=120,weight='bold')
    ax1[1].text(4.1,1e6,'ii. O$\mathregular{_{2}}$', fontsize=120,weight='bold')
    # plt.axhline(y=1, color ='lightgrey', linestyle='--',linewidth=2)
    
    # plt.xlim(2,2.1)
    # plt.legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=40)
    ax1[0].legend(handles=[plot_colmass,plot_colCO2,plot_colH2O], fontsize=100,loc='center right')
    ax1[1].legend(handles=[plot_colmass,plot_colCO2,plot_colO2], fontsize=100,loc='center right')
    # ax1[1].legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=100)
    ax1[0].set_title('$\mathbf{a. No \; runaway \; greenhouse \; surface \; melting}$' '\n' + str(float( '%g' % (gwater))) + ' m GEL, '+ str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))), fontsize=110) 
    # ax1[0].set_title('$\mathbf{b. With\;runaway\;greenhouse\;surface\;melting}$' '\n' + str(float( '%g' % (gwater))) + ' m GEL, '+ str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))), fontsize=110) 
    
    
    ax1[0].annotate('', xy=(3.5,10),  xycoords='data',
            xytext=(3.5, 25), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[0].annotate('', xy=(2.5,10),  xycoords='data',
            xytext=(2.5, 25), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[0].annotate('', xy=(1.5,10),  xycoords='data',
            xytext=(1.5, 25), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[0].annotate('', xy=(0.5,10),  xycoords='data',
            xytext=(0.5, 25), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    
    ax1[1].annotate('', xy=(3.5,10),  xycoords='data',
            xytext=(3.5, 22), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[1].annotate('', xy=(2.5,10),  xycoords='data',
            xytext=(2.5, 22), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[1].annotate('', xy=(1.5,10),  xycoords='data',
            xytext=(1.5, 22), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    ax1[1].annotate('', xy=(0.5,10),  xycoords='data',
            xytext=(0.5, 22), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
    
#%%

if plot_number == 2:
    
    file_list = glob.glob("C:\\Users\\sasha\\Box\\third_year\\venus_earlyhab\\volccalc\\VolcGases-master\\results\\no_melting\\old\\NM_tstart*gwater300_H2O0_CO2300_f_ext10_f_volc10.csv")
    
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
        
        f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),header=None)
        
        results = f1.to_numpy()
        
        i = results[:,0]
        i = i.astype(int)
        i1 = i[-1]
        
        tstart = results[0,1]
        t = secingyr*np.linspace(tstart,4.5,num=(np.max(i1)+1),endpoint=True)
        f_ext = results[0,3]
        f_volc = results[0,4]
        
        # print('ext,volc',f_ext,f_volc)
        # b_err = results[0,5]
        # z_erupt = results[:,6]
        # q = results[0,7]
        gwater = results[0,8]
        mCO2tot = results[0,9]
        mH2Otot = results[0,10]
        if mH2Otot == 0:
            mH2Otot = 0.001e-2
        # FMQ = results[0,11]
        # Pres = results[:,12]
        moles_atm = results[:,13]
        moles_CO2 = results[:,14]
        moles_H2O = results[:,15]
        molesH2O_loss = results[:,16]
        moles_O2 = results[:,17]
        Rmm_atm = results[:,18]
        T_surf = results[:,19]
        # ash = results[:,20]
        # lava = results[:,21]
        # magma = results[:,22]
        # nt = results[:,23]
        
        
        z_erupt = (((moles_CO2 - np.append(0,moles_CO2[0:-1]))*Rmm_CO2)/(g*rho*mCO2tot))/((t - np.append(0,t[0:-1]))/secingyr)
        CO2_mass_er = np.cumsum(z_erupt)*(t[1]-t[0])*rho*mCO2tot*g/secingyr
        z_erupt[T_surf>0]==0
        
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
        
        if tstart == 1.5:
            # plot_colmass, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=15,color='black',linestyle='--')
            plot_colCO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=15,color ='red',linestyle='--')
            # plot_colCO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),CO2_mass_er[0:i1],label='CO\N{SUBSCRIPT TWO}',linewidth=15,color ='red',linestyle='--')
            plot_colH2O, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=15,color ='midnightblue',linestyle='--')
            plot_colO2, = ax2.plot(4.5-(t[0:i1]/(1.e9*365.*60.*60.*24.)),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=15,color ='lightseagreen',linestyle='--')
        
        ax2.annotate('', xy=(3.5,10),  xycoords='data',
            xytext=(3.5, 22), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
        ax2.annotate('', xy=(2.5,10),  xycoords='data',
            xytext=(2.5, 25), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
        ax2.annotate('', xy=(1.5,10),  xycoords='data',
            xytext=(1.5, 22), textcoords='data',
            arrowprops=dict(facecolor='teal', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
        ax2.annotate('', xy=(-0.15,10),  xycoords='data',
            xytext=(-0.15, 25), textcoords='data',
            arrowprops=dict(facecolor='midnightblue', shrink=0,width=20,headwidth=50,headlength=60),
            horizontalalignment='left', verticalalignment='bottom',
            )
        
    #%%
    
    ax2.axhline(y=(100.e-6 * moles_atm[i1] * Rmm_H2O),color ='midnightblue', linestyle='--',linewidth=15)
    ax2.axhline(y=(50.e-6 * moles_atm[i1] * Rmm_O2), color ='teal', linestyle='--',linewidth=15)
    
    ax2.text(1.4,30,'Modern Venus H$\mathregular{_{2}}$O', fontsize=90,weight='bold',color='midnightblue')
    ax2.text(1.3,14,'Modern Venus O$\mathregular{_{2}}$', fontsize=90,weight='bold',color='teal')
    
    ax2.text(3.95,1e6,'$\mathregular{t_{0}}=4.0$', fontsize=90,color='grey')
    ax2.text(2.95,1e6,'$\mathregular{t_{0}}=3.0$', fontsize=90,color='grey')
    # plt.axhline(y=1, color ='lightgrey', linestyle='--',linewidth=2)
    
    # plt.xlim(2,2.1)
    # plt.legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=40)
    ax2.legend(handles=[plot_colCO2,plot_colH2O,plot_colO2], fontsize=100,loc='center right')
    ax2.set_title('$\mathbf{a. No \;runaway\;greenhouse\;surface\;melting}$' '\n' + str(float('%g' % (gwater))) + ' m GEL, '+ str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, '   + str(float('%.1g' % (1.e2*mH2Otot))) +  ' wt% H$\mathregular{_{2}}$O,' '\n' ' f$\mathregular{_{ext}}$=' + str(float('%.1g' % (f_ext))) + ', f$\mathregular{_{volc}}$ = '  + str(float('%.1g' % (f_volc))), fontsize=100) 

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
        
        f1 = pd.read_csv(x,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),header=None)
        
        results = f1.to_numpy()
        
        i = results[:,0]
        i = i.astype(int)
        i1 = i[-1]
        
        tstart = results[0,1]
        t = secingyr*np.linspace(tstart,4.5,num=(np.max(i1)+1),endpoint=True)
        f_ext = results[0,3]
        f_volc = results[0,4]
        
        # print('ext,volc',f_ext,f_volc)
        # b_err = results[0,5]
        z_erupt = results[:,6]
        # q = results[0,7]
        gwater = results[0,8]
        mCO2tot = results[0,9]
        mH2Otot = results[0,10]
        if mH2Otot == 0:
            mH2Otot = 0.001e-2
        # FMQ = results[0,11]
        # Pres = results[:,12]
        moles_atm = results[:,13]
        moles_CO2 = results[:,14]
        print(moles_CO2[-1])
        moles_CO2 = np.cumsum(z_erupt*((t[1]-t[0])/secingyr)*rho*mCO2tot*8.87/Rmm_CO2)
        moles_H2O = results[:,15]
        molesH2O_loss = results[:,16]
        moles_O2 = results[:,17]
        Rmm_atm = results[:,18]
        T_surf = results[:,19]
        # ash = results[:,20]
        # lava = results[:,21]
        # magma = results[:,22]
        # nt = results[:,23]
        
        
        
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