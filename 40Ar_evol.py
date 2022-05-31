# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 10:49:44 2022

@author: sasha
"""

import numpy as np
import math
import pandas as pd
from matplotlib import pyplot as plt


results = np.zeros((1,13))
nsteps = 10000 #number of timesteps in model

plotting = 0 #if 1- plot figs, if 0 - don't. 
#DO NOT SET plotting = 1 IF EXPLORING LARGE PARAMETER SPACE. EARTH WILL ENTER THE VENUS ZONE BEFORE RUNS COMPLETE.

reservoir = 0 #if reservoir = 1, assume no re-mixing of main mantle and recycled lithosphere
exp = 1 #1 if calculating 40Ar for POST habitable era evoltion, 0 if calculating 40Ar during HABITABLE ERA 

secingyr = 1.e9 * (365.*60.*60.*24.) #seconds in 1 Gyr

#%% Model constants

# Atmospheric species
Rmm_CO2 = 0.04401 #kg mol^-1
Rmm_H2O = 0.01801528 #kg mol^-1
Rmm_H = 0.001 #kg mol^-1 
Rmm_O2 = 0.032 #kg mol^-1 
Rmm_O = Rmm_O2/2
Avo = 6.023e23 #Avogadro's number
R = 8.31 #universal gas constant 
pres_Ar_max = (1.61+0.54)*1e16 #upper limit on modern mass 40Ar in Venus' atmosphere (Kaula 1999, O'Rourke & Korenaga 2015)


# Venus parameters
age_sun = 4.5 #Gyr
a_venus = 108.21e6 * 1.e3 #semimajor axis m
AU = 0.7 
g = 8.87 #m s^-2
G = 6.67e-11 #gravitational constant
kb = 1.38e-23 #boltzmann constant
A_venus = 4.6e14 #m2
R_venus = 6.0518e6 #m
R_core = 3.331e6 #m (Spohn 1991 Icarus, implemented in O'Rourke & Korenaga 2015)
M_venus = 4.867e24 #kg
rho_melt = 2800 #kg m^-3 density of basalt
rho_mantle = 4000 #kg m^-3
M_mantle = rho_mantle*((4/3)*math.pi)*(R_venus**3 - R_core**3)

lambda_Ar = 0.0581 #Gyr ^ -1, decay constant for 40K -> 40Ar
lambda_Ca = 0.4962 #Gyr ^ -1, decay constant for 40K -> Ca
lambda_K = lambda_Ar + lambda_Ca #combined decay constant
 
D_K = 1. #partition coefficient of K into melt (O'Rourke & Korenaga 2015)
dzdt_pres = 1e9 #modern Venus volcanism rate estimate
P_final = 90e5 #Pa, present day Venus CO2 pressure
b_err = 1 #e-folding timescale for volcanism in Gyr (as for oxygen loss model)
rad_K = 1.165e-4 # 40K/K - O'Rourke & Korenaga (2015)

#%% #Setup arrays for gridsearch model, uncomment commented lines to perform run for single chosen parameter combination

d_crust_arr = 1e3*np.array([10,25,50,75,100]) #m, thickness of crust
# d_crust_arr = 1e3*np.array([50])

frac_pre_arr = np.array([0.0,0.1,0.3,0.5,1.0]) #fraction of radiogenic 40Ar in mantle degassed BEFORE 40Ar model begins 
# frac_pre_arr = np.array([0.0])

CO2_arr = np.array([300,500,1000,2000])*1e-6 #concentration of CO2 in average Venus melts
# CO2_arr = np.array([2000e-6])
label_CO2 = np.array(['300','500','1000','2000']) #labels for plots 
# label_CO2 = np.array(['2000'])

frac_volc_arr = np.array([0.1,0.5,0.9,1.0]) #fraction of Venus' modern atmosphere degassed through post-habitable era magmatism
# frac_volc_arr = np.array([0.5])

frac_ext_arr = np.array([0.1,0.3,0.5,1.0]) #extrusive volcanism fraction
# frac_ext_arr = np.array([1.0])

K_U_arr = [7220, 13800] #Ratio of K to U in primordial mantle; Venus: 7220 \pm 1220 (Kaula 1999), Earth: 13,800 \pm 1300 (Arevalo et al. 2009) - Overall - O'Rourke & Korenaga (2015)
# K_U_arr = [7220]

conc_U_arr = np.array([21,17,13])*1e-9 #Primordial mantle U concentration, Kaula 1999
# conc_U_arr = np.array([21])*1e-9 

phi_arr = np.array([0.025,0.05,0.1]) #volume averaged melt fraction for producing surface melts
# phi_arr = np.array([0.05])  

#%% Set up model times for post-hab and hab-era evolutions. Times all in Gyr since Venus formation. (4.5-t gives time in Ga)
if exp == 1:
    tstart_arr = np.array([0.5,1.5,3.0])
    # tstart_arr = np.array([3.0])
    # tstart_arr = np.array([0.5])
    age_sun_arr = np.array([4.5])
elif exp == 0:
    tstart_arr = np.array([0.1])
    age_sun_arr = np.array([0.5,1.5,3.0])
    # age_sun_arr = np.array([0.5])



#%% Set up plots, if plotting
if plotting == 1:
    
    cmap = plt.get_cmap('viridis')
    colors = cmap(np.linspace(0, 1, len(CO2_arr)))
    
    fig0,ax0 = plt.subplots()
    ax0.set_ylim((14,17.5))
    sty = np.array(['-',':'])
    
    fig1, ax1 = plt.subplots(2,1,figsize=(50,40))
    ax1[1].axhline(y=(1.61+0.54)*1e16,linestyle = '--',linewidth=15)
    ax1[1].set_yscale('log')
    ax1[1].set_ylim((1e14,1e17))
    ax1[0].tick_params(axis='both', which='major', labelsize=100)
    ax1[1].tick_params(axis='both', which='major', labelsize=100)
    
    fig2, ax2 = plt.subplots(2,1,figsize=(50,40))

    ax2[0].set_yscale('log')
    ax2[1].set_yscale('log')
    ax2[0].tick_params(axis='both', which='major', labelsize=100)
    ax1[1].tick_params(axis='both', which='major', labelsize=100)

    fig3,ax3 = plt.subplots(figsize=(50,40))
    ax3.tick_params(axis='both', which='major', labelsize=100)
    
    fig4,ax4 = plt.subplots(figsize=(50,40))
    ax4.set_yscale('log')
    ax4.tick_params(axis='both', which='major', labelsize=100)
    
    fig5,ax5= plt.subplots(figsize=(50,40))
    ax5.set_yscale('log')
    ax3.tick_params(axis='both', which='major', labelsize=100)


#%% Grid search and timestepping
for i11 in range(0,np.size(K_U_arr)):
    for i10 in range(0,np.size(conc_U_arr)):
        for i9 in range(0,np.size(phi_arr)):
            for i8 in range(0,np.size(tstart_arr)):
                
                for i7 in range(0,np.size(age_sun_arr)):
                    for i6 in range(0,np.size(d_crust_arr)):
                        for i5 in range(0,np.size(frac_pre_arr)):
                            for i4 in range(0,np.size(frac_volc_arr)):
                                for i3 in range(0,np.size(frac_ext_arr)):
                                    for i2 in range(0,np.size(CO2_arr)):
                                        mCO2tot = CO2_arr[i2]
                                        conc_U = conc_U_arr[i10]
                                        phi_bar = phi_arr[i9]
                                        tstart = tstart_arr[i8]
                                        age_sun = age_sun_arr[i7]
                                        d_crust = d_crust_arr[i6]
                                        frac_pre = frac_pre_arr[i5]
                                        frac_volc = frac_volc_arr[i4]
                                        frac_ext = frac_ext_arr[i3]
                                        K_U = K_U_arr[i11]
                                        
                                        t_arr = np.linspace(tstart,age_sun,num=nsteps) #model time in Gyr
                                        t = t_arr*secingyr #model time in seconds
                                        
                                        #Set up arrays to record model output
                                        K_man = np.zeros(nsteps) #Mantle 40K concentration
                                        K_to_m = np.zeros(nsteps) #Mass 40K returned to mantle at each timestep through recycling at base of crust
                                        K_to_l = np.zeros(nsteps) #Mass 40K removed to crust at each timestep via mantle melting
                                        #Intitial mantle 40K concentration
                                        K_man[0] = (conc_U*K_U*rad_K)*math.exp(lambda_K*(4.5 - tstart))
                                        conc_0 = K_man[0]
                                        abs_0 = (conc_U*K_U*rad_K)*math.exp(lambda_K*(4.5))
                                        
                                        crust_tot = np.zeros(nsteps)
                                        
                                        Ar_atm = np.zeros(nsteps)
                                        Ar_man = np.zeros(nsteps)
                                        Ar_crust = np.zeros(nsteps)
                                        Ar_melt = np.zeros(nsteps)
                                        Ar_man[0] = abs_0*(1 - math.exp(-lambda_Ar*tstart))*(1-frac_pre)
                                        Ar_atm[0] = abs_0*(1 - math.exp(-lambda_Ar*tstart))*(frac_pre)*M_mantle
                                        
                                        #Set up depth array for tracking 40K concentration in crust
                                        dz = np.arange(0,d_crust,1)
                                        K_dz = np.zeros(np.size(dz)) #Start with no K in lithosphere - unrealistic, but allows quantification of post-hab era only
                                        
                                        #Calculate total crustal production needed to degas up to Venus' modern atmospheric pressure during model
                                        M_CO20 = P_final*(1-frac_volc)/g
                                        M_CO21 = P_final/g
                                        z_lava = (M_CO21 - M_CO20)/(mCO2tot*rho_melt)
                                        
                                        if exp == 1:
                                            if frac_ext < 1:
                                                dzt_p = (1.e9*dzdt_pres/A_venus)/frac_ext #m per Gyr erupted per unit surface area
                                            elif frac_ext == 1:
                                                dzt_p = (1.e9*dzdt_pres/A_venus)
                                            
                                            #Parameters for exponential function to calculate crustal production rate dz(t) = m*exp(-b*t) + nn
                                            nn = 0.
                                            m = (dzt_p - (z_lava*-b_err))/np.exp(-b_err*tstart)
                                            
                                            #Calculate melt production rate over time
                                            MP = m*np.exp(-b_err*(t/secingyr))
                                            
                                        if exp == 0: #Assume constant crustal production rate over course of habitable era
                                            MP = np.ones(nsteps)*z_lava/(age_sun - tstart)
                                            
                                        for i1 in range(1,nsteps):
                                            
                                            dt = (t_arr[i1]-t_arr[i1-1]) #timestep size
                                            
                                            K_man[i1] = K_man[i1-1] - K_man[i1-1]*(1-np.exp(-dt*lambda_K)) #radioactive decay of 40K in mantle
                                            
                                            Ar_man[i1] = Ar_man[i1-1] + (lambda_Ar/lambda_K)*K_man[i1-1]*(1-np.exp(-dt*lambda_K)) #Accumulation of 40Ar due to 40K decay
                                            
                                            Ar_crust[i1] = rho_melt * A_venus * np.trapz(K_dz,dz) *(lambda_Ar/lambda_K)*(1-np.exp(-dt*lambda_K)) #Ar added to atmosphere by 40K decay in crust
                                            
                                            K_dz = K_dz*math.exp(-lambda_K*dt) #radioactive decay of 40K in crust at each depth
                                            
                                            z_new = MP[i1]*dt #added thickness of crust (assume all at surface)
                                            crust_tot[i1] = crust_tot[i1-1] + z_new #total thickness of crust
                                            
                                            Ar_melt[i1] = z_new * rho_melt * A_venus * Ar_man[i1]/phi_bar #Total mass of Ar added to atmosphere by mantle melting (assume degassed instantaneously; O'Rourke & Korenaga 2015)
                                            
                                            if Ar_melt[i1] > M_mantle * Ar_man[i1]:
                                                Ar_melt[i1] = M_mantle * Ar_man[i1]
                                                
                                            dz = dz + z_new #update depth array
                                            K_to_l[i1] = z_new * rho_melt * (K_man[i1]/phi_bar) #mass 40K per column removed from mantle by eruptions
                                            
                                            
                                            if z_new >= 1.:
                                                #Append points to top of crust 40K concentration array to record added 40K bearing crust
                                                new_crust = np.arange(0,z_new,1)
                                                dz = np.append(new_crust,dz)
                                                
                                                new_K = np.ones(np.size(new_crust)) * K_to_l[i1]/(z_new*rho_melt) #concentration of 40K in new crust
                                                K_dz = np.append(new_K,K_dz)
                                                
                                                #Remove crust deeper than set crustal thickness from base of crust and add 40K back to mantle
                                                to_del = np.where(dz>=d_crust) #find points beneath thickness of lithosphere (fixed) - assume recycled
                                                
                                                K_to_m[i1] = np.mean(K_dz[to_del]) * z_new * rho_melt #mass 40K per column returned to mantle by recycling
                                                
                                                #Ignore re-mixing of 40K into mantle if assuming recycled crust enters separate mantle reservoir
                                                if reservoir == 1:
                                                    K_to_m[i1] = 0.
                                                
                                                dz = np.delete(dz,to_del)
                                                K_dz = np.delete(K_dz,to_del)
                                            
                                            #Update all 40K and 40Ar inventories at end of timestep
                                            K_man[i1] = K_man[i1] + (K_to_m[i1] - K_to_l[i1])*A_venus/M_mantle #assumes instantaneous mixing throughout mantle 
                                            Ar_man[i1] = Ar_man[i1] - Ar_melt[i1]/M_mantle
                                            Ar_atm[i1] = Ar_atm[i1-1] + Ar_melt[i1] + Ar_crust[i1]
                                        
                                        #At end of run, determine whether final mass of 40Ar in atmosphere <modern upper limit.
                                        #If yes, mark run as "success"
                                        if Ar_atm[-1]>pres_Ar_max:
                                            success = 0
                                        elif Ar_atm[-1]<=pres_Ar_max:
                                            success = 1
                                        
                                        #Save parameter combination and results
                                        add_row = np.array([tstart,age_sun,exp,phi_bar,mCO2tot,frac_pre,frac_volc,frac_ext,conc_U,K_U,d_crust,Ar_atm[-1],success])
                                        results = np.vstack((results,add_row))
                                    #%%
                                        if plotting == 1:
                                            ax1[1].plot(t_arr,Ar_atm,c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax1[0].plot(t_arr,MP/1e3,c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax1[0].set_title('$\mathregular{z_{lith}}$ = ' + str(float('%g' % (d_crust/1e3))) 
                                                             + ' km, $\mathregular{f_{melt}}$ = '+ str(float('%g' % (phi_bar))) 
                                                             + ', $\mathregular{f_{ext}}$ = '+ str(float('%g' % (frac_ext)))
                                                             + ', $\mathregular{f_{volc}}$ = '+ str(float('%g' % (frac_volc)))
                                                             + ', $\mathregular{f_{pre}}$ = '+ str(float('%g' % (frac_pre)))
                                                             , fontsize=100) #', $\mathregular{t_{start}}$ = '+ str(float('%g' % (4.5-tstart))) + 'Ga')
                                            ax1[0].legend(fontsize=100,loc='center right')
                                            
                                            
                                            ax1[0].set_ylabel('Melt production' '\n' 'rate (km/Gyr)', fontsize=100,weight='bold')
                                            ax1[1].set_ylabel('Atmospheric '' \n' '$\mathregular{^{40}}$Ar (kg)', fontsize=100,weight='bold')
                                            ax1[1].set_xlabel('Time (Gyr)', fontsize=100,weight='bold')
                                        
                                            ax2[0].plot(t_arr,K_to_l,c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax2[0].plot(t_arr,K_to_m,linestyle=':',c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax2[1].set_ylabel('$\mathregular{^{40}}$K (kg)')
                                            ax2[0].legend(fontsize=100,loc='center right')
                                            ax2[0].set_title('$\mathregular{z_{lith}}$ = ' + str(float('%g' % (d_crust/1e3))) 
                                                             + ' km, $\mathregular{f_{melt}}$ = '+ str(float('%g' % (phi_bar))) 
                                                             + ', $\mathregular{f_{ext}}$ = '+ str(float('%g' % (frac_ext)))
                                                             + ', $\mathregular{f_{volc}}$ = '+ str(float('%g' % (frac_volc)))
                                                             + ', $\mathregular{f_{pre}}$ = '+ str(float('%g' % (frac_pre)))
                                                             , fontsize=100) #', $\mathregular{t_{start}}$ = '+ str(float('%g' % (4.5-tstart))) + 'Ga')
                                            
                                            ax2[1].plot(t_arr,Ar_melt,c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax2[1].plot(t_arr,Ar_crust,linestyle=':',c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax2[1].set_ylabel('$\mathregular{^{40}}$Ar (kg)', fontsize=100,weight='bold')
                                            ax2[1].set_xlabel('Time (Gyr)', fontsize=100,weight='bold')
                                            
                                            ax2[0].set_ylabel('K (kg)')
                                            ax2[1].set_ylabel('Ar (kg)')
                                            ax2[1].set_xlabel('Time (Gyr)')
                                            
                                            ax2[0].set_title('$\mathregular{z_{lith}}$ = ' + str(float('%g' % (d_crust/1e3))) + ' km, $\mathregular{f_{melt}}$ = '+ str(float('%g' % (phi_bar)))) #', $\mathregular{t_{start}}$ = '+ str(float('%g' % (4.5-tstart))) + 'Ga')
                                            
                                            ax3.plot(t_arr,crust_tot/1e3,c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax3.set_xlabel('Time (Gyr)', fontsize=100,weight='bold')
                                            ax3.set_ylabel('Total crustal production (km)', fontsize=100,weight='bold')
                                            ax3.legend(fontsize=100,loc='center right')
                                            ax3.set_title('$\mathregular{z_{lith}}$ = ' + str(float('%g' % (d_crust/1e3))) 
                                                             + ' km, $\mathregular{f_{melt}}$ = '+ str(float('%g' % (phi_bar))) 
                                                             + ', $\mathregular{f_{ext}}$ = '+ str(float('%g' % (frac_ext)))
                                                             + ', $\mathregular{f_{volc}}$ = '+ str(float('%g' % (frac_volc)))
                                                             + ', $\mathregular{f_{pre}}$ = '+ str(float('%g' % (frac_pre)))
                                                             , fontsize=100) #', $\mathregular{t_{start}}$ = '+ str(float('%g' % (4.5-tstart))) + 'Ga')
                                            
                                            ax4.plot(t_arr,K_man,c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax4.plot(t_arr,Ar_man,linestyle = ':',c=colors[i2],label=label_CO2[i2],linewidth=15)
                                            ax4.set_xlabel('Time (Gyr)', fontsize=100,weight='bold')
                                            ax4.set_ylabel('Mantle concentration', fontsize=100,weight='bold')
                                            ax4.legend(fontsize=100,loc='center right')
                                            ax4.set_title('$\mathregular{z_{lith}}$ = ' + str(float('%g' % (d_crust/1e3))) 
                                                             + ' km, $\mathregular{f_{melt}}$ = '+ str(float('%g' % (phi_bar))) 
                                                             + ', $\mathregular{f_{ext}}$ = '+ str(float('%g' % (frac_ext)))
                                                             + ', $\mathregular{f_{volc}}$ = '+ str(float('%g' % (frac_volc)))
                                                             + ', $\mathregular{f_{pre}}$ = '+ str(float('%g' % (frac_pre)))
                                                             , fontsize=100) #', $\mathregular{t_{start}}$ = '+ str(float('%g' % (4.5-tstart))) + 'Ga')
                                        
# %% Save output to .csv
data = {'tstart':results[:,0],
        'tend':results[:,1],
        'exp':results[:,2],
        'phi':results[:,3],
        'CO2':results[:,4],
        'f_pre':results[:,5],
        'f_volc':results[:,6],
        'f_ext':results[:,7],
        'conc_U':results[:,8],
        'K_U':results[:,9],
        'd_crust':results[:,10],
        '40_Ar':results[:,11],
        'success':results[:,12]}
        
df = pd.DataFrame(data)
df.drop_duplicates(keep=False,inplace=True)
# print(df)

if reservoir == 1:
    
    if exp == 1:
        title_str = '40_Ar_exp_dlith_KU_nomix.csv'
    elif exp == 0:
        title_str = '40_Ar_hab_dlith_KU_nomix.csv'
        
elif reservoir == 0:
    
    if exp == 1:
        title_str = '40_Ar_exp_dlith_KU_mix.csv'
    elif exp == 0:
        title_str = '40_Ar_hab_dlith_KU_mix.csv'


df.to_csv(title_str,header = False) 