# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 09:43:16 2022

@author: sasha
"""
from scipy.io import loadmat
import numpy as np
import time
import math
import random
from VolcGases.functions import solve_gases
from scipy import optimize,interpolate
from scipy.interpolate import LinearNDInterpolator
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
import sys
# import pathlib
# from pathlib import Pathimport 
import pandas as pd
from datetime import datetime
from modular_functions_clean_constantvolc import timestep,surface_T,stdrk4
import sklearn
from sklearn.neighbors import KNeighborsRegressor,NearestNeighbors,RadiusNeighborsRegressor

#%%


results = np.zeros(13)

P_final = 90.e5

frac_volc_arr = np.array([0.1,0.5,0.9,1.0])
tstart_arr = np.array([0.5,1.5,3.])
GEL_arr = np.array([1000.])


mH2Otot = 0.1e-2
gwater_GEL = 500. #m
mCO2tot = 1000e-6
frac_volc = 1.
# frac_ext = 0.5 

FMQ = 0.
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
T_atm = 300.
dzdt_pres = 1e9


#%%
secingyr = 1.e9 * (365.*60.*60.*24.)



# Atmospheric species
Rmm_CO2 = 0.04401 #kg mol^-1
Rmm_H2O = 0.01801528 #kg mol^-1
Rmm_H = 0.001 #kg mol^-1 
Rmm_O2 = 0.032 #kg mol^-1 
Rmm_O = Rmm_O2/2
Avo = 6.023e23 #Avogadro's number
R = 8.31 #universal gas constant 

# Venus parameters
age_sun = 4.5 #Gyr
a_venus = 108.21e6 * 1.e3 #semimajor axis m
AU = 0.7 
g = 8.87 #m s^-2
G = 6.67e-11 #gravitational constant
kb = 1.38e-23 #boltzmann constant
A_venus = 4.6e14 #m2
R_venus = 6.0518e6 #m
M_venus = 4.867e24 #kg
rho = 2800 #kg m^-3 density of basalt
max_water = 0.0125 #maximum mass fraction water that can be oxidized by ash
colEO = 270.*1.e5/(9.81*Rmm_H2O) #column moles water for 1 EO

# For degassing calculations: VolcGases input parameters for FMQ = 0
A = 25738.
B = 9.
C = 0.092

#import Gillmann et al. 2020 data for nonthermal O escape
gillmann = np.genfromtxt('gillmann2020.csv',delimiter=',')
gill_t = gillmann[:,0]/1.e3
gill_molesO2 = gillmann[:,3]

#Turbet et al 2020 surface T analytical function coefficients
c1 = 3.401
c2 = 1.501e-1
c3 = -3.146e-2
c4 = 4.702e-2
c5 = -4.911e-3
c6 = 8.519e-3
c7 = -1.467e-2
c8 = -7.091e-3
c9 = -7.627e-3
c10 = 8.348e-3
k1 = 2.668
k2 = 1.019
k3 = 1.099
k4 = 4.683e-1
k5 = 7.664e-1
k6 = 4.224e-1

y_T = (math.log10(g) - k3)/k4

# import Turbet et al. 2019 for backup T_surf calculations during runaway GH
turbet = np.genfromtxt('turbet.csv',delimiter=',')
turbet_f = turbet[:,0]
turbet_w = np.log10(turbet[:,1])
turbet_T = turbet[:,2]
rbfi = LinearNDInterpolator(list(zip(turbet_f, turbet_w)), turbet_T)

# Heat diffusion problem set-up for surface melt layer during runaway GH
T_bg = 373. #K, temperature at end of habitable period
k = 4. #W m^-1 #K^-1
cp = 1.e3 #J kg^-1 K^-1
kappa = k/(rho*cp)  
lat_melt = 4e5 #J kg^-1
nz = 500

# Matrix for heat diffusion calculations during runaway greenhouse
matrix = np.zeros([nz,nz])   
for i in range(1,nz-1):
    matrix[i][i-1] = 1.
    matrix[i][i] = -2.
    matrix[i][i+1] = 1.
    # print(matrix[i][:])

matrix[nz-1][nz-2] = 1.
matrix[nz-1][nz-1] = -1.

matrix = csr_matrix(matrix)

#%%

z_end = 50.e3 # max depth of interest for surface melt layer
z = np.linspace(0,z_end,nz)
P_z = z*rho*g
dz = z[1]-z[0]

# Import basalt solidus from alphaMELTS calculations
solidus = np.genfromtxt('basalt_sol.csv',delimiter=',',skip_header=1)

# Import melt fraction data from alphaMELTS results
melting = np.genfromtxt('melt_f_old.csv',delimiter=',',skip_header=1)
melt_T = melting[:,1] #K
melt_P = melting[:,0]*1e5 #Pa
melt_f = melting[:,3] #fraction
# Find fit to melt fraction as a function of temperature from alphaMELTS results
# melt_function = interpolate.interp2d(melt_T,melt_P,melt_f,kind='linear')

melting = pd.read_csv('./melt_f_old.csv', usecols=[0,1,2,3])
melting = melting.dropna(axis=0).sort_values(['Pressure (bars)', 'Temperature (K)'])

mf_10k = interpolate.interp1d(melt_T[melt_P==(10000*1e5)],melt_f[melt_P==(10000*1e5)],kind='linear')
mf_15k = interpolate.interp1d(melt_T[melt_P==(15000*1e5)],melt_f[melt_P==(15000*1e5)],kind='linear')
arr_10_15 = np.array([10000.,15000.])*1e5

knn = KNeighborsRegressor(n_neighbors=3,weights='distance',p=2)
knn.fit(melting[['Pressure (bars)', 'Temperature (K)']].values, melting['melt fraction'].values)


#%%
luminosity = np.genfromtxt('luminosity.csv',delimiter=',')
solart = luminosity[:,0] #Gyr
solarL = luminosity[:,1]/(AU**2.) #relative to Earth present day 
L = interpolate.interp1d(solart,solarL,kind='linear')

n_run = 100
Pgrid = np.linspace(0.001,200,num=n_run)

degas = 0
    
magma_yn = 1.

# Calculate stable timestep size 
# nsteps = 5*int((16/5)*((t2 - t1)/secingyr)/(1e-6))
#%%
for i5 in range(0,np.size(GEL_arr)):
    gwater_GEL = GEL_arr[i5]
    
    
        
    
    for i4 in range(2,np.size(tstart_arr)):
        #start & end times in Gyr
        if gwater_GEL == 500:
            if i4 < 2:
                i4 = int(2)
        
        t_start = tstart_arr[i4]
        
        
        
        if gwater_GEL < 500.:
            
            nsteps = int(10000000)
            t_end = t_start+(5*nsteps*1e-9)
            
        if gwater_GEL > 300.:
            nsteps = int(100000000)
            left = 4.5 - t_start
            t_end = t_start+(left*10*nsteps*1e-9)
            

        t1 = t_start*secingyr
        t2 = t_end*secingyr
        
        
        
        for i3 in range(1,np.size(frac_volc_arr)):
            if gwater_GEL == 700:
                if i3 < 2:
                    i3 = int(2)
        
            frac_volc = frac_volc_arr[i3]
            
            print('tstart, fvolc, GEL',t_start,frac_volc,gwater_GEL)
            
        
            t,t_sun,Fxuv_venus,FH_ref,Fxuv_E = timestep(t_start,t_end,nsteps,epsilon,magma_yn)
            
            # Set up arrays to store model output during time evolution:
            moles_final = P_final/(g*Rmm_CO2) # Number of moles CO2 in Venus' modern atmosphere, could work for alternative atmospheres if P_final changes
            
            moles_CO2 = np.zeros(nsteps)
            moles_CO2[0] = ((1.-frac_volc)*P_final)/(g*Rmm_CO2) # Initialize atmosphere with CO2 from before end of habitable era, i.e. any CO2 not volcanically degassed in this model
            
            moles_H2O = np.zeros(nsteps)
            moles_H2O[0] = gwater_GEL*1000./Rmm_H2O # Initialize atmosphere with column H2O determined by habitable era water inventory - assumed all enters atmosphere
            moles_H2Oloss = np.zeros(nsteps) # cumulative moles H2O lost
            oceans = np.zeros(nsteps)
            
            moles_O2 = np.zeros(nsteps)
            moles_O2[0] = 0. # Atmosphere has no initial O2 because H2O disocciation and H escape have not yet occurred
            moles_atm = np.zeros(nsteps)
            moles_atm[0] = moles_CO2[0] + moles_H2O[0] + moles_O2[0] #total moles in atmosphere
            
            Rmm_atm = np.zeros(nsteps)
            Rmm_atm[0] = (moles_CO2[0]*Rmm_CO2 + moles_H2O[0]*Rmm_H2O + moles_O2[0]*Rmm_O2)/moles_atm[0]
            
            Pres = np.zeros(nsteps) # Pa
            Pres[0] = moles_atm[0] * Rmm_atm[0] * g # Calculate pressure from existing atmospheric constituents
            
            # Store total # moles O2 lost to different O sinks:
            moles_O2_magma = np.zeros(nsteps)
            moles_O2_nt = np.zeros(nsteps)
            moles_O2_ash = np.zeros(nsteps)
            moles_O2_lava = np.zeros(nsteps)
            O2_ox_melt = np.zeros(nsteps)
            
            # dz = np.zeros(nsteps) # Eruption rate (km Gyr^-1)
            
            T_surf = np.zeros(nsteps) # Surface temperature K - only explicitly calculated during runaway greenhouse
            z_melt = np.zeros(nsteps) # Depth of equivalent melt layer added at individual timestep, m
            T_melt = np.zeros(nsteps) # Basalt solidus at the surface as a function of pressure, K
            z_meltlayer_tot = np.zeros(nsteps) # Cumulative depth of melt throughout runaway greenhouse
            
            z_meltlayer = np.zeros(nsteps) # Actual equiv. depth of melt at each timestep
            
            meltfrac0 = np.zeros(nz) # For melt fraction as a function of depth at each timestep during surface melting
            ox_frac0 = np.zeros(nz) # Cumulative fraction oxidised material with depth at each timestep during melting
            latent = np.zeros(nz) # Latent heat release/consumption with depth at each timestep during melting
            melt_arr_new = np.zeros(nz)
            melt_arr_diff = np.zeros(nz)
            melt_arr_old = np.zeros(nz)
            melt_arr_cum = np.zeros(nz)
            melt_arr_cumdiff = np.zeros(nz)
            
            #%% ADDED TO TRACK H2O DISSOLUTION IN MELT
            diss_H2O = np.zeros(nz) # Mass fraction dissolved H2O in melt with depth
            
            # Set up arrays to track surface melting variables
            T_surf = np.zeros(nsteps)
            z_meltlayer = np.zeros(nsteps)
            
            Rmm_atm = np.zeros(nsteps)
            Rmm_atm[0] = (moles_CO2[0]*Rmm_CO2 + moles_H2O[0]*Rmm_H2O + moles_O2[0]*Rmm_O2)/moles_atm[0]
            Pres[0] = moles_atm[0] * Rmm_atm[0] * g
            
            # 1D heat diffusion in crust and melting
            dz = z[1]-z[0]
            
                
            T_arr = T_bg*np.ones(nz)  
            T_arr = T_arr + q*z/k
            T_arr0 = np.zeros(nsteps)
                
            meltfrac = np.zeros(nz)
            # meltfrac_out = np.zeros(nz)
            # ox_frac = np.zeros(nz)
            # ox_frac_new = np.zeros(nz)
            max_ox = np.zeros(nsteps) # max # moles oxygen removed per timestep
            ox_prop = np.zeros(nsteps) # fraction of oxygen in atmosphere removed by oxidation
            z_meltlayer_tot = np.zeros(nsteps)
            moles_O2_esc = np.zeros(nsteps)
            
            moles_diss = np.zeros(nsteps)
            delta_moles_diss_H2O = np.zeros(nsteps)
            
            # fig0,ax0 = plt.subplots()
            # fig0m,ax0m = plt.subplots()
            
            
            endmagma = 0
            
                
            
            for i1 in range(1,nsteps-1):
            
                oceans[i1] = moles_H2O[i1-1]/colEO
            
                if oceans[i1] <= 1.e-5:
                    T_surf[i1] = 700.
                    print('T_surf too low - set to 700')
                    break
                
                # Find surface temperature as a function of water in atmosphere
                T_surf[i1],z_melt = surface_T(oceans[i1],Pres[i1-1],t[i1],moles_H2O[i1-1],moles_atm[i1-1])
                # print('T_surf',T_surf[i1])
                
                # Set temperature at top of column of crust to surface temperature for diffusion calculation
                T_arr[0] = T_surf[i1]
                
                if T_surf[i1]<1273 and i1>=1000:
                # Terminate surface melting timestep sequence when surface below basalt solidus
                    if np.mean(T_surf[i1-1000:i1]) < 1273:
                        
                        test_oceans = (moles_H2O[i1-1] + moles_diss[i1-1])/colEO
                        
                        T_test,z_melt_test = surface_T(test_oceans,Pres[i1-1],t[i1],(moles_H2O[i1-1]+moles_diss[i1-1]),(moles_atm[i1-1]+moles_diss[i1-1]))
                        
                        if T_test < 1273:
                            z_meltlayer[i1] = 0
                    
                            # Update variables for return to main script
                            # endmagma = 1.
                            moles_H2O_out = moles_H2O[i1-1] 
                            moles_O2_out = moles_O2[i1-1]    
                            moles_H2Oloss = np.sum(moles_H2Oloss)
                            ash = np.sum(moles_O2_ash)
                            lava = np.sum(moles_O2_lava)
                            magma = np.sum(moles_O2_magma)
                            nonthermal = np.sum(moles_O2_nt)
                            
                             
                            
                            print('T_surf too low for melting, and its not coming back up')
                            break
            
                if z_meltlayer_tot[i1] < (0.99*z_end):
            
                    dTdz = (1/dz) * matrix.dot(T_arr)
                    dTdz = dTdz - (latent/k)
                    dTdz[nz-1] = dTdz[nz-1] + q/k 
                    T_new = T_arr + np.multiply((kappa * (t[i1]-t[i1-1])/ (dz)), dTdz)
                    T_arr=T_new
                    T_arr0[i1] = T_arr[0]
                    # if i1 % (nsteps/12) == 0:
                        # pltt_arr = ax0.plot(z,T_arr)
                        
                    T_arr[0]=T_surf[i1]
                
                # Terminate run if surface temperature is too high - usually a numerical stability error generated by timestep size
                if np.max(T_arr) > 3500.:
                    print('ERROR: UNSTABLE','max T',np.max(T_arr))
                    break
            
                # Calculate surface melt fraction to determine if surface melting still occuring
                # print('max T in array',np.max(T_arr))
                
                # mf1 = melt_function(T_arr[0],(P_z[0]+Pres[i1]))
                test_melt = float(knn.predict(np.array([[(P_z[0]+Pres[i1])/1.e5,T_arr[0]]])))
                # test_melt = mf1.item() 
                
                
                # If no melting, terminate surface melting sequence
                # if test_melt <= 0.01:
                    
                #     z_meltlayer[i1] = 0
                #     print('end of surface melting, time',(t1/secingyr)*1e9,'T_surf',T_surf[i1])
                    
                #     endmagma = 1.
                #     moles_H2O_out = moles_H2O[i1-1]
                #     moles_O2_out = moles_O2[i1-1]    
                #     moles_H2Oloss = np.sum(moles_H2Oloss)
                #     ash = np.sum(moles_O2_ash)
                #     lava = np.sum(moles_O2_lava)
                #     magma = np.sum(moles_O2_magma)
                #     nonthermal = np.sum(moles_O2_nt)
                    
                #     break
                
                # If melting, continue surface melting sequence
                if endmagma == 0:
                    
                    endmagma = 0.
                    
                    a = np.array([(P_z+Pres[i1])/1.e5,T_arr])
                    b = np.fliplr(np.rot90(a,axes=(1,0)))
                    melt_arr_new = (knn.predict(b))
                    
                    # if i1 % nsteps/100000 == 0:
                    #     plt.plot(b[:,0],melt_arr_new)
                    
                    # break
                    
                    for i2 in range(np.size(P_z>10000*1e5),np.size(nz)):
                        if T_arr[i2]>1300:
                                
                            if (P_z[i2]+Pres[i1])>=10000*1e5:
                                mf10 = mf_10k(T_arr[i2])
                                melt_10 = mf10.item()
                                mf15 = mf_10k(T_arr[i2])
                                melt_15 = mf15.item()
                                
                                split = np.abs(arr_10_15 - (P_z[i2]+Pres[i1]))
                                frac = 5000/split
                                
                                melt_arr_new[i2]=((mf10*split[0]) + (mf15*split[1]))/np.sum(frac)
                                
                            if melt_arr_new[i2]<0:
                                melt_arr_new[i2]=0
                            elif melt_arr_new[i2]>1:
                                melt_arr_new[i2]=1
                    # melt_arr_new[T_arr>1200] = melt_function(T_arr[T_arr>1200],P_z[T_arr>1200])
                    
                    # melt_arr_diff = melt_arr_new - melt_arr_old
                    # print('added melt',np.mean(melt_arr_diff))
                    
                    for i2 in range(0,nz-1):
                        if melt_arr_cum[i2]<melt_arr_new[i2]:
                            melt_arr_cumdiff[i2] = melt_arr_new[i2] - melt_arr_cum[i2]
                    
                    z_meltlayer[i1] = np.mean(melt_arr_new)*z_end
                    
                    # oxidizable_melt = (melt_arr_new - ox_frac)
                    # oxidizable_melt[oxidizable_melt<0] = 0.
                    
                    max_ox[i1] = (((0.5*rho*z_end*np.mean(melt_arr_new)) *max_water)/Rmm_H2O) - moles_O2_magma[i1-1]
                    
                    if max_ox[i1] < 0.:
                        max_ox[i1] = 0.
                    # print('max_ox',max_ox[i1])
                    
                    if melt_arr_new[0] < 0.01:
                        max_ox[i1] = 0
                    
                    #%% water dissolution
                    partial_P_H2O = (moles_H2O[i1-1]/moles_atm[i1-1])*Pres[i1-1]
                    mass_frac_H2O = 3.44e-8 * ((partial_P_H2O)**0.74)
                    
                    delta_moles_diss_H2O[i1] = mass_frac_H2O*(z_meltlayer[i1] - z_meltlayer[i1-1])*rho/Rmm_H2O
                    # print('dissolved H2O moles this timestep',delta_moles_diss_H2O)
                    
                    # if z_meltlayer[i1] < z_meltlayer[i1-1]:
                        
                        # print('shrinking melt layer',z_meltlayer[i1] - z_meltlayer[i1-1])
                        
                    if z_meltlayer[i1] < 0:
                        print('negative melt layer',z_meltlayer[i1])
                        break
                    
                    if moles_H2O[i1-1] > delta_moles_diss_H2O[i1]:
                        # print('more H2O in atm')
                        moles_H2O[i1-1] = moles_H2O[i1-1] - delta_moles_diss_H2O[i1]
                    elif moles_H2O[i1-1] < delta_moles_diss_H2O[i1]:
                        print('less H2O in atm')
                        mass_frac_H2O = mass_frac_H2O * (moles_H2O[i1-1]/delta_moles_diss_H2O[i1])
                        delta_moles_diss_H2O[i1] = mass_frac_H2O*np.mean(melt_arr_diff)*z_end*rho*mass_frac_H2O/Rmm_H2O
                        moles_H2O[i1-1] = moles_H2O[i1-1] - delta_moles_diss_H2O[i1]
                        
                    moles_diss[i1] = moles_diss[i1-1] + delta_moles_diss_H2O[i1]
                    
                    moles_diss_max = mass_frac_H2O *rho* z_meltlayer[i1]/Rmm_H2O
                    
                    if moles_diss[i1] > moles_diss_max:
                        moles_H2O[i1-1] = moles_H2O[i1-1] + (moles_diss[i1] - moles_diss_max)
                        moles_diss[i1] = moles_diss[i1] - (moles_diss[i1] - moles_diss_max)
                        
                        # break
                #%% O2 loss
                    
                
                
                #%% H2O dissolution
                
                #%%
                #     # if i1 % (nsteps/2) == 0:
                #     # Determine depths at which melting is occuring
                #     for i2 in range(1,nz-1):
                #         if T_arr[i2] < 1200.:
                #             test_melt = 0.
                #             break
                #         elif T_arr[i2] >= 1200.:
                #             test_melt = melt_function(T_arr[i2])
            
                #         if test_melt <= 0.:
                #             test_melt = 0.
                #             break
                #         if test_melt >= 1.:
                #             test_melt = 1.
                            
                #         if i2 == 0 and test_melt < 0.01:
                #             print('no melting at surface T')
                #             break
                #         melt_arr_new[T_arr>1200] = melt_function(T_arr[T_arr>1200])
                #         ax0m.plot(z,melt_arr)
                        
                        
                #         # Find cumulative melt fraction 
                #         cummelt = meltfrac_out[i2]
                        
                #         # Find melt fraction added at timestep i1
                #         if cummelt > 1:
                #             exit('Error: cumulative melt fraction >1')
                        
                #         if cummelt == 1.:
                #             meltfrac[i2] = 0.
                            
                #         elif test_melt - cummelt <= 0.:
                #             meltfrac[i2] = 0.
                            
                #         elif test_melt - cummelt > 0.:
                #             meltfrac[i2] = test_melt - cummelt
                            
                #         meltfrac_out[i2] = cummelt + meltfrac[i2]
            
                #     else:   
                #         meltfrac = meltfrac*0.
            
                #     melt_ox = meltfrac - ox_frac
                    
                #     for i2 in range(0,nz-1):
                        
                #         if melt_ox[i2]<= 0.:
                #             melt_ox[i2] = 0.
                    
                #     # Track cumulative thickness of surface melt layer at each timestep
                #     z_meltlayer[i1] = np.mean(melt_ox)*(z_end-dz)
                #     z_meltlayer_tot[i1] =  z_meltlayer_tot[i1] + z_meltlayer[i1]
                #%%
                    
                    # Latent heat consumed/released at timestep
                    # if i1 % (nsteps/2) == 0:
                    latent = melt_arr_diff * (z_end/nz) * lat_melt * rho/(t[i1]-t[i1-1])
            
                    if z_meltlayer_tot[i1]>(z_end-dz):
                        print('melt layer too thick!',np.sum(z_meltlayer[0:i1]),'melt layer',z_meltlayer[0:i1])
                        
                    # Calculate O loss to melt oxidation
                    # O2_ox_melt[i1] = (0.5*rho*z_meltlayer[i1]*max_water)/Rmm_H2O
                    
                # Set up variables for RK4 integration
                if i1 == 1:
                    rk1_CO2 = 0.
                    rk1_H2O = 0.
                    rk1_O2 = 0.
                    rk1_moles_O2_nt = 0.
                elif i1 > 1: 
                    rk1_CO2 = moles_CO2[i1-1] - moles_CO2[i1-2]
                    rk1_H2O = moles_H2O[i1-1] - moles_H2O[i1-2]
                    rk1_O2 = moles_O2[i1-1] - moles_O2[i1-2]
                    rk1_moles_O2_nt = moles_O2_nt[i1-1] - moles_O2_nt[i1-2]
                
                thalf = t[i1-1] + 0.5*(t[i1] - t[i1-1])
            
                rk2_CO2,rk2_H2O,rk2_O2,rk2_Pres,rk2_moles_O2_ash,rk2_moles_O2_lava,rk2_moles_O2_nt = stdrk4(degas,0,0,FH_ref[i1],T,T_atm,0,0,0,(moles_CO2[i1-1] + (0.5*rk1_CO2)),(moles_H2O[i1-1] + (0.5*rk1_H2O)),(moles_O2[i1-1] + (0.5*rk1_O2)),t[i1-1],thalf,t1,0,0,0)
                rk3_CO2,rk3_H2O,rk3_O2,rk3_Pres,rk3_moles_O2_ash,rk3_moles_O2_lava,rk3_moles_O2_nt = stdrk4(degas,0,0,FH_ref[i1],T,T_atm,0,0,0,(moles_CO2[i1-1] + (0.5*rk2_CO2)),(moles_H2O[i1-1] + (0.5*rk2_H2O)),(moles_O2[i1-1] + (0.5*rk2_O2)),t[i1-1],thalf,t1,0,0,0)
                rk4_CO2,rk4_H2O,rk4_O2,rk4_Pres,rk4_moles_O2_ash,rk4_moles_O2_lava,rk4_moles_O2_nt = stdrk4(degas,0,0,FH_ref[i1],T,T_atm,0,0,0,(moles_CO2[i1-1] + rk3_CO2),(moles_H2O[i1-1] + rk3_H2O),(moles_O2[i1-1] + rk3_O2),t[i1-1],t[i1],t1,0,0,0)
                
                moles_CO2[i1] = moles_CO2[i1-1] + ((1/6) * (rk1_CO2 + (2*rk2_CO2) + (2*rk3_CO2) + rk4_CO2))
                moles_H2O[i1] = moles_H2O[i1-1] + ((1/6) * (rk1_H2O + (2*rk2_H2O) + (2*rk3_H2O) + rk4_H2O))
                moles_O2[i1] = moles_O2[i1-1] + ((1/6) * (rk1_O2 + (2*rk2_O2) + (2*rk3_O2) + rk4_O2))
                
                
                moles_O2_esc[i1] = (moles_O2[i1] - moles_O2[i1-1])
                # print('moles_O2 from atm escape of H',(moles_O2_esc[i1]))
            
                moles_O2_ash[i1] = 0.
                moles_O2_lava[i1] = 0.
                moles_O2_nt[i1] = moles_O2_nt[i1-1] + ((1/6) * (rk1_moles_O2_nt + (2*rk2_moles_O2_nt) + (2*rk3_moles_O2_nt) + rk4_moles_O2_nt))
                
                if endmagma == 1:
                    print('endmagma==1')
                    break
                
                # # If melt layer deeper than 50km, assume enough oxidizeable material present and in contact with atmosphere to remove all O2 at each timestep - assumes mixing in melt
                # # A 50km melt layer is able to remove 3000 m GEL water which exceeds initial H2O used in this model
                # if z_meltlayer_tot[i1] >= (z_end-dz):
                #     O2_ox_melt[i1] = moles_O2[i1]
                #     moles_O2[i1] = 0.
                #     ox_frac = (ox_frac * 0.) + 1.
                #     print('all oxidising!')
                    
                # # Find O2 lost to melt layer oxidation at each timestep
                # elif z_meltlayer_tot[i1] < (z_end-dz):
                #     if O2_ox_melt[i1] >= moles_O2[i1]:
                #         moles_O2_magma[i1] = moles_O2[i1]
                #         moles_O2[i1] = 0.
                        
                #         for i2 in range(0,nz-1):
                #             ox_frac_new = (moles_O2[i1]/O2_ox_melt[i1]) * meltfrac[i2]
            
                #             if (ox_frac_new + ox_frac[i2])>= meltfrac_out[i2]:
                #                 ox_frac_new = (meltfrac_out[i2] - ox_frac[i2])
                #                 ox_frac[i2] = meltfrac_out[i2]
                #             elif (ox_frac_new + ox_frac[i2]) < meltfrac_out[i2]:
                #                 ox_frac[i2] = ox_frac0[i2] + ox_frac_new
            
                #     else:
                #         moles_O2_magma[i1] = O2_ox_melt[i1]
                #         moles_O2[i1] = moles_O2[i1] - O2_ox_melt[i1]
                        
                #         for i2 in range(0,nz-1):
                #             ox_frac_new = meltfrac[i2]
                            
                #             if (ox_frac_new + ox_frac[i2])>= meltfrac_out[i2]:
                #                 ox_frac_new = (meltfrac_out[i2] - ox_frac[i2])
                #                 ox_frac[i2] = meltfrac_out[i2]
                #             elif (ox_frac_new + ox_frac[i2]) < meltfrac_out[i2]:
                #                 ox_frac[i2] = ox_frac0[i2] + ox_frac_new
                
                if max_ox[i1]>0:
                    ox_prop[i1] = (moles_O2[i1]/max_ox[i1]) 
                    # print('ox prop',ox_prop[i1])
                if max_ox[i1]<0:
                    max_ox[i1]=0
                    ox_prop[i1]=0
                    
                
                
                # ox_frac_new = ox_frac + (ox_prop[i1]*melt_arr_new)
                # print('new ox_frac',ox_frac_new[0:10],'old',ox_frac[0:10])
                
                # for i2 in range(0,nz-1):
                #     if ox_frac_new[i2]>1:
                #         ox_frac_new[i2] = 1.
            
                
                if max_ox[i1] > moles_O2[i1]:
                    moles_O2_magma[i1] = moles_O2_magma[i1-1] + moles_O2[i1]
                    moles_O2[i1] = 0.
                elif max_ox[i1] < moles_O2[i1]:
                    moles_O2_magma[i1] = moles_O2_magma[i1-1] + (ox_prop[i1]*max_ox[i1])
                    moles_O2[i1] = moles_O2[i1] - (ox_prop[i1]*max_ox[i1])
                
                # #Update melt array w/ depth
                # for i2 in range(0,nz-1):
                #     if melt_arr_new[i2]>melt_arr_cum[i2]:
                #         # print('new melting, update cumulative array',i2)
                #         melt_arr_cum[i2] = melt_arr_new[i2]
                
                # if i1 == 1000:                                                                                                                       
                #     break
                
                melt_arr_old = melt_arr_new
                # ox_frac = ox_frac_new
                
                # Update atmosphere
                moles_atm[i1] = moles_CO2[i1] + moles_H2O[i1] + moles_O2[i1]
                Rmm_atm[i1] = (moles_CO2[i1]*Rmm_CO2 + moles_H2O[i1]*Rmm_H2O + moles_O2[i1]*Rmm_O2)/moles_atm[i1]
                
                # Update total pressure
                Pres[i1] = moles_atm[i1] * Rmm_atm[i1] * g
                
                
            #%%
            results_add = np.array([t_start,t[i1-1]/secingyr,Pres[i1-1],moles_atm[i1-1],Rmm_atm[i1-1],
                                                moles_CO2[i1-1],moles_H2O[i1-1],moles_diss[i1-1],T_surf[i1-1],
                                                moles_O2[i1-1],np.sum(moles_O2_nt),np.sum(moles_O2_magma),
                                                np.max(z_meltlayer)])
                            
            results = np.vstack([results,results_add])
            
            gwater_string = str(int(gwater_GEL))
            f_volc_string = str(int(10*frac_volc))
            CO2_string = str(int(1e6*mCO2tot))
            t_string = str(int(10*t_start))
            
            title = ("withdiss_t" + t_string + "_gwater" + gwater_string + "_f_volc" + f_volc_string + ".csv")
            
            data = {'tstart':t_start,
                            'f_volc':frac_volc,
                            'gwater':gwater_GEL,
                            'time':t[0:i1],
                            'Pressure':Pres[0:i1],
                            'moles_atm':moles_atm[0:i1],
                            'moles_CO2':moles_CO2[0:i1],
                            'moles_H2O':moles_H2O[0:i1],
                            'moles_O2':moles_O2[0:i1],
                            'nonthermal':moles_O2_nt[0:i1],
                            'oxidation':moles_O2_magma[0:i1],
                            'Rmm_atm':Rmm_atm[0:i1],
                            'T_surf':T_surf[0:i1],
                            'dissolved':moles_diss[0:i1],
                            'z_melt':z_meltlayer[0:i1]}
                        
            df = pd.DataFrame(data)
            
            df.to_csv(title,header = True)
            # #%% what's going on
            # fig1, ax1= plt.subplots ()
            # plot_Tsurf = plt.plot(t[0:i1]/secingyr,T_surf[0:i1])
            # plt.ylim((1000,1800))
            
            # #%%
            
            # fig2, ax2 = plt.subplots(figsize=(20,15))
            # plot_colmass, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_atm[0:i1]*Rmm_atm[0:i1],label='Total',linewidth=4,color='black')
            # plot_colCO2, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_CO2[0:i1]*Rmm_CO2,label='CO\N{SUBSCRIPT TWO}',linewidth=4,color ='indianred')
            # plot_colH2O, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_H2O[0:i1]*Rmm_H2O,label='H\N{SUBSCRIPT TWO}O',linewidth=4,color ='slateblue')
            # plot_colO2, = plt.plot(t[0:i1]/(1.e9*365.*60.*60.*24.),moles_O2[0:i1]*Rmm_O2,label='O\N{SUBSCRIPT TWO}',linewidth=4,color ='lightseagreen')
            # plt.yscale("log")
            # plt.axhline(y=(100.e-6 * moles_atm[i1] * Rmm_H2O),color ='slateblue', linestyle='--',linewidth=4)
            # plt.axhline(y=(50.e-6 * moles_atm[i1] * Rmm_O2), color ='lightseagreen', linestyle='--',linewidth=4)
            # # plt.axhline(y=1, color ='black', linestyle='--',linewidth=2)
            # plt.ylim((1.e-6,1.e8))
            # # plt.xlim(2,2.1)
            # plt.legend(handles=[plot_colmass,plot_colCO2,plot_colH2O,plot_colO2], fontsize=40)
            # plt.title(str(float('%.1g' % (1.e6*mCO2tot))) + ' ppm CO$\mathregular{_{2}}$, ' + str(float('%.1g' % (1.e2*mH2Otot)))+ ' wt% H$\mathregular{_{2}}$O' , fontsize=40) 
            # ax2.set_ylabel('Column Mass (kg)', fontsize=40)
            # ax2.set_xlabel('Time (Gyr)', fontsize=40)
            # plt.gca().set_ylim(bottom=1e0)
            # ax2.tick_params(axis='both', which='major', labelsize=30)
            
            # #%%
            # fig3, ax3= plt.subplots()
            # plot_meltlayer = plt.plot(t[0:i1]/secingyr,z_meltlayer[0:i1])
            
            # #%%
            # fig4, ax4= plt.subplots()
            # plot_meltlayer = plt.plot(t[0:i1]/secingyr,moles_O2_magma[0:i1])
            
            # #%%
            
            # plt.plot(t[0:-1]/secingyr,(z_meltlayer[1:] - z_meltlayer[0:-1]))
            # plt.yscale("log")
            
            
            # #%%
            # plt.plot(t/secingyr,T_surf)
            # plt.axhline(y=(1283), linestyle='--',linewidth=1)
            # plt.ylim((1100,1600))