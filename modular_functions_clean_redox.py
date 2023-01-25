# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 09:16:21 2022

@author: sasha
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 12:50:57 2021

@author: aowarren
"""

from scipy.io import loadmat
import numpy as np
import time
import math
from VolcGases.functions import solve_gases
from scipy import optimize,interpolate
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import Rbf
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
import sys
import pathlib
from pathlib import Path
import random

threshold_time = (1*24*60*60) #+ (11.5*60*60) #time after which code should be terminated and saved!

secingyr = 1.e9 * (365.*60.*60.*24.)

# Atmospheric species
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
A = 25738
B = 9
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

z_end = 50.e3 # max depth of interest for surface melt layer
z = np.linspace(0,z_end,nz)

# Import basalt solidus from alphaMELTS calculations
solidus = np.genfromtxt('basalt_sol.csv',delimiter=',',skip_header=1)

# Import melt fraction data from alphaMELTS results
melting = np.genfromtxt('melt_f.csv',delimiter=',',skip_header=1)
melt_T = melting[:,1] #K
melt_f = melting[:,3] #fraction
# Find fit to melt fraction as a function of temperature from alphaMELTS results
melt_function = interpolate.interp1d(melt_T,melt_f,kind='linear')

#%%
luminosity = np.genfromtxt('luminosity.csv',delimiter=',')
solart = luminosity[:,0] #Gyr
solarL = luminosity[:,1]/(AU**2.) #relative to Earth present day 
L = interpolate.interp1d(solart,solarL,kind='linear')

n_run = 100
Pgrid = np.linspace(0,200,num=n_run)

mult_arr = np.linspace(0.7,1.3,3)

time_meltplot = np.linspace(1.5,1.51,num=11)

#%%

def porespace(tstart,q,phi0):
    # This function finds the pore closure depth and global equivalent layer depth of groundwater
    # For use with script "gwater.py"
    
    k = 3.
    A = 6.12
    Q = 276.e3
    n = 3.05
    rho_crust = 2800. #kg m^-3
    R = 8.31
    g = 8.87
    Ts = 400. #K
    
    t = tstart*secingyr #convert to seconds
    
    dz = 1.
    z = np.arange(0.,5.e4,dz)
    P = z*g*rho_crust
    
    PMPa = P/1.e6 # convert pressure to megapascals
    
    porosity = np.zeros(len(z))
    T = np.zeros(len(z)) + Ts
    
    T = T + z*q/k
    eta = ((PMPa**(1.-n))/A)*np.exp(Q/(R*T))
    porosity = phi0*np.exp((-P*t)/eta)

    z_closure = z[np.amin(np.where(porosity<phi0/math.exp(1)))]
    
    phi_mean = np.mean(porosity[0:np.amin(np.where(porosity<phi0/math.exp(1)))])
    
    gwater_GEL = z_closure*phi_mean
    
    return(gwater_GEL)
#%%

def time_evol(tstart,nsteps,metamorphic_degas,frac_volc,P_final,b_err,dzdt_pres,gwater_GEL,T,epsilon,T_atm,mCO2tot,mH2Otot,FMQ,frac_ext,q): 
    # This function solves for the time evolution of Venus' atmosphere:
        
    start_time = time.time()    
    
    # Calculate look-up table of H2O and CO2 solubilities using VolcGasses (Wogan et al. 2020)
    # Used to find H2O degassed at each timestep, and total volume of degassed H2O and CO2 at each timestep to determine whether volcanism is explosive
    sH2O = np.zeros(n_run) #solubility of H2O, kg/mol
    sCO2 = np.zeros(n_run) #solubility of CO2, kg/mol
    CO_CO2_rat = np.zeros(n_run)
    CH4_CO2_rat = np.zeros(n_run)
    H2_H2O_rat = np.zeros(n_run)
    
    for i1p in range(1,n_run):
        Pn = Pgrid[i1p]        
        log_FMQ = (-A/T+B+C*(Pn-1)/T)
        f_O2 = 10**(log_FMQ+FMQ)
        x = 0.01550152865954013 #defined as mol/g of magma
        
        P_H2O,P_H2,P_CO2,P_CO,P_CH4,alphaG,x_CO2,x_H2O = solve_gases(T,Pn,f_O2,mCO2tot,mH2Otot)
        
        #CO2 and H2O solubility as mass fraction
        H2O_degas = (1000*alphaG*x*(1/(1-alphaG))*P_H2O/Pn)
        H2_degas = (1000*alphaG*x*(1/(1-alphaG))*P_H2/Pn)
        CO2_degas = (1000*alphaG*x*(1/(1-alphaG))*P_CO2/Pn)
        CO_degas = (1000*alphaG*x*(1/(1-alphaG))*P_CO/Pn)
        CH4_degas = (1000*alphaG*x*(1/(1-alphaG))*P_CH4/Pn)
        
        
        sH2O[i1p] = x*x_H2O*Rmm_H2O*1000 #factor of 1000 converts from kg/mol to g/mol for RMM
        sCO2[i1p] = x*x_CO2*Rmm_CO2*1000
        CO_CO2_rat[i1p] = CO_degas/(CH4_degas+CO_degas+CO2_degas)
        CH4_CO2_rat[i1p] = CH4_degas/(CH4_degas+CO_degas+CO2_degas)
        H2_H2O_rat[i1p] = H2_degas/(H2_degas+H2O_degas)

    if frac_volc == 0.:
        metamorphic_degas = 1
        
#%%
    # Generate file name components to search for existing runs, if needed
    gwater_string = str(int(gwater_GEL))
    h2o_string = str(int(1000*mH2Otot))
    f_ext_string = str(int(10*frac_ext))
    f_volc_string = str(int(10*frac_volc))
    CO2_string = str(int(1e6*mCO2tot))
    t_string = str(int(10*tstart))
    FMQ_string = str(int(FMQ))
    
    # Check for existing calculations with input parameters (useful for re-starting time-consuming runs from part way through)
    title_fail = ("results/no_melt_change/ZFAIL_NMtstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + "_FMQ-" + FMQ_string + ".csv")
    title_complete = ("results/no_melt_change/NMtstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + "_FMQ-" + FMQ_string + ".csv")
            
    file = pathlib.Path(title_complete)
    filef = pathlib.Path(title_fail)
    if filef.exists ():
        print('already done this run!')
        title = title_fail
    if file.exists ():
        print('already done this run!')
        title = title_complete
    else:
        title_partial = ("results/partial/partial_NMtstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + "_FMQ-" + FMQ_string + ".csv")
        
        file = pathlib.Path(title_partial)
        
        if file.exists ():
            print('existing run found, continuing from here!')
            title = title_partial
        else:
            title_magma = ("magmaocean_only/magmaonly_tstart" + t_string + 'water' + gwater_string + "_f_volc" + f_volc_string + "_FMQ-" + FMQ_string + '.csv')
            file = pathlib.Path(title_magma)
            if file.exists ():
                print('magma ocean calculations found, continuing from here!')
                title = title_magma
            else:
                title = str("none")
    

#%%
    # Sequence for restarting from existing partial runs:
    title = str('none')
    
    if "none" in title:
        print('no existing runs, start from scratch!')
        begin = 1 # begin may be >1 for runs recovered from partially completed run files
        
        t_init = tstart*secingyr
        i1_init = 0
        # Set up arrays to store model output during time evolution:
        moles_final = P_final/(g*Rmm_CO2) # Number of moles CO2 in Venus' modern atmosphere, could work for alternative atmospheres if P_final changes
        
        moles_CO2 = np.zeros(nsteps)
        moles_CO2[0] = ((1.-frac_volc)*P_final)/(g*Rmm_CO2) # Initialize atmosphere with CO2 from before end of habitable era, i.e. any CO2 not volcanically degassed in this model
        
        moles_CO = np.zeros(nsteps)
        moles_H2 = np.zeros(nsteps)
        moles_CH4 = np.zeros(nsteps)
        
        moles_H2O = np.zeros(nsteps)
        moles_H2O[0] = gwater_GEL*1000./Rmm_H2O # Initialize atmosphere with column H2O determined by habitable era water inventory - assumed all enters atmosphere
        oceans = np.zeros(nsteps)
        
        moles_O2 = np.zeros(nsteps)
        moles_O2[0] = 0. # Atmosphere has no initial O2 because H2O disocciation and H escape have not yet occurred
        moles_atm = np.zeros(nsteps)
        moles_atm[0] = moles_CO2[0] + moles_H2O[0] + moles_O2[0] + moles_CO[0] + moles_H2[0] #total moles in atmosphere
        
        Rmm_atm = np.zeros(nsteps)
        Rmm_atm[0] = (moles_CO2[0]*Rmm_CO2 + moles_H2O[0]*Rmm_H2O + moles_O2[0]*Rmm_O2 + moles_CO[0]*Rmm_CO + moles_H2[0]*Rmm_H2 + moles_CH4[0]*Rmm_CH4)/moles_atm[0]
        
        Pres = np.zeros(nsteps) # Pa
        Pres[0] = moles_atm[0] * Rmm_atm[0] * g # Calculate pressure from existing atmospheric constituents
        
        # Store total # moles O2 lost to different O sinks:
        moles_O2_magma = np.zeros(nsteps)
        moles_O2_nt = np.zeros(nsteps)
        moles_O2_ash = np.zeros(nsteps)
        moles_O2_lava = np.zeros(nsteps)
        moles_O2_red = np.zeros(nsteps)
        
        dz = np.zeros(nsteps) # Eruption rate (km Gyr^-1)
        
        T_surf = np.zeros(nsteps) # Surface temperature K - only explicitly calculated during runaway greenhouse
        z_melt = np.zeros(nsteps) # Depth of equivalent melt layer added at individual timestep, m
        T_melt = np.zeros(nsteps) # Basalt solidus at the surface as a function of pressure, K
        z_meltlayer_tot = np.zeros(nsteps) # Cumulative depth of melt throughout runaway greenhouse
        
        meltfrac0 = np.zeros(nz) # For melt fraction as a function of depth at each timestep during surface melting
        ox_frac0 = np.zeros(nz) # Cumulative fraction oxidised material with depth at each timestep during melting
        latent = np.zeros(nz) # Latent heat release/consumption with depth at each timestep during melting
        
    else:
        startup = np.genfromtxt(title,delimiter=',',skip_header=0)
        
        # For loading existing run data from file so it can be used to resume run 
        
        if np.ndim(startup)>1:
            i1_init = np.nan_to_num(startup[:,0])
            si1 = int(np.max(i1_init))+1
            
            t_init = startup[0:si1,2]
            z_erupt_init = startup[0:si1,6]
            Pres_init = startup[0:si1,12]
            moles_atm_init = startup[0:si1,13]
            moles_CO2_init = startup[0:si1,14]
            moles_CO_init = startup[0:si1,15]
            moles_H2_init = startup[0:si1,16]
            moles_H2O_init = startup[0:si1,17]
            moles_O2_init = startup[0:si1,18]
            Rmm_atm_init = startup[0:si1,19]
            T_atm_init = startup[0:si1,20]
            ash_init = startup[0:si1,21]
            lava_init = startup[0:si1,22]
            magma_init = startup[0:si1,23]
            nonthermal_init = startup[0:si1,24]
            reducing_init = startup[0:si1,25]
            z_meltlayer_tot_init = startup[0:si1,26]
            moles_CH4_init = startup[0:si1,27]
    
        elif np.ndim(startup)==1:
            i1_init = startup[0]
            t_init = startup[2]
            z_erupt_init = startup[6]
            Pres_init = startup[12]
            moles_atm_init = startup[13]
            moles_CO2_init = startup[14]
            moles_CO_init = startup[15]
            moles_H2_init = startup[16]
            moles_H2O_init = startup[17]
            moles_O2_init = startup[18]
            Rmm_atm_init = startup[19]
            T_atm_init = startup[20]
            ash_init = startup[21]
            lava_init = startup[22]
            magma_init = startup[23]
            nonthermal_init = startup[24]
            reducing_init = startup[25]
            z_meltlayer_tot_init = startup[26]
            moles_CH4_init = startup[27]
            
            si1 = 1
            
    
        tstart = np.max(t_init)/secingyr

        moles_final = P_final/(g*Rmm_CO2)
        Pres = np.zeros(nsteps)
        moles_CO2 = np.zeros(nsteps)
        moles_CO2[0:si1] = moles_CO2_init
        
        moles_CO = np.zeros(nsteps)
        moles_CO[0:si1] = moles_CO_init
        moles_CH4 = np.zeros(nsteps)
        moles_CH4[0:si1] = moles_CH4_init
        moles_H2 = np.zeros(nsteps)
        moles_H2[0:si1] = moles_H2_init
        
        moles_H2O = np.zeros(nsteps)
        moles_H2O[0:si1] = moles_H2O_init
        oceans = np.zeros(nsteps)
        moles_O2 = np.zeros(nsteps)
        moles_O2[0:si1] = moles_O2_init
        moles_atm = np.zeros(nsteps)
        moles_atm[0:si1] = moles_atm_init
        Rmm_atm = np.zeros(nsteps)
        Rmm_atm[0:si1] = Rmm_atm_init
        Pres[0] = moles_atm[0] * Rmm_atm[0] * g
        Pres[0:si1] = Pres_init
        
        moles_O2_magma = np.zeros(nsteps)
        moles_O2_nt = np.zeros(nsteps)
        moles_O2_ash = np.zeros(nsteps)
        moles_O2_lava = np.zeros(nsteps)
        moles_O2_red = np.zeros(nsteps)
        moles_O2_magma[0:si1] = magma_init
        moles_O2_nt[0:si1] = nonthermal_init
        moles_O2_lava[0:si1] = lava_init
        moles_O2_ash[0:si1] = ash_init
        moles_O2_red[0:si1] = reducing_init
        
        dz = np.zeros(nsteps)
        dz[0:si1] = z_erupt_init
        
        T_surf = np.zeros(nsteps)
        T_surf[0:si1] = T_atm_init
        z_melt = np.zeros(nsteps)
        T_melt = np.zeros(nsteps)
        z_meltlayer_tot = np.zeros(nsteps)
        z_meltlayer_tot[0:si1] = z_meltlayer_tot_init
        meltfrac0 = np.zeros(nz)
        ox_frac0 = np.zeros(nz)
        latent = np.zeros(nz)
        
        begin = si1
    
    if begin >= 0.99*nsteps:
        # If model loaded from previous file is already complete, return output to save as completed run
        print('run already completed')
        return(t_init,i1_init,Pres,moles_atm,Rmm_atm,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,moles_O2,dz,T_surf,moles_O2_ash,moles_O2_lava,moles_O2_magma,moles_O2_nt,moles_O2_red,z_meltlayer_tot)

    #%%

    magma_yn = 0. # All models start with magma_yn = 0 i.e. no surface melting. When magma_yn = 1 (calculated later), model will run at smaller timesteps
    # Find time steps as a function of model start and whether or not surface melting is occuring
    t,t_sun,Fxuv_venus,FH_ref = timestep(tstart,age_sun,nsteps,epsilon,magma_yn)

    if metamorphic_degas == 0:
        nn,m = eruption_rate(nsteps,moles_CO2[0],frac_ext,P_final,mCO2tot,age_sun,tstart,b_err,dzdt_pres,frac_volc)
        dz_check = m*np.exp(-b_err*t_sun) + nn
        total_erupt = np.sum(dz_check[1:nsteps]*((t[1:nsteps]-t[0:nsteps-1])/secingyr))
        
    elif metamorphic_degas == 1:
        # If all modern CO2 in Venus' atmosphere assumed to be degassed before/at end of habitable era:
        dz = dz*0.

   #%%     
    endmagma = 0.
    
    for i1 in range(begin,nsteps):
        
        # calculate upper limit on total erupted thickness based on constraints from internal heating rate and models
        z_max = (4.5-tstart)*1.e9*((4/3) * 175.e9)/(4.6e14)

        # discontinue run if not possible to degas remaining needed CO2 without violating upper limit on eruption rate
        if total_erupt > z_max:
            
            print('z_erupt too large, abort run',total_erupt,'z_max',z_max)
            
            moles_atm[i1:nsteps] = moles_atm[i1]
            moles_H2O[i1:nsteps] = moles_H2O[i1]
            moles_CO2[i1:nsteps] = moles_CO2[i1]
            moles_CO[i1:nsteps] = moles_CO[i1]
            moles_CH4[i1:nsteps] = moles_CH4[i1]
            moles_H2[i1:nsteps] = moles_H2[i1]
            moles_O2[i1:nsteps] = moles_O2[i1]
            
            moles_O2_ash[i1:nsteps] = moles_O2_ash[i1]
            moles_O2_lava[i1:nsteps] = moles_O2_lava[i1]
            moles_O2_magma[i1:nsteps] = moles_O2_magma[i1]
            moles_O2_nt[i1:nsteps] = moles_O2_nt[i1]
            moles_O2_red[i1:nsteps] = moles_O2_red[i1]
            
            return(t,i1,Pres,moles_atm,Rmm_atm,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,moles_O2,dz,T_surf,moles_O2_ash,moles_O2_lava,moles_O2_magma,moles_O2_nt,moles_O2_red,z_meltlayer_tot,z_max)
            
            break
        
        # Calculate crustal production at current timestep
        dz[i1] = m*np.exp(-b_err*(t[i1-1]/secingyr)) + nn
        
        # end run if eruption rate calculation has failed and generated -ve values (this shouldn't happen unless something is wrong with your input times!)
        if dz[i1] < 0:
            
            print('ERROR: z_erupt negative, abort run',total_erupt)
            
            moles_atm[i1:nsteps] = moles_atm[i1]
            moles_H2O[i1:nsteps] = moles_H2O[i1]
            moles_CO2[i1:nsteps] = moles_CO2[i1]
            moles_CO[i1:nsteps] = moles_CO[i1]
            moles_CH4[i1:nsteps] = moles_CH4[i1]
            moles_H2[i1:nsteps] = moles_H2[i1]
            moles_O2[i1:nsteps] = moles_O2[i1]
            moles_O2_ash[i1:nsteps] = moles_O2_ash[i1]
            moles_O2_lava[i1:nsteps] = moles_O2_lava[i1]
            moles_O2_magma[i1:nsteps] = moles_O2_magma[i1]
            moles_O2_nt[i1:nsteps] = moles_O2_nt[i1]
            moles_O2_red[i1:nsteps] = moles_O2_nt[i1]
            
            return(t,i1,Pres,moles_atm,Rmm_atm,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,moles_O2,dz,T_surf,moles_O2_ash,moles_O2_lava,moles_O2_magma,moles_O2_nt,moles_O2_red,z_meltlayer_tot,z_max)
            
            break
        
        # convert column H2O in atmosphere to Earth oceans to calculate surface temperature using Turbet et al. 2020
        oceans[i1] = moles_H2O[i1-1]/colEO
        

        if oceans[i1] > 1.e-5:
            # Find surface temperature as a function of atmospheric water
            T_surf[i1],z_melt[i1] = surface_T(oceans[i1],Pres[i1-1],t[i1],moles_H2O[i1-1],moles_atm[i1-1])
            
            if T_surf[i1] < 1283.3:
                endmagma = 1.
            if T_surf[i1] > 1283.3:
                endmagma = 0.
        elif oceans[i1] <= 1.e-5:
            z_melt[i1] = 0.
            T_surf[i1] = T_surf[i1-1]
        
        # artificially turn off surface melting - this is necessary for no melt runs!
        z_melt[i1] = 0.
        
        # stop calculating depth of melt layer beyond 50000 m (sufficient to remove O2 from 1 Earth Ocean through oxidation, i.e. 3x max model H2O input)
        if z_meltlayer_tot[i1] >= 50000:
            z_melt[i1] = 0.
            endmagma = 1.
            
            # re-calculated total erupted volume needed to erupt CO2 from end of surface melting to end of model
            nn,m = eruption_rate(nsteps,moles_CO2[i1-1],frac_ext,P_final,mCO2tot,age_sun,(t[i1]/secingyr),b_err,dzdt_pres,frac_volc)
            
        
        if z_melt[i1] > 0. and endmagma == 0.:
            degas = 0.
            dz[i1] = 0.
            moles_CO2[i1] = moles_CO2[i1-1]
            moles_H2O[i1] = moles_H2O[i1-1]
            
            #break into small timestep loops to deal with surface melting
            moles_H2O[i1],moles_O2[i1],tf,T_surf[i1],meltfrac0,moles_O2_ash[i1],moles_O2_lava[i1],dmoles_O2_magma,moles_O2_nt[i1],moles_O2_red[i1],T_melt,latent,ox_frac0,z_meltlayer_tot[i1],endmagma = meltloss(t[i1-1],t[i1],moles_H2O[i1],moles_CO2[i1],moles_CO[i1],moles_CH4[i1],moles_H2[i1],
                                                                                  moles_O2[i1],T_atm,kb,q,meltfrac0,epsilon,T_melt,latent,ox_frac0,z_meltlayer_tot[i1-1],endmagma)
                                                                                  
                                                                                 
            # Update atmospheric O2 after losses to surface melt layer
            moles_O2_magma[i1] = moles_O2_magma[i1-1] + dmoles_O2_magma

            t[i1] = tf

            dz_check = m*np.exp(-b_err*t_sun) + nn
            total_erupt = np.sum(dz_check[1:nsteps]*((t[1:nsteps]-t[0:nsteps-1])/secingyr))

            if total_erupt > z_max:
                
                print('ERROR: z_erupt too large, abort run',total_erupt)
                
                moles_atm[i1:nsteps] = moles_atm[i1]
                moles_H2O[i1:nsteps] = moles_H2O[i1]
                moles_CO2[i1:nsteps] = moles_CO2[i1]
                moles_CO[i1:nsteps] = moles_CO[i1]
                moles_CH4[i1:nsteps] = moles_CH4[i1]
                moles_H2[i1:nsteps] = moles_H2[i1]
                moles_O2[i1:nsteps] = moles_O2[i1]
                moles_O2_ash[i1:nsteps] = moles_O2_ash[i1]
                moles_O2_lava[i1:nsteps] = moles_O2_lava[i1]
                moles_O2_magma[i1:nsteps] = moles_O2_magma[i1]
                moles_O2_nt[i1:nsteps] = moles_O2_nt[i1]
                moles_O2_red[i1:nsteps] = moles_O2_red[i1]
                
                return(t,i1,Pres,moles_atm,Rmm_atm,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,moles_O2,dz,T_surf,moles_O2_ash,moles_O2_lava,moles_O2_magma,moles_O2_nt,moles_O2_red,z_meltlayer_tot,z_max)
                
                break
        
        # Switch from surface melting sequence to regular model run
        if z_melt[i1] == 0. or endmagma == 1.:
            
            rk1_CO2 = 0.
            rk1_CO = 0.
            rk1_CH4 = 0.
            rk1_H2 = 0.
            rk1_H2O = 0.
            rk1_O2 = 0.
            rk1_moles_O2_ash = 0.
            rk1_moles_O2_lava = 0.
            rk1_moles_O2_nt = 0.
            rk1_moles_O2_red = 0.

            degas = 1

            if np.mean(meltfrac0)*z_end >= 50000: #track atmosphere but DON'T add volatiles through volcanism
                degas = 0.
                

            # Set up values for RK4 integration if on 1st step of run
            if i1 == 1:
                
                rk1_CO2 = 0.
                rk1_CO = 0.
                rk1_CH4 = 0.
                rk1_H2 = 0.
                rk1_H2O = 0.
                rk1_O2 = 0.
                rk1_moles_O2_ash = 0.
                rk1_moles_O2_lava = 0.
                rk1_moles_O2_nt = 0.
                rk1_moles_O2_red = 0.
                
                
            if i1 > 1: 
                rk1_CO2 = moles_CO2[i1-1] - moles_CO2[i1-2]
                rk1_CO = moles_CO[i1-1] - moles_CO[i1-2]
                rk1_CH4 = moles_CH4[i1-1] - moles_CH4[i1-2]
                rk1_H2 = moles_H2[i1-1] - moles_H2[i1-2]
                rk1_H2O = moles_H2O[i1-1] - moles_H2O[i1-2]
                rk1_O2 = moles_O2[i1-1] - moles_O2[i1-2]
                rk1_moles_O2_ash = moles_O2_ash[i1-1] - moles_O2_ash[i1-2]
                rk1_moles_O2_lava = moles_O2_lava[i1-1] - moles_O2_lava[i1-2]
                rk1_moles_O2_nt = moles_O2_nt[i1-1] - moles_O2_nt[i1-2]
                rk1_moles_O2_red = moles_O2_red[i1-1] - moles_O2_red[i1-2]
                
                
            # RK4 integration
            thalf = t[i1-1] + 0.5*(t[i1] - t[i1-1])
            trat = (t[i1]-t[i1-1])/(thalf-t[i1-1])
        
            rk2_CO2,rk2_CO,rk2_CH4,rk2_H2,rk2_H2O,rk2_O2,rk2_Pres,rk2_moles_O2_ash,rk2_moles_O2_lava,rk2_moles_O2_nt,rk2_moles_O2_red = stdrk4(degas,frac_ext,b_err,FH_ref[i1],T,T_atm,FMQ,mCO2tot,mH2Otot,(moles_CO2[i1-1] + (0.5*rk1_CO2)),(moles_CO[i1-1] + (0.5*rk1_CO)),(moles_CH4[i1-1] + (0.5*rk1_CH4)),(moles_H2[i1-1] + (0.5*rk1_H2)),(moles_H2O[i1-1] + (0.5*rk1_H2O)),(moles_O2[i1-1] + (0.5*rk1_O2)),t[i1-1],thalf,tstart,dz[i1],sH2O,sCO2,CO_CO2_rat,CH4_CO2_rat,H2_H2O_rat)
            rk3_CO2,rk3_CO,rk3_CH4,rk3_H2,rk3_H2O,rk3_O2,rk3_Pres,rk3_moles_O2_ash,rk3_moles_O2_lava,rk3_moles_O2_nt,rk3_moles_O2_red = stdrk4(degas,frac_ext,b_err,FH_ref[i1],T,T_atm,FMQ,mCO2tot,mH2Otot,(moles_CO2[i1-1] + (0.5*rk2_CO2)),(moles_CO[i1-1] + (0.5*rk2_CO)),(moles_CH4[i1-1] + (0.5*rk2_CH4)),(moles_H2[i1-1] + (0.5*rk2_H2)),(moles_H2O[i1-1] + (0.5*rk2_H2O)),(moles_O2[i1-1] + (0.5*rk2_O2)),t[i1-1],thalf,tstart,dz[i1],sH2O,sCO2,CO_CO2_rat,CH4_CO2_rat,H2_H2O_rat)
            rk4_CO2,rk4_CO,rk4_CH4,rk4_H2,rk4_H2O,rk4_O2,rk4_Pres,rk4_moles_O2_ash,rk4_moles_O2_lava,rk4_moles_O2_nt,rk4_moles_O2_red = stdrk4(degas,frac_ext,b_err,FH_ref[i1],T,T_atm,FMQ,mCO2tot,mH2Otot,(moles_CO2[i1-1] + rk3_CO2),(moles_CO[i1-1] + rk3_CO),(moles_CH4[i1-1] + rk3_CH4),(moles_H2[i1-1] + rk3_H2),(moles_H2O[i1-1] + rk3_H2O),(moles_O2[i1-1] + rk3_O2),t[i1-1],t[i1],tstart,dz[i1],sH2O,sCO2,CO_CO2_rat,CH4_CO2_rat,H2_H2O_rat)
            
            moles_CO2[i1] = moles_CO2[i1-1] + ((1/6) * (rk1_CO2 + (2*rk2_CO2*trat) + (2*rk3_CO2*trat) + rk4_CO2))
            moles_CO[i1] = moles_CO[i1-1] + ((1/6) * (rk1_CO + (2*rk2_CO*trat) + (2*rk3_CO*trat) + rk4_CO))
            moles_CH4[i1] = moles_CH4[i1-1] + ((1/6) * (rk1_CH4 + (2*rk2_CH4*trat) + (2*rk3_CH4*trat) + rk4_CH4))
            moles_H2[i1] = moles_H2[i1-1] + ((1/6) * (rk1_H2 + (2*rk2_H2*trat) + (2*rk3_H2*trat) + rk4_H2))
            moles_H2O[i1] = moles_H2O[i1-1] + ((1/6) * (rk1_H2O + (2*rk2_H2O*trat) + (2*rk3_H2O*trat) + rk4_H2O))
            moles_O2[i1] = moles_O2[i1-1] + ((1/6) * (rk1_O2 + (2*rk2_O2*trat) + (2*rk3_O2*trat) + rk4_O2))
            
            moles_O2_ash[i1] = moles_O2_ash[i1-1] + ((1/6) * (rk1_moles_O2_ash + (2*rk2_moles_O2_ash*trat) + (2*rk3_moles_O2_ash*trat) + rk4_moles_O2_ash))
            moles_O2_lava[i1] = moles_O2_lava[i1-1] + ((1/6) * (rk1_moles_O2_lava + (2*rk2_moles_O2_lava*trat) + (2*rk3_moles_O2_lava*trat) + rk4_moles_O2_lava))
            moles_O2_nt[i1] = moles_O2_nt[i1-1] + ((1/6) * (rk1_moles_O2_nt + (2*rk2_moles_O2_nt*trat) + (2*rk3_moles_O2_nt*trat) + rk4_moles_O2_nt))
            moles_O2_red[i1] = moles_O2_red[i1-1] + ((1/6) * (rk1_moles_O2_red + (2*rk2_moles_O2_red*trat) + (2*rk3_moles_O2_red*trat) + rk4_moles_O2_red))
            
            
            moles_O2_magma[i1] = moles_O2_magma[i1-1]
            
        # Update atmosphere before next timestep
        moles_atm[i1] = moles_CO2[i1] + moles_CO[i1] + moles_H2[i1] + moles_H2O[i1] + moles_O2[i1] + moles_CH4[i1]
        Rmm_atm[i1] = (moles_CO2[i1]*Rmm_CO2 + moles_H2O[i1]*Rmm_H2O + moles_O2[i1]*Rmm_O2 + moles_CO[i1]*Rmm_CO + moles_H2[i1]*Rmm_H2 + moles_CH4[i1]*Rmm_CH4)/moles_atm[i1]
        
        # Update total pressure
        Pres[i1] = moles_atm[i1] * Rmm_atm[i1] * g
        # print('Pres',Pres[i1]/1e5,'moles_atm',moles_atm[i1],'Rmm_atm',Rmm_atm[i1])

        # Save partial run if model is close to maximum allowed run time        
        check_time = time.time()
        t_diff = check_time - start_time
        if t_diff >= threshold_time:
            print('time is up, lets save our progress!')
            return(t,i1,Pres,moles_atm,Rmm_atm,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,moles_O2,dz,T_surf,moles_O2_ash,moles_O2_lava,moles_O2_magma,moles_O2_nt,moles_O2_red,z_meltlayer_tot)
            
        # Quit when final present day moles CO2 in atmosphere reached
        if metamorphic_degas == 0:
            if moles_CO2[i1] >= moles_final:
                print('its over!')
                moles_atm[i1:nsteps] = moles_atm[i1]
                moles_H2O[i1:nsteps] = moles_H2O[i1]
                moles_CO2[i1:nsteps] = moles_CO2[i1]
                moles_CO[i1:nsteps] = moles_CO[i1]
                moles_H2[i1:nsteps] = moles_H2[i1]
                moles_O2[i1:nsteps] = moles_O2[i1]
                
                moles_O2_ash[i1:nsteps] = moles_O2_ash[i1]
                moles_O2_lava[i1:nsteps] = moles_O2_lava[i1]
                moles_O2_magma[i1:nsteps] = moles_O2_magma[i1]
                moles_O2_nt[i1:nsteps] = moles_O2_nt[i1]
                moles_O2_red[i1:nsteps] = moles_O2_red[i1]
                
                return(t,i1,Pres,moles_atm,Rmm_atm,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,moles_O2,dz,T_surf,moles_O2_ash,moles_O2_lava,moles_O2_magma,moles_O2_nt,moles_O2_red,z_meltlayer_tot,z_max)
                
                break
                
            
    return(t,i1,Pres,moles_atm,Rmm_atm,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,moles_O2,dz,T_surf,moles_O2_ash,moles_O2_lava,moles_O2_magma,moles_O2_nt,moles_O2_red,z_meltlayer_tot,z_max)

#%%
def timestep(tstart,t_end,nsteps,epsilon,magma):
    # Calculate timestep sizes for model depending on whether surface melting is occuring or not, and find XUV flux at all timesteps
    tau_sun = 0.1 #gyrs #gyrs (100 Myr after CAIs), end of habitable era

    if magma == 1:
        nsteps = nsteps
        t = secingyr*np.linspace(tstart,t_end,nsteps)
    elif magma == 0:
        t = secingyr*np.linspace(tstart,t_end,nsteps)
    
    t_sun = t/(1.e9 * (365.*60.*60.*24.)) #use for calculating XUV
    beta = 1.24 #exponent for Sun evolution used in Zahnle & Casting Cosmic Shorelines 2017

    #calculate XUV flux at each timestep:
    Fxuv_venus = 1.6e30*((t_sun/tau_sun)**(-beta)) * (1.e-7/A_venus)/(a_venus**2) #factor of -7 to convert from ergs to J, W m^-2 s^-1, A_Venus to give energy/time/unit area
    
    # Reference H escape rate (from Tian 2015)
    FH_ref = (epsilon * (R_venus + 200.e3)**3 * Fxuv_venus)/(4*math.pi*G*M_venus)
    
    return(t,t_sun,Fxuv_venus,FH_ref)

#%%
def eruption_rate(nsteps,moles_CO2,frac_ext,P_final,mCO2tot,t_sun,tstart,b_err,dzdt_pres,frac_volc):
    # Calculate eruption rate required to degas all of remaining CO2 to reach present day atmosphere and match modern Venus eruption rates
    # Assumes eruption rate decreases exponentially over time 
    
    M_CO20 = P_final*(1-frac_volc)/g
    M_CO21 = P_final/g
    
    z_lava = (M_CO21 - M_CO20)/(mCO2tot*rho)
    

    # For eruption rate purposes, t_0 is start of model run, so adjust time accordingly
    te = 4.5 - tstart

    if frac_ext < 1:
        dzt_p = (1.e9*dzdt_pres/A_venus)/frac_ext #m per Gyr erupted per unit surface area
    elif frac_ext == 1:
        dzt_p = (1.e9*dzdt_pres/A_venus)
    
    #Parameters for exponential function dz(t) = m*exp(-b*t) + nn
    nn = 0.
    m = (dzt_p - (z_lava*-b_err))/np.exp(-b_err*tstart)
    
    return(nn,m)
    
#%%
def surface_T(oceans,Pres,t,moles_H2O,moles_atm):
    
    # Calculate surface temperature based on atmospheric water vapor during runaway greenhouse
    T_ep = 80.*(math.log10(oceans)**2.) + 700.*(math.log10(oceans)) + 1700.
    
    if T_ep > 800.:
        fL = L(t/secingyr)
        flux = fL.item()
        
        if oceans >= 0.1:
            
            P_H2O = 1e-5* Pres * (moles_H2O/moles_atm)
            x_T = (math.log10(P_H2O) - k1)/k2
            z_T = (math.log10(flux) - k5)/k6
            
            logT = c1 + c2*x_T + c3*y_T + c4*z_T + c5*x_T**2 + c6*x_T*y_T + c7*y_T**2 + c8*z_T**2 + c9*y_T**3 + c10*z_T**3
            
            T_surf = 10**(logT)
            
        # Backup calculation - interpolates between points extracted from graph used to find logT equation in Turbet 2020
        else:
            Oc_arr = np.log10([0.1,1,10])
            T_arr = np.zeros(3)
            for ib in range(0,np.size(T_arr)):
                T_arr[ib] = rbfi(flux,Oc_arr[ib])
                
            T_fit = np.polyfit(Oc_arr,T_arr,deg=1)

            T_surf = T_fit[1] + T_fit[0]*math.log10(oceans)
    else:
        T_surf = 700. # Surface temperature not explicitly calculated outside of surface melting regime

    if T_surf < 1323.3:
        test_melt = 0.
    elif T_surf > 1323.3 and T_surf<3000.0:
        test_melt = melt_function(T_surf)
        # print('T_surf',T_surf)
    elif T_surf >= 3000.0:
        test_melt = 1.0
    
    z_melt = 0.
    if test_melt <= 0.:
        test_melt = 0.
        z_melt = 0.
    elif test_melt > 0.:
        z_melt = 1.
        
    return(T_surf,z_melt)

#%%
def degassing(Pres,T,FMQ,mCO2tot,mH2Otot,moles_CO20,moles_CO0,moles_CH40,moles_H20,moles_H2O0,t1,t2,dz,sH2O,sCO2,CO_CO2_rat,CH4_CO2_rat,H2_H2O_rat):
    # Calculate H2O and CO2 degassed at each timestep using Wogan et al. 2020
    # This function REQUIRES VolcGasses
    
    P = Pres/1.e5 # bar
    solH2O = np.interp(P,Pgrid,sH2O)
    solCO2 = np.interp(P,Pgrid,sCO2)
    CO_CO2 = np.interp(P,Pgrid,CO_CO2_rat)
    CH4_CO2 = np.interp(P,Pgrid,CH4_CO2_rat)
    H2_H2O = np.interp(P,Pgrid,H2_H2O_rat)
    
    
    # Find exsolved mass fraction of CO2 and H2O
    dCO2 = mCO2tot - solCO2
    dH2O = mH2Otot - solH2O
    
    # print('solCO2',solCO2,'mCO2tot',mCO2tot)
    
    
    if dH2O < 0:
        dH2O = 0.
    if dCO2 < 0:
        dCO2 = 0.
    
    # Volume fraction of gas in erupting material 
    # First calculate moles of each gas
    mol_CO2_degas = rho* dCO2/Rmm_CO2 * (1-(CO_CO2+CH4_CO2))
    mol_CO_degas = rho* dCO2/Rmm_CO2 * CO_CO2
    mol_CH4_degas = rho* dCO2/Rmm_CO2 * CH4_CO2
    
    
    mol_H2_degas = rho* dH2O/Rmm_H2O * H2_H2O
    mol_H2O_degas = rho* dH2O/Rmm_H2O * (1-H2_H2O)
    mol_tot_degas = mol_CO2_degas + mol_H2O_degas + mol_CO_degas + mol_H2_degas
    
    Vdegas = (mol_tot_degas*R*T)/Pres
    
    # Find volume fraction gas in erupted material
    V_frac = Vdegas/(Vdegas + 1.)
    
    # Find total erupted thickness for timestep (note that 'dz' input is already midpoint between dz[t1] and dz[t2] from main loop)
    dz_C = dz*((t2-t1)/(secingyr))
    
    # print('dz',dz)

    # Calculate degassed moles CO2 and H2O before atmospheric escape
    moles_CO2 = moles_CO20 + (dz_C*mol_CO2_degas)
    moles_CO = moles_CO0 + (dz_C*mol_CO_degas)
    moles_CH4 = moles_CH40 + (dz_C*mol_CH4_degas)
    moles_H2 = moles_H20 + (dz_C*mol_H2_degas)
    moles_H2O = moles_H2O0 + (dz_C*mol_H2O_degas)

    return(moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,V_frac)

#%%
def Hloss(moles_H2O,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_O2,T_atm,FH_ref,t1,t2,kb,Rmm_atm):
    # Calculate XUV and diffusion limited H escape
    
    # Calculate background atmosphere properties
    moles_atm = moles_H2O + moles_CO2 + moles_O2 + moles_CO + moles_CH4 + 0.5*moles_H2O + moles_H2
    moles_atm_heavy = 0.5*moles_H2O + moles_CO2 + moles_O2 + moles_CO + moles_CH4
    mean_matm_heavy = (Rmm_O2*(0.5*moles_H2O + moles_O2) + Rmm_CO2*moles_CO2 + Rmm_CO*moles_CO + Rmm_CH4*moles_CH4)/moles_atm_heavy
    
    # Mass ratio of H and O for calculating whether O dragged along
    q = Rmm_O/Rmm_H
    
    # Mixing ratios of H and O
    mixing_rat_H = (2.*moles_H2O + 2*moles_H2)/(moles_atm_heavy + (2.*moles_H2O + 2*moles_H2))
    mixing_rat_O = (2.*moles_O2)/(moles_atm_heavy + (2.*moles_H2O))
    M_Oloss_column = 0.
    
    # Diffusion coefficient for O in H
    b =   100. * 4.8e17 * (T_atm**0.75) #m^-1 s^-1, from Tian 2015
    
    # Calculate parameters for XUV limited escape
    if mixing_rat_H>1.e-6:

        # Calculate critical mass that can be dragged along by escaping H
        mc = (((mixing_rat_H + mixing_rat_O*(Rmm_O/Rmm_H)**2)/(mixing_rat_H + mixing_rat_O*(Rmm_O/Rmm_H)))*(Rmm_H/Avo)) + ((((kb*T_atm*FH_ref)/(b*g)))/(mixing_rat_H + mixing_rat_O*(Rmm_O/Rmm_H)))
        
        # If mc > mO, treat atmosphere as mixture of O and H only and find O dragged along
        if mc >= (Rmm_O/Avo):
            print('ESCAPING O')
            
            Y = (kb*T_atm*FH_ref)/((q-1.)*b*g*mixing_rat_H*(Rmm_H/Avo))
            eta = (Y-1.)/(Y * (mixing_rat_H/(1. - mixing_rat_H)) + q)
            
            M_Hloss_column_XUV = (t2-t1) *(1./(1.+(eta*q)))*FH_ref
            M_Oloss_column = (t2-t1) * (Rmm_O/Rmm_H)*(eta/(1.+(eta*q)))*FH_ref
            
        elif mc < (Rmm_O/Avo):
            # If only H escaping, find XUV limited H flux
            M_Hloss_column_XUV = (t2-t1) * FH_ref

        # Binary diffusion of H in background atmosphere:
        fCO2 = moles_CO2/moles_atm
        fCO = moles_CO2/moles_atm
        fCH4 = moles_CH4/moles_atm
        fO = mixing_rat_O
        
        # Combine binary diffusion coefficients for H in background CO2, CO, and O
        b_H_CO2 = 100 * 8.4e17 * T_atm**0.6
        b_H_CO = 100 * 6.5e17 * T_atm**0.7
        b_H_O = b
        b_H_atm = (b_H_CO2 * fCO2 + b_H_CO * fCO + b_H_O*fO) / (fCO2 + fCO + fO)

        # Calculate diffusion limited escape
        M_Hloss_column_diff = (t2-t1) * (Rmm_H/Avo)* mixing_rat_H * ((G * M_venus )/(R_venus**2)) * ((mean_matm_heavy - Rmm_H)/Avo) * (b_H_atm/(kb*T_atm))
        
        # Compare diffusion and XUV limited escape to implement lowest limit
        if M_Hloss_column_diff < M_Hloss_column_XUV:
            M_Hloss_column = M_Hloss_column_diff
            M_Oloss_column = 0.

        elif M_Hloss_column_diff >= M_Hloss_column_XUV:
            M_Hloss_column = M_Hloss_column_XUV

    else:
        M_Hloss_column = 0.
       
    # H2O loss in moles 
    moles_H2loss = (M_Hloss_column/Rmm_H)*0.5 #equivalent moles of H2 lost through H escape
    
    
    # Update atmosphere variables and ensure nothing can drop below 0

    if moles_H2loss >= moles_H2:
        moles_H2Oloss = moles_H2loss - moles_H2
        moles_H2 = 0.
    elif moles_H2loss < moles_H2:
        moles_H2 = moles_H2 - moles_H2loss
        moles_H2Oloss = 0.
  
    if moles_H2Oloss >= moles_H2O:
        moles_O2gain = moles_H2O*0.5 - (0.5*M_Oloss_column/Rmm_O)
        moles_H2Oloss = moles_H2O
        moles_H2O = 0.
    else:
        moles_H2O = moles_H2O - moles_H2Oloss
        moles_O2gain = moles_H2Oloss*0.5 - (0.5*M_Oloss_column/Rmm_O)
        
    moles_O2 = moles_O2 + moles_O2gain
    if moles_O2 < 0.:
        moles_O2gain = -moles_O2
        moles_O2 = 0.
    
    return(moles_H2O,moles_O2,moles_H2)

#%%
def meltloss(t1,t2,moles_H2O0,moles_CO20,moles_CO0,moles_CH40,moles_H20,moles_O20,T_atm,kb,q,meltfrac0,epsilon,T,latent,ox_frac0,z_meltlayer_tot0,endmagma): #,rbfi,solidus,gill_t,gill_molesO2):
    # Timestepping sequence for case where surface melting is occuring 
    
    degas = 0
    
    magma_yn = 1.
    
    # Calculate stable timestep size
    nsteps = 5*int((16/5)*((t2 - t1)/secingyr)/(1e-6))
    
    t,t_sun,Fxuv_venus,FH_ref = timestep((t1/secingyr),(t2/secingyr),nsteps,epsilon,magma_yn)

    Pres = np.zeros(nsteps)
        
    # Set up arrays to track atmosphere
    moles_CO2 = np.zeros(nsteps)
    moles_CO2[0] = moles_CO20
    moles_CO = np.zeros(nsteps)
    moles_CO[0] = moles_CO0
    moles_CH4 = np.zeros(nsteps)
    moles_CH4[0] = moles_CH40
    moles_H2 = np.zeros(nsteps)
    moles_H2[0] = moles_H20
    moles_H2O = np.zeros(nsteps)
    moles_H2O[0] = moles_H2O0
    oceans = np.zeros(nsteps)
    moles_O2 = np.zeros(nsteps)
    moles_O2[0] = moles_O20
    moles_atm = np.zeros(nsteps)
    moles_atm[0] = moles_CO2[0] + moles_H2O[0] + moles_O2[0] + moles_CO + moles_H2
    
    # Set up arrays to track O loss mechanisms
    moles_O2_magma = np.zeros(nsteps)
    moles_O2_nt = np.zeros(nsteps)
    moles_O2_ash = np.zeros(nsteps)
    moles_O2_lava = np.zeros(nsteps)
    moles_O2_red = np.zeros(nsteps)
    O2_ox_melt = np.zeros(nsteps)
    
    # Set up arrays to track surface melting variables
    T_surf = np.zeros(nsteps)
    z_meltlayer = np.zeros(nsteps)
    
    Rmm_atm = np.zeros(nsteps)
    Rmm_atm[0] = (moles_CO2[0]*Rmm_CO2 + moles_H2O[0]*Rmm_H2O + moles_O2[0]*Rmm_O2 + moles_CO[0]*Rmm_CO + moles_H2[0]*Rmm_H2 + moles_CH4*Rmm_CH4)/moles_atm[0]
    Pres[0] = moles_atm[0] * Rmm_atm[0] * g
    
    # 1D heat diffusion in crust and melting
    dz = z[1]-z[0]
    
    if T[100]== 0:
        
        T = T_bg*np.ones(nz)  
        T = T + q*z/k
        
    meltfrac = np.zeros(nz)

    meltfrac_out = meltfrac0
    ox_frac = ox_frac0
    z_meltlayer_tot =  z_meltlayer_tot0
    
    for i1 in range(1,nsteps):

        oceans[i1] = moles_H2O[i1-1]/colEO

        if oceans[i1] <= 1.e-5:
            T_surf[i1] = 700.
            break
        
        # Find surface temperature as a function of water in atmosphere
        T_surf[i1],z_melt = surface_T(oceans[i1],Pres[i1-1],t[i1],moles_H2O[i1-1],moles_atm[i1-1])
        
        # Set temperature at top of column of crust to surface temperature for diffusion calculation
        T[0] = T_surf[i1]
        
        # Terminate surface melting timestep sequence when surface above basalt solidus
        if T_surf[i1] < 1323.3:
            z_meltlayer[i1] = 0

            # Update variables for return to main script
            endmagma = 1.
            moles_H2O_out = moles_H2O[i1-1]
            moles_O2_out = moles_O2[i1-1]    
            ash = np.sum(moles_O2_ash)
            lava = np.sum(moles_O2_lava)
            magma = np.sum(moles_O2_magma)
            nonthermal = np.sum(moles_O2_nt)
            reducing = np.sum(moles_O2_red)
            
            return(moles_H2O_out,moles_O2_out,t[i1],T_surf[i1],meltfrac_out,ash,lava,magma,nonthermal,reducing,T,latent,ox_frac,z_meltlayer_tot0,endmagma)

        if z_meltlayer_tot < (0.99*z_end):

            dTdz = (1/dz) * matrix.dot(T)
            dTdz = dTdz - (latent/k)
            dTdz[nz-1] = dTdz[nz-1] + q/k 
            T_new = T + np.multiply((kappa * (t[i1]-t[i1-1])/ (dz)), dTdz)
            T=T_new
        
        # Terminate run if surface temperature is too high - usually a numerical stability error generated by timestep size
        if np.max(T) > 3500.:
            print('ERROR: UNSTABLE','max T',np.max(T))
            exit()

        # Calculate surface melt fraction to determine if surface melting still occuring
        mf1 = melt_function(T[0])
        test_melt = mf1.item() 
        
        # If no melting, terminate surface melting sequence
        if test_melt <= 0.01:
            
            z_meltlayer[i1] = 0
            print('end of surface melting, time',(t1/secingyr)*1e9,'T_surf',T_surf[i1])
            
            endmagma = 1.
            moles_H2O_out = moles_H2O[i1-1]
            moles_O2_out = moles_O2[i1-1]    
            ash = np.sum(moles_O2_ash)
            lava = np.sum(moles_O2_lava)
            magma = np.sum(moles_O2_magma)
            nonthermal = np.sum(moles_O2_nt)
            reducing = np.sum(moles_O2_red)
            
            return(moles_H2O_out,moles_O2_out,t[i1],T_surf[i1],meltfrac_out,ash,lava,magma,nonthermal,reducing,T,latent,ox_frac,z_meltlayer_tot0,endmagma)
            
        # If melting, continue surface melting sequence
        elif test_melt > 0.01:
            
            endmagma = 0.
            
            if i1 % (nsteps/2) == 0:
                # Determine depths at which melting is occuring
                for i2 in range(0,nz-1):
                    if T[i2] < 1200.:
                        test_melt = 0.
                        break
                    elif T[i2] >= 1200.:
                        test_melt = melt_function(T[i2])
        
                    if test_melt <= 0.:
                        test_melt = 0.
                        break
                    if test_melt >= 1.:
                        test_melt = 1.
                        
                    if i2 == 0 and test_melt < 0.01:
                        return(moles_H2O_out,moles_O2_out,t[i1],T_surf[i1],meltfrac_out,ash,lava,magma,nonthermal,reducing,T,latent,ox_frac,z_meltlayer_tot0,endmagma)
                    
                    # Find cumulative melt fraction 
                    cummelt = meltfrac_out[i2]
                    
                    # Find melt fraction added at timestep i1
                    if cummelt > 1:
                        exit('Error: cumulative melt fraction >1')
                    
                    if cummelt == 1.:
                        meltfrac[i2] = 0.
                        
                    elif test_melt - cummelt <= 0.:
                        meltfrac[i2] = 0.
                        
                    elif test_melt - cummelt > 0.:
                        meltfrac[i2] = test_melt - cummelt
                        
                    meltfrac_out[i2] = cummelt + meltfrac[i2]

            else:   
                meltfrac = meltfrac*0.

            melt_ox = meltfrac - ox_frac
            
            for i2 in range(0,nz-1):
                
                if melt_ox[i2]<= 0.:
                    melt_ox[i2] = 0.
            
            # Track cumulative thickness of surface melt layer at each timestep
            z_meltlayer[i1] = np.mean(melt_ox)*(z_end-dz)
            z_meltlayer_tot =  z_meltlayer_tot + z_meltlayer[i1]
            
            # Latent heat consumed/released at timestep
            if i1 % (nsteps/2) == 0:
                latent = meltfrac * (z_end/nz) * lat_melt * rho/((t2-t1)/(nsteps/i1))

            if z_meltlayer_tot>(z_end-dz):
                print('melt layer too thick!',np.sum(z_meltlayer[0:i1]),'melt layer',z_meltlayer[0:i1])
                
                exit()
            
            # Calculate O loss to melt oxidation
            O2_ox_melt[i1] = (0.5*rho*z_meltlayer[i1]*max_water)/Rmm_H2O
            
        # Set up variables for RK4 integration
        if i1 == 1:
            rk1_CO2 = 0.
            rk1_H2O = 0.
            rk1_CO = 0.
            rk1_CH4 = 0.
            rk1_H2 = 0.
            rk1_O2 = 0.
            rk1_moles_O2_nt = 0.
        elif i1 > 1: 
            rk1_CO2 = moles_CO2[i1-1] - moles_CO2[i1-2]
            rk1_CO2 = moles_CO[i1-1] - moles_CO[i1-2]
            rk1_CH4 = moles_CH4[i1-1] - moles_CH4[i1-2]
            rk1_H2O = moles_H2[i1-1] - moles_H2[i1-2]
            rk1_H2O = moles_H2O[i1-1] - moles_H2O[i1-2]
            rk1_O2 = moles_O2[i1-1] - moles_O2[i1-2]
            rk1_moles_O2_nt = moles_O2_nt[i1-1] - moles_O2_nt[i1-2]
        
        thalf = t[i1-1] + 0.5*(t[i1] - t[i1-1])
        trat = (t[i1]-t[i1-1])/(thalf-t[i1-1])
        
        rk2_CO2,rk2_CO,rk2_CH4,rk2_H2,rk2_H2O,rk2_O2,rk2_Pres,rk2_moles_O2_ash,rk2_moles_O2_lava,rk2_moles_O2_nt,rk2_moles_O2_red = stdrk4(degas,0,0,FH_ref[i1],T,T_atm,0,0,0,(moles_CO2[i1-1] + (0.5*rk1_CO2)),(moles_CO[i1-1] + (0.5*rk1_CO)),(moles_CH4[i1-1] + (0.5*rk1_CH4)),(moles_H2[i1-1] + (0.5*rk1_H2)),(moles_H2O[i1-1] + (0.5*rk1_H2O)),(moles_O2[i1-1] + (0.5*rk1_O2)),t[i1-1],thalf,t1,0,0,0,0,0)
        rk3_CO2,rk3_CO,rk3_CH4,rk3_H2,rk3_H2O,rk3_O2,rk3_Pres,rk3_moles_O2_ash,rk3_moles_O2_lava,rk3_moles_O2_nt,rk3_moles_O2_red = stdrk4(degas,0,0,FH_ref[i1],T,T_atm,0,0,0,(moles_CO2[i1-1] + (0.5*rk2_CO2)),(moles_CO[i1-1] + (0.5*rk2_CO)),(moles_CH4[i1-1] + (0.5*rk2_CH4)),(moles_H2[i1-1] + (0.5*rk2_H2)),(moles_H2O[i1-1] + (0.5*rk2_H2O)),(moles_O2[i1-1] + (0.5*rk2_O2)),t[i1-1],thalf,t1,0,0,0,0,0)
        rk4_CO2,rk4_CO,rk4_CH4,rk4_H2,rk4_H2O,rk4_O2,rk4_Pres,rk4_moles_O2_ash,rk4_moles_O2_lava,rk4_moles_O2_nt,rk4_moles_O2_red = stdrk4(degas,0,0,FH_ref[i1],T,T_atm,0,0,0,(moles_CO2[i1-1] + rk3_CO2),(moles_CO[i1-1] + rk3_CO),(moles_CH4[i1-1] + rk3_CH4),(moles_H2[i1-1] + rk3_H2),(moles_H2O[i1-1] + rk3_H2O),(moles_O2[i1-1] + rk3_O2),t[i1-1],t[i1],t1,0,0,0,0,0)

        moles_CO2[i1] = moles_CO2[i1-1] + ((1/6) * (rk1_CO2 + (2*rk2_CO2*trat) + (2*rk3_CO2*trat) + rk4_CO2))
        moles_CO[i1] = moles_CO[i1-1] + ((1/6) * (rk1_CO + (2*rk2_CO*trat) + (2*rk3_CO*trat) + rk4_CO))
        moles_CH4[i1] = moles_CH4[i1-1] + ((1/6) * (rk1_CH4 + (2*rk2_CH4*trat) + (2*rk3_CH4*trat) + rk4_CH4))
        moles_H2[i1] = moles_H2[i1-1] + ((1/6) * (rk1_H2 + (2*rk2_H2*trat) + (2*rk3_H2*trat) + rk4_H2))
        moles_H2O[i1] = moles_H2O[i1-1] + ((1/6) * (rk1_H2O + (2*rk2_H2O*trat) + (2*rk3_H2O*trat) + rk4_H2O))
        moles_O2[i1] = moles_O2[i1-1] + ((1/6) * (rk1_O2 + (2*rk2_O2*trat) + (2*rk3_O2*trat) + rk4_O2))
        
        moles_O2_ash[i1] = 0.
        moles_O2_lava[i1] = 0.
        moles_O2_nt[i1] = moles_O2_nt[i1-1] + ((1/6) * (rk1_moles_O2_nt + (2*rk2_moles_O2_nt*trat) + (2*rk3_moles_O2_nt*trat) + rk4_moles_O2_nt))
        moles_O2_red[i1] = 0.
        
        
        if endmagma == 1:
            break
        
        # If melt layer deeper than 50km, assume enough oxidizeable material present and in contact with atmosphere to remove all O2 at each timestep - assumes mixing in melt
        # A 50km melt layer is able to remove 3000 m GEL water which exceeds initial H2O used in this model
        if z_meltlayer_tot >= (z_end-dz):
            O2_ox_melt[i1] = moles_O2[i1]
            moles_O2[i1] = 0.
            ox_frac = (ox_frac * 0.) + 1.
            print('all oxidising!')
            
        # Find O2 lost to melt layer oxidation at each timestep
        elif z_meltlayer_tot < (z_end-dz):
            if O2_ox_melt[i1] >= moles_O2[i1]:
                moles_O2_magma[i1] = moles_O2[i1]
                moles_O2[i1] = 0.
                
                for i2 in range(0,nz-1):
                    ox_frac_new = (moles_O2[i1]/O2_ox_melt[i1]) * meltfrac[i2]

                    if (ox_frac_new + ox_frac[i2])>= meltfrac_out[i2]:
                        ox_frac_new = (meltfrac_out[i2] - ox_frac[i2])
                        ox_frac[i2] = meltfrac_out[i2]
                    elif (ox_frac_new + ox_frac[i2]) < meltfrac_out[i2]:
                        ox_frac[i2] = ox_frac0[i2] + ox_frac_new

            else:
                moles_O2_magma[i1] = O2_ox_melt[i1]
                moles_O2[i1] = moles_O2[i1] - O2_ox_melt[i1]
                
                for i2 in range(0,nz-1):
                    ox_frac_new = meltfrac[i2]
                    
                    if (ox_frac_new + ox_frac[i2])>= meltfrac_out[i2]:
                        ox_frac_new = (meltfrac_out[i2] - ox_frac[i2])
                        ox_frac[i2] = meltfrac_out[i2]
                    elif (ox_frac_new + ox_frac[i2]) < meltfrac_out[i2]:
                        ox_frac[i2] = ox_frac0[i2] + ox_frac_new
                        
        # Update atmosphere
        moles_atm[i1] = moles_CO2[i1] + moles_H2O[i1] + moles_O2[i1] + moles_CO[i1] + moles_H2[i1] + moles_CH4[i1]
        Rmm_atm[i1] = (moles_CO2[i1]*Rmm_CO2 + moles_H2O[i1]*Rmm_H2O + moles_O2[i1]*Rmm_O2 + moles_CO[i1]*Rmm_CO + moles_H2[i1]*Rmm_H2 + moles_CH4[i1]*Rmm_CH4)/moles_atm[i1]
        
        # Update total pressure
        Pres[i1] = moles_atm[i1] * Rmm_atm[i1] * g
    
    # Prepare outputs for return to main script
    moles_H2O_out = moles_H2O[i1-1]
    moles_O2_out = moles_O2[i1-1]    
    ash = np.sum(moles_O2_ash)
    lava = np.sum(moles_O2_lava)
    magma = np.sum(moles_O2_magma)
    nonthermal = np.sum(moles_O2_nt)
    reducing = np.sum(moles_O2_red)

    return(moles_H2O_out,moles_O2_out,t[i1],T_surf[i1],meltfrac_out,ash,lava,magma,nonthermal,reducing,T,latent,ox_frac,z_meltlayer_tot,endmagma)

#%%
def oxloss(V_frac,frac_ext,dz_C,t1,t2,moles_O2):
    # Calculate atmospheric O2 lost to oxidizing basalt through volcanism
    
    # Check if volcanism explosive
    if V_frac >= 0.75:
        # If explosive, oxidation is 100% efficient
        O2_ox_volc = 0.5*frac_ext*(max_water*rho*dz_C*((t2-t1)/secingyr)/Rmm_H2O)

    # Oxidise lava produced in effusive eruptions based on Gillmann et al. 2020 Methods (1cm thick oxidized layer on 1m average lava flows)
    else:
        diff_z = 0.01*dz_C*((t2-t1)/secingyr)*frac_ext
        O2_ox_volc = (0.5*rho*diff_z*max_water)/Rmm_H2O
    
    # Updates moles O2 remaining in atmosphere after losses
    if O2_ox_volc >= moles_O2:
        moles_O2 = 0.
    else:
        moles_O2 = moles_O2 - O2_ox_volc
        
    return(moles_O2)

#%%
def nonthermal(t1,t2,moles_O2,Rmm_atm,T_atm,moles_H2O,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_atm):
    # Find oxygen lost via non-thermal escape by interpolating from Gillmann 2020 results
    
    O_1 = np.interp((t1/secingyr),gill_t,gill_molesO2)
    O_2 = np.interp((t2/secingyr),gill_t,gill_molesO2)
    dOdt = O_1 - O_2
    
    # O2 cannot be gained from nonthermal O loss
    if dOdt < 0:
        print('Interpolation error from Gillmann 2020 - negative O loss')
        dOdt = 0
    
    # Calculate diffusion limit for O2 loss
    mixing_rat_O = moles_O2/moles_atm
    mean_matm_notO = (moles_CO2*Rmm_CO2 + moles_H2O*Rmm_H2O + moles_CO*Rmm_CO + moles_H2*Rmm_H2 + moles_CH4*Rmm_CH4)/(moles_CO2 + moles_H2O + moles_CO + moles_CH4 + moles_H2)
    b_O_CO2 = 100 *(1.52e18) * math.sqrt(T_atm) * math.sqrt((0.001/Rmm_atm) + (0.001/Rmm_O))
    M_Oloss_column_diff = (t2-t1) * (Rmm_O/Avo)* mixing_rat_O * ((G * M_venus )/(R_venus**2)) * ((mean_matm_notO - Rmm_O)/Avo) * (b_O_CO2/(kb*T_atm))
    moles_Odiff = (M_Oloss_column_diff/Rmm_O)*0.5
    
    # Apply most limiting rate out of nonthermal losses and diffusion-limited loss
    if moles_Odiff < dOdt:
        dOdt = moles_Odiff
    
    # Update moles O2 remaining in atmosphere
    if dOdt >= moles_O2:
        moles_O2 = 0.
    else:
        moles_O2 = moles_O2 - dOdt

    return(moles_O2)

#%%
def stdrk4(degas,frac_ext,b_err,FH_ref,T,T_atm,FMQ,mCO2tot,mH2Otot,moles_CO20,moles_CO0,moles_CH40,moles_H20,moles_H2O0,moles_O20,tstart,tend,tstart_model,dz,sH2O,sCO2,CO_CO2_rat,CH4_CO2_rat,H2_H2O_rat):
    # Compute values for RK4 numerical integration
    
    # Initial atmosphere at beginning of timestep
    moles_atm0 = moles_CO20 + moles_H2O0 + moles_O20 + moles_CO0 + moles_H20
    Rmm_atm0 = (moles_CO20*Rmm_CO2 + moles_H2O0*Rmm_H2O + moles_O20*Rmm_O2 + moles_CO0*Rmm_CO + moles_H20*Rmm_H2)/moles_atm0
    Pres0 = moles_atm0 * Rmm_atm0 * g
    
    moles_O2_red = 0.

    if degas == 1:
        # Include volcanic degassing in contributions to atmospheric species during timestep
        moles_CO2,moles_CO,moles_CH4,moles_H2,moles_H2O,V_frac = degassing(Pres0,T,FMQ,mCO2tot,mH2Otot,moles_CO20,moles_CO0,moles_CH40,moles_H20,moles_H2O0,tstart,tend,dz,sH2O,sCO2,CO_CO2_rat,CH4_CO2_rat,H2_H2O_rat)
        moles_O2 = moles_O20
        
        if moles_O2 < 2*moles_CH4:
            moles_CH4 = moles_CH4 - 0.5*moles_O2
            if moles_CH4 < 0:
                moles_CH4 = 0.
            moles_CO2 = moles_CO2 + 0.5*moles_O2
            moles_H2O = moles_H2O + moles_O2
            moles_O2_red = moles_O2_red + moles_O2
            moles_O2 = 0.
        elif moles_O2 >= 4*moles_CH4:
            moles_O2 = moles_O2 - 2*moles_CH4
            if moles_O2 < 0:
                moles_O2 = 0.
            moles_CO2 = moles_CO2 + moles_CH4
            moles_H2O = moles_H2O + 2*moles_CH4
            moles_O2_red = moles_O2_red + 2*moles_CH4
            moles_CH4 = 0.
        
        if moles_O2 < 0.5*moles_CO:
            moles_CO = moles_CO - 2*moles_O2
            if moles_CO < 0:
                moles_CO = 0.
            moles_CO2 = moles_CO2 + 2*moles_O2
            moles_O2_red = moles_O2_red + moles_O2
            moles_O2 = 0.
        elif moles_O2 >= 0.5*moles_CO:
            moles_O2 = moles_O2 - 0.5*moles_CO
            if moles_O2 < 0:
                moles_O2 = 0.
            moles_CO2 = moles_CO2 + 0.5*moles_CO
            moles_O2_red = moles_O2_red + 0.5*moles_CO
            moles_CO = 0.
            
        if moles_O2 < 0.5 * moles_H2:
            moles_H2 = moles_H2 - 2*moles_O2
            if moles_H2 < 0:
                moles_H2 = 0.
            moles_H2O = moles_H2O + 2*moles_O2
            moles_O2_red = moles_O2_red + moles_O2
            moles_O2 = 0.
        if moles_O2 >= 0.5*moles_H2:
            moles_O2 = moles_O2 - 0.5*moles_H2
            if moles_O2 < 0:
                moles_O2 = 0.
            moles_O2_red = moles_O2_red + 0.5*moles_H2
            moles_H2 = 0.
        
        
        
    elif degas == 0:
        # No volcanic degassing - i.e. surface melt layer present
        moles_CO2 = moles_CO20
        moles_CO = moles_CO0
        moles_CH4 = moles_CH40
        moles_H2 = moles_H20
        moles_H2O = moles_H2O0
        moles_O2 = moles_O20
        V_frac = 0.
    
    
    # Update atmosphere after degassing - oxidize reducing species in atmosphere, starting with CO:
    
    
    moles_atm = moles_CO2 + moles_H2O + moles_O2 + moles_CO + moles_H2 + moles_CH4
    Rmm_atm = (moles_CO2*Rmm_CO2 + moles_H2O*Rmm_H2O + moles_O2*Rmm_O2 + moles_CO*Rmm_CO + moles_H2*Rmm_H2 + moles_CH4*Rmm_CH4)/moles_atm
    
    # Find H loss and O2 accumulated
    moles_H2O,moles_O2,moles_H2 = Hloss(moles_H2O,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_O2,T_atm,FH_ref,tstart,tend,kb,Rmm_atm)
    
    # Update atmosphere after H loss
    moles_atm = moles_CO2 + moles_H2O + moles_O2 + moles_CO + moles_H2 + moles_CH4
    Rmm_atm = (moles_CO2*Rmm_CO2 + moles_H2O*Rmm_H2O + moles_O2*Rmm_O2 + moles_CO*Rmm_CO + moles_H2*Rmm_H2 + moles_CH4*Rmm_CH4)/moles_atm
    
    # O2 losses to oxidizing volcanic products
    moles_O2_volc = oxloss(V_frac,frac_ext,dz,tstart,tend,moles_O2) #find moles of O2 left in atmosphere after volcanic products oxidized
    
    if degas == 0:
        moles_O2_lava = 0
        moles_O2_ash = 0
    elif degas == 1:
        if V_frac >= 0.75:
            # Explosive volcanism occuring - record O losses to ash ox
            moles_O2_ash = moles_O2 - moles_O2_volc
            moles_O2_lava = 0
        elif V_frac < 0.75:
            # Effusive volcanism occuring - record O losses to lava ox
            moles_O2_ash = 0
            moles_O2_lava = moles_O2 - moles_O2_volc
    
    # Update atmosphere
    moles_O2 = moles_O2_volc
    moles_atm = moles_CO2 + moles_H2O + moles_O2 + moles_CO + moles_H2 + moles_CH4
    Rmm_atm = (moles_CO2*Rmm_CO2 + moles_H2O*Rmm_H2O + moles_O2*Rmm_O2 + moles_CO*Rmm_CO + moles_H2*Rmm_H2 + moles_CH4*Rmm_CH4)/moles_atm
    
    # Non-thermal O losses
    moles_O2_space = nonthermal(tstart,tend,moles_O2,Rmm_atm,T_atm,moles_H2O,moles_CO2,moles_CO,moles_CH4,moles_H2,moles_atm)
    moles_O2_nt = moles_O2 - moles_O2_space
    moles_O2 = moles_O2_space
    
    # Update atmosphere before next timestep
    moles_atm = moles_CO2 + moles_H2O + moles_O2 + moles_CO + moles_H2 + moles_CH4
    Rmm_atm = (moles_CO2*Rmm_CO2 + moles_H2O*Rmm_H2O + moles_O2*Rmm_O2 + moles_CO*Rmm_CO + moles_H2*Rmm_H2 + moles_CH4*Rmm_CH4)/moles_atm
    Pres = moles_atm * Rmm_atm * g
    dmoles_CO2 = moles_CO2 - moles_CO20
    dmoles_CO = moles_CO - moles_CO0
    dmoles_CH4 = moles_CH4 - moles_CH40
    dmoles_H2 = moles_H2 - moles_H20
    dmoles_H2O = moles_H2O - moles_H2O0
    dmoles_O2 = moles_O2 - moles_O20
    dPres = Pres - Pres0
    
    return(dmoles_CO2,dmoles_CO,dmoles_CH4,dmoles_H2,dmoles_H2O,dmoles_O2,dPres,moles_O2_ash,moles_O2_lava,moles_O2_nt,moles_O2_red) #,dmoles_atm,dRmm_atm,dPres