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
import time
import random
import pandas as pd
from datetime import datetime
from modular_functions_clean import time_evol

now = datetime.now()
dt_string = now.strftime("%d-%m-%Y_%H%M")
nsteps = 100000 #number of timesteps

#Atmospheric species
Rmm_CO2 = 0.04401 #kg mol^-1
Rmm_H2O = 0.01801528 #kg mol^-1
Rmm_H = 0.001 #kg mol^-1 
Rmm_O2 = 0.032 #kg mol^-1 
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


# rand_start = random.randint(0,795) # Use if running a large number of runs to select random runs to start from to minimize chance of identical runs occuring simultaneously
# for i0 in range(rand_start,796):
    
for i0 in range(0,1): # For individual runs with manually set parameters

    print(i0)
    
    ## Load parameters for runs from 'gridsearch_all.csv'
    # param = np.genfromtxt('gridsearch_all.csv',delimiter=',',invalid_raise = False)
    # p_count = np.size(param)
    
    # while p_count < 7: # if p_count file being accessed by another iteration of script, this will tell code to wait 1 second then try to read again to prevent an error
    #     time.sleep(1)
    #     param = np.genfromtxt('gridsearch_all.csv',delimiter=',',invalid_raise = False)
    #     p_count = np.size(param)

    
    # param = np.genfromtxt('gridsearch_all.csv',delimiter=',',invalid_raise = False,usecols=(0,1,2,3,4,5))

    ## Extract 
    # tstart = param[i0,0]
    # mH2Otot = round(param[i0,1] * 1e-2,3)
    # gwater_GEL = param[i0,2]
    # mCO2tot = param[i0,3] * 1e-6
    # frac_volc = param[i0,4]
    # frac_ext = param[i0,5]
    
    tstart = 0.5
    mH2Otot = 0.001 * 1e-2
    gwater_GEL = 100
    mCO2tot = 500 * 1e-6
    frac_volc = 1.0
    frac_ext = 1.0
    
    print('tstart',tstart,'q',q,'epsilon',epsilon,'gwater_GEL',gwater_GEL,'CO2',mCO2tot,'H2O',mH2Otot,'FMQ',FMQ,'fext',frac_ext,
          'fvolc',frac_volc,'b_er=',b_err,'Tmelt',T,'T_atm',T_atm)
    #%%
  
    # Run tmospheric evolution model
    t,i1,Pres,moles_atm,Rmm_atm,moles_CO2,moles_H2O,moles_O2,dz,molesH2O_loss,T_surf,ash,lava,magma,nt,z_meltlayer_tot,z_max = time_evol(tstart,nsteps,metamorphic_degas,
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
    
    i1 = int(np.max(i1))
    
    if total_erupt > z_max: # Saves FAILED runs where erupted volume needed to produce all of Venus' modern CO2 exceeds mantle processing upper limit
        title2 = ("test_plotZFAIL_tstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + ".csv")
    else:
        if i1 >= int(nsteps-1000): # Saves runs as complete if they have reached end of nsteps
            title2 = ("test_plottstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + ".csv")
        else: # Saves partial runs to be picked up later
            title2 = ("test_plotpartial_tstart" + t_string + "_gwater" + gwater_string + "_H2O" + h2o_string + "_CO2" + CO2_string + "_f_ext" + f_ext_string + "_f_volc" + f_volc_string + ".csv")
    
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
            'moles_H2O':moles_H2O[0:i1],
            'lost_H2O': molesH2O_loss[0:i1],
            'moles_O2':moles_O2[0:i1],
            'Rmm_atm':Rmm_atm[0:i1],
            'T_surf':T_surf[0:i1],
            'ash':ash[0:i1],
            'lava':lava[0:i1],
            'magma':magma[0:i1],
            'nonthermal':nt[0:i1],
            'z_melt':z_meltlayer_tot[0:i1]}
        
    df = pd.DataFrame(data)
    # print(df)
    df.to_csv(title2,header = False) 
    
    del df
    del data