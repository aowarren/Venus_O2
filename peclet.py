# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 10:20:10 2022

@author: sasha
"""


import numpy as np
import math
from matplotlib import pyplot as plt
from VolcGases.functions import solve_gases

hyp = np.genfromtxt('ven_hyp.csv',delimiter=',') #Square root cumulative area percentage of Venus as a function of elevation (Fig 2 Rosenblatt et al. 1994 GRL)
frac = np.flip(1-((hyp[:,0]**2)/100.)) #Convert quare root cumulative area percentage to area fraction
topo_km = np.flip(hyp[:,1]) #Extract elevation (given as planetary radius) from hyposometry data
topo_m = (topo_km - np.min(topo_km))*1000. #Convert elevation to m and subtract lowest point on Venus from all others 
A_venus = 4.6e14 #m2

vol = topo_m * frac #Volume of planet beneath each elevation = topo_m*frac*A_Venus, dividing by A_Venus gives water inventory as GEL = topo_m*frac

#Set up grid to find H2O and CO2 solubilities as a function of pressure
n_run = 100
Pgrid = np.linspace(0.001,200,num=n_run) 
# For degassing calculations: VolcGasses input parameters for FMQ = 0
#A - C are constants used as inputs for VolcGasses
A = 25738.
B = 9.
C = 0.092
FMQ = 0. #Melt fO2
T = 1000. + 273. #melt temperature, K
Rmm_H2O = 0.018 #kg mol^-1


#%%

eclogite = np.genfromtxt('eclogite_transition.csv',delimiter=',') #P-T data for eclogite transition (Bucher & Frey 2002 "Metamorphic Rocks" in Petrogenesis of Metamorphic Rocks)
T_eclo = 273 + eclogite[:,0] #convert to K
GPa_eclo = eclogite[:,1]

serp = np.genfromtxt('serpent_dehyd2.csv',delimiter=',') #P-T data for serpentinite dehydration (Hacker et al. 2003 JGR: Solid Earth)
T_de = 273+ serp[:,0] #convert to K
GPa_de = serp[:,1] 
mCO2 = 300e-6 #melt CO2 concentration for all runs
t = 0.5 #end of habitable era in Gyrs after Venus formation (4.5 - t gives time in Ga)
secingyr = 60*60*24*365*1e9 #seconds in 1 Gyr
g = 8.87 #m s^-2
rho = 2800 #kg m^-3, Venus crust and melts
k = 4. #W m^-1 #K^-1
cp = 1.e3 #J kg^-1 K^-1
alpha = k/(rho*cp)
T_surf = 300. #K
T_m = 1300+273 #K, T at base of lithosphere defined by 1300C isotherm
deltaT = T_m - T_surf #K, T difference between surface & base of lithosphere


#Set up axes for future plots
fig1, ax1 = plt.subplots(figsize=(25,20))
figb, bx1 = plt.subplots(figsize=(25,20))

ax1.set_ylim((0,5))
ax1.set_xlim((T_surf,T_m))

ax1.set_ylabel('Pressure (GPa)',weight='bold')
ax1.set_xlabel('Temperature (K)',weight='bold')

bx1.set_ylim((0,2))
bx1.set_xlim((T_surf,1350))

bx1.set_ylabel('Pressure (GPa)',weight='bold', fontsize=50)
bx1.set_xlabel('Temperature (K)',weight='bold', fontsize=50)

bx1.tick_params(axis='both', which='major', labelsize=40)


#%%

D_arr = 1e3 * np.array([25,50,100,150,200]) #Array of lithosphere thicknesses in m
A_arr = 1e-6 * np.array([1,2,5]) # micro W m^-3, array of volumetric heating rates in crust
PCO2 = 80e5 #Total degassed CO2 during habitable era, Pa. Equivalent to (1-f_volc)*PCO2_pres

#Setup for plotting
cmap = plt.get_cmap('gnuplot_r') 
colors = cmap(np.linspace(0, 1, len(D_arr)+1))
styles = ['-.','--',':']


#Calculate analytical solution to heat transfer eq. for advection and diffusion ONLY (O'Reilly & Davies 1981 GRL)
for i1 in range(0,np.size(D_arr)):
    
    D = D_arr[i1] # m 
    
    z_crust = PCO2 / (g * mCO2 * rho) #thickness of crust produced needed to degas PCO2. Assumes all melt degasses (intrusions shallow)
    v = z_crust/(t*secingyr) #convert total crust production to advection velocity (i.e. downward velocity of buried lithosphere due to new crust produced at/near surface)
    lamb = (D*v*rho*cp)/k #cluster terms into single expression for convenience (following O'Reilly & Davies 1981 GRL)
    
    z = np.linspace(0,D,num=100,endpoint=True) #Set up array of points through lithosphere to calculate T-depth profile
    P = z*rho*g #Convert depth to pressure
    
    wiggle = z/D #cluster terms into single expression for convenience (following O'Reilly & Davies 1981 GRL)
    
    T_z = T_surf + (deltaT*(np.exp(lamb*wiggle) - 1.))/(np.exp(lamb) - 1.) #analytical solution for T-depth profile for advection & diffusion w/o heat production (O'Reilly & Davies 1981 GRL)
    
    bx1.plot(T_z,P/1e9,color=colors[i1+1],linewidth=7,label = str('z$_{lith}$ = ' + str(float('%g' % (D/1e3))) + ' km,\nH = 0 $\mu$W m$^{-3}$'))
    
    #%%
    
#Calculate analytical solution to heat transfer eq. for advection, diffusion, AND internal heat production (O'Reilly & Davies 1981 GRL, Kankanamge & Moore 2019 JGR:Planets)
for i2 in range(0,np.size(A_arr)):    
    for i1 in range(0,np.size(D_arr)):
        D = D_arr[i1] # m 
        A = A_arr[i2] # W m^-3
        
        z_crust = PCO2 / (g * mCO2 * rho)
        v = z_crust/(t*secingyr)
        
        z = np.linspace(0,D,num=100,endpoint=True)
        P = z*rho*g
        
        l = alpha/v #cluster terms into single expression for convenience

        #Calculate constants for analytical solution (Kankanamge & Moore 2019 JGR:Planets)
        c1 = (v/alpha)* ((deltaT - (A*D/(rho*cp*v))) / (np.exp(v*D/alpha ) - 1))
        c2 = T_surf - c1*alpha/v
        
        #Analytical solution
        T_z = (A*z/(rho*cp*v)) + c1*(alpha/v)*np.exp(v*z/alpha) + c2
        
        if np.max(T_z) <= T_m: #Exclude unphysical solutions where T in lithosphere exceeds 1300C (defines base of lithosphere)
            bx1.plot(T_z,P/1e9,color=colors[i1+1],linewidth=7,linestyle =styles[i2], label = str('z$_{lith}$ = ' + str(float('%g' % (D/1e3))) + ' km,\nH = ' + str(float('%g' % (A*1e6))) + ' $\mu$W m$^{-3}$'))

#%%
#Setup for plots
ax1.invert_yaxis()
bx1.invert_yaxis()
bx1.legend(bbox_to_anchor=(1, 1), fontsize=30)

bx1.plot(T_eclo,GPa_eclo,color = 'yellowgreen',linewidth=7)
bx1.plot(T_de,GPa_de,color = 'dodgerblue',linewidth=7)
bx1.text(990,0.3,'Serpentinite' '\n' 'dehydration',color='dodgerblue', fontsize=40)
bx1.text(940,1.52,'Eclogite' '\n' 'formation',color='yellowgreen',rotation='-15', fontsize=40)

bx1.set_title(str(float('%g' % (PCO2/1e5))) + ' bar CO$_{2}$ degassed by ' + str(float('%g' % (4.5 - t))) + ' Ga',weight='bold',fontsize=50)
    
#%%
GEL = 3000. #m, initial water inventory at start of habitable era
rho_w = 1000. #kg m^{-3}, density of water
m_w = GEL*rho_w #kg, mass of total water column
ven_hab_t = 0.1 #Gyr, start time of habitable era on Venus in Gyr after Venus formation (4.5 - t gives time in Ga)
    
hyd_arr = np.logspace(-3,math.log10(0.03),num=5,endpoint=True) #array for mass fraction water in hydrated crust (Schaefer & Sasselov 2015 ApJ)
f_mH2O = np.array([0.001,0.1,0.2,0.5,0.7,1]) #array for wt% water in melts on Venus

f_ext = 1 #extrusive volcanism fraction

#Setup for plots
cmap = cmap = plt.get_cmap('viridis_r')
colors = cmap(np.linspace(0, 1, len(f_mH2O)))
cmap = cmap = plt.get_cmap('plasma_r')
colors_hyd = cmap(np.linspace(0, 1, len(hyd_arr)))


#%%
fig3, ax3 = plt.subplots(2,1,figsize=(25,30))
fig3.tight_layout(pad=20.0)
ax3[0].set_ylim((10,5*GEL))
ax3[0].set_yscale("log")
ax3[0].axhline(y=500,linestyle='--',color='k',linewidth = 7)
ax3[0].axhline(y=1000,linestyle=':',color='k',linewidth = 7)
# ax3[0].axhline(y=500/1000,linestyle='--',color='dodgerblue',linewidth = 7)
ax3[0].set_xlabel('melt H$\mathregular{_{2}}$O concentration (wt%)',weight='bold',fontsize=40)
ax3[0].set_ylabel('GEL water (m)',weight='bold',fontsize=40)

ax3[1].set_ylim((100,5*GEL))
ax3[1].set_yscale("log")
ax3[1].axhline(y=500,linestyle='--',color='k',linewidth = 7)
ax3[1].axhline(y=1000,linestyle=':',color='k',linewidth = 7)
# ax3[1].axhline(y=500/1000,linestyle='--',color='dodgerblue',linewidth = 7)
ax3[1].set_xlabel('melt H$\mathregular{_{2}}$O concentration (wt%)',weight='bold',fontsize=40)
ax3[1].set_ylabel('GEL water (m)',weight='bold',fontsize=40)

# ax3[0].text(0,0.54,'500 m', fontsize=40,weight='bold',color='dodgerblue')
ax3[0].text(0,520,'500 m', fontsize=40,weight='bold')
ax3[0].text(0,1100,'1000 m', fontsize=40,weight='bold')

# ax3[0].annotate('', xy=(0.17,0.2),  xycoords='data',
# xytext=(0.17, 0.5), textcoords='data',
# arrowprops=dict(color='dodgerblue', shrink=0,width=10,headwidth=20,headlength=40),
# horizontalalignment='left', verticalalignment='bottom',
# )

ax3[0].annotate('', xy=(0.1,300),  xycoords='data',
xytext=(0.1, 500), textcoords='data',
arrowprops=dict(color='k', shrink=0,width=10,headwidth=20,headlength=40),
horizontalalignment='left', verticalalignment='bottom',
)

# ax3[1].text(0.9,0.54,'500 m', fontsize=40,weight='bold',color='dodgerblue')
ax3[1].text(0,520,'500 m', fontsize=40,weight='bold')
ax3[1].text(0,1100,'1000 m', fontsize=40,weight='bold')

# ax3[1].annotate('', xy=(1.02,0.2),  xycoords='data',
# xytext=(1.02, 0.5), textcoords='data',
# arrowprops=dict(color='dodgerblue', shrink=0,width=10,headwidth=20,headlength=40),
# horizontalalignment='left', verticalalignment='bottom',
# )

ax3[1].annotate('', xy=(0.9,300),  xycoords='data',
xytext=(0.9, 500), textcoords='data',
arrowprops=dict(color='k', shrink=0,width=10,headwidth=20,headlength=40),
horizontalalignment='left', verticalalignment='bottom',
)


for i4 in range(0,np.size(hyd_arr)):
    
    hyd = hyd_arr[i4] 
    rem = v*hyd*rho #kg s^{-1}, removal rate of hydrated crust by burial
    
    #%%
    fig2, ax2 = plt.subplots()
    ax2.set_ylim((100,GEL))
    ax2.set_yscale("log")
    ax2.axhline(y=500,linestyle='--',color='k')
    ax2.axhline(y=1000,linestyle=':',color='k')
    # ax2.axhline(y=500/1000,linestyle='--',color='dodgerblue')
    ax2.set_xlabel('Time (Gyr)',weight='bold')
    ax2.set_ylabel('GEL water (m)',weight='bold')
    
    # ax2.text(4.5-0.1,0.54,'500 m', weight='bold',color='dodgerblue')
    ax2.text(4.5-0.15,520,'500 m',weight='bold')
    ax2.text(4.5-0.15,1010,'1000 m',weight='bold')
    
    # ax2.annotate('', xy=(4.5-0.101,0.2),  xycoords='data',
    # xytext=(4.5-0.101, 0.5), textcoords='data',
    # arrowprops=dict(color='dodgerblue', shrink=0,width=2,headwidth=4,headlength=5),
    # horizontalalignment='left', verticalalignment='bottom',
    # )
    
    ax2.annotate('', xy=(4.5-0.15,300),  xycoords='data',
    xytext=(4.5-0.15, 500), textcoords='data',
    arrowprops=dict(color='k', shrink=0,width=2,headwidth=4,headlength=5),
    horizontalalignment='left', verticalalignment='bottom',
    )
    #%%
    
    f_rep_arr = f_mH2O/(100*hyd) #Ratio of melt H2O concentration to max. mass fraction water in hydrated crust
    
    time_arr = np.linspace(ven_hab_t,t,num=1000)*secingyr #Time array for water inventory evolution
    t_evol = np.zeros(np.size(time_arr)) #Array to record water inventory at each timestep
    t_evol[0] = m_w #Set initial water inventory 
    
    #Arrays to record water inventory evol. model end inventory and mean inventory
    scatter_end = np.zeros(np.size(f_rep_arr)) 
    scatter_mean = np.zeros(np.size(f_rep_arr))
    
    for i2 in range(0,np.size(f_rep_arr)):
        f_rep = f_rep_arr[i2]
        
        sH2O = np.zeros(n_run) #solubility of H2O in melt reaching surface, kg/mol
        for i1p in range(1,n_run):
            Pn = Pgrid[i1p]        
            log_FMQ = (-A/T+B+C*(Pn-1)/T)
            f_O2 = 10**(log_FMQ+FMQ)
            x = 0.01550152865954013 #defined as mol/g of magma
            
            #Calculate H2O solubility in melt as a function of pressure for each point in Pgrid
            P_H2O,P_H2,P_CO2,P_CO,P_CH4,alphaG,x_CO2,x_H2O = solve_gases(T,Pn,f_O2,mCO2,f_rep*hyd)
            
            #H2O solubility as mass fraction
            sH2O[i1p] = x*x_H2O*Rmm_H2O*1000
            
        for i5 in range(0,1):
            if i5 == 0:
                f_ext = 1
            if i5 == 1:
                f_ext = 0.1
            for i3 in range(1,np.size(time_arr)):
                Pres =  10 #bar, fixed value <10 bar. For evolving CO2 pressure can use: (time_arr[i3]/time_arr[-1])*(PCO2/1e5) 
                #interpolate Pgrid solubility calculations to find H2O solubility for current atmospheric pressure
                solH2O = np.interp(Pres,Pgrid,sH2O)
                
                #Use melt H2O concentration and H2O solubility to find degassed H2O:
                deH2O = (f_rep*hyd) - solH2O
                if deH2O < 0:
                    deH2O = 0

                f_de = deH2O/hyd #ratio of mass fraction water added by degassing to mass fraction water removed by hydrated crust

                #find fraction of surface area covered by water inventory using hyposmetry data:
                if t_evol[i3-1]/rho_w <= 0:
                    prop = 0.
                    t_evol[i3-1] = 0.
                elif t_evol[i3-1]/rho_w > 10e3:
                    prop = 1.
                else:
                    prop = np.interp((t_evol[i3-1]/rho_w),vol,frac) 
                
                #Calculate water inventory at end of timestep due to water added by volcanism and removed by hydration
                t_evol[i3] = t_evol[i3-1] + (time_arr[i3] - time_arr[i3-1])*f_ext*((rem*f_de) - (prop*rem))
                if t_evol[i3] <= 0.:
                    t_evol[i3] = 1e-5
                
        
            t_evol_m = t_evol / rho_w #convert water inventory to m
            
            scatter_end[i2] = t_evol_m[i3] #record final water inventory at end of habitable era
            scatter_mean[i2] = np.mean(t_evol_m) #record mean water inventory over course of habitable era

            if i5 == 0:
                ax2.plot(4.5 - (time_arr/secingyr),t_evol_m,label=str(f_mH2O[i2]),color=colors[i2])
            if i5 == 1:
                ax2.plot(time_arr/secingyr,t_evol_m,linestyle=':')
            
            ax2.legend(bbox_to_anchor=(1.1, 1),title = 'Melt H$\mathregular{_{2}}$O \nconcentration (wt%)')
            ax2.set_title('End of habitable era = ' + str(4.5-t) + ' Ga,\nHydrated mass fraction = ' + '{0:.1g}'.format(hyd*100) + ' wt%,\n' + str(mCO2*1e6) + ' ppm CO$\mathregular{_{2}}$',weight='bold')
    
    ax2.invert_xaxis()
    #%%%
    ax3[0].set_title('a. End habitable era water inventory\n('+ str(int(GEL)) + ' m GEL, ' + str(4.5-t) + ' Ga, ' + str(mCO2*1e6) + ' ppm CO$\mathregular{_{2}}$)',weight='bold',fontsize=50)
    
    ax3[0].plot(f_mH2O,scatter_end,'-o',color=colors_hyd[i4],linewidth = 7,markersize=25,label=str(' f$\mathregular{_{hyd}}$ = ' + str('{0:.1g}'.format(100*hyd_arr[i4])) + ' wt%'))
    ax3[0].legend(loc='lower right',fontsize=35)
    ax3[1].set_title('b. Mean habitable era water inventory',weight='bold',fontsize=50)
    
    ax3[1].plot(f_mH2O,scatter_mean,':s',color=colors_hyd[i4],linewidth = 7,markersize=25, label=str('{0:.1g}'.format(100*hyd_arr[i4])))
    
    ax3[0].tick_params(axis='both', which='major', labelsize=40)
    ax3[1].tick_params(axis='both', which='major', labelsize=40)