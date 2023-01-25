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
from matplotlib.ticker import FormatStrFormatter

#%%
#To reproduce main text figures using newly generated data, you MUST run gen_data_forplots.py followed by statistics.py BEFORE maintext_plotting.py
#or, download existing data from paper at LINK

# exec(open('statistics.py').read())

#Import .csv file with last line results (generated using gen_data_forplots.py) for all baseline model runs OR with runaway greenhouse surface melting with 40Ar results added:
    
# 
# results_pd = pd.read_csv('C:\\Users\\sasha\\Box\\Clean Code\\results\\data_plots_melt_nd.csv_Ar.csv')


#If baseline, set melt = 0, if RGH surface melting, melt = 1
melt = 0

if melt == 1:
    results_pd = pd.read_csv('C:\\Users\\sasha\\Box\\Clean Code\\results\\data_plots_melt_wd.csv_Ar.csv')
elif melt == 0:
    results_pd = pd.read_csv('C:\\Users\\sasha\\Box\\Clean Code\\results\\data_plots_nomelt_0.csv_Ar.csv')

results = results_pd.to_numpy()

# mantle_redox = -1.

# title_m = 'TESTstatistics_melt.npz'
# title_nm = 'TESTstatistics_nomelt.npz'
# title_m = 'nd_statistics_melt.npz'
# title_nm = 'statistics_nomelt.npz'
title_m = 'wd_statistics_melt_CO.npz'
title_nm = 'statistics_nomelt_CO.npz'

#%%
melt_npz = np.load(title_m)
nomelt_npz = np.load(title_nm)

g_totals_m = np.nan_to_num(melt_npz['g_m'])
w_totals_m = np.nan_to_num(melt_npz['h_m'])
c_totals_m = np.nan_to_num(melt_npz['c_m'])
v_totals_m = np.nan_to_num(melt_npz['v_m'])
t_totals_m = np.nan_to_num(melt_npz['t_m'])
e_totals_m = np.nan_to_num(melt_npz['e_m'])
f_totals_m = np.nan_to_num(melt_npz['f_m'])
totals_m = np.nan_to_num(melt_npz['all_m'])

g_totals_nm = np.nan_to_num(nomelt_npz['g_nm'])
w_totals_nm = np.nan_to_num(nomelt_npz['h_nm'])
c_totals_nm = np.nan_to_num(nomelt_npz['c_nm'])
v_totals_nm = np.nan_to_num(nomelt_npz['v_nm'])
t_totals_nm = np.nan_to_num(nomelt_npz['t_nm'])
e_totals_nm = np.nan_to_num(nomelt_npz['e_nm'])
f_totals_nm = np.nan_to_num(nomelt_npz['f_nm'])
totals_nm = np.nan_to_num(nomelt_npz['all_nm'])

g_totals_Ar_m = np.nan_to_num(melt_npz['g_m_Ar'])
w_totals_Ar_m = np.nan_to_num(melt_npz['h_m_Ar'])
c_totals_Ar_m = np.nan_to_num(melt_npz['c_m_Ar'])
v_totals_Ar_m = np.nan_to_num(melt_npz['v_m_Ar'])
t_totals_Ar_m = np.nan_to_num(melt_npz['t_m_Ar'])
e_totals_Ar_m = np.nan_to_num(melt_npz['e_m_Ar'])
f_totals_Ar_m = np.nan_to_num(melt_npz['f_m_Ar'])


g_totals_Ar_nm = np.nan_to_num(nomelt_npz['g_nm_Ar'])
w_totals_Ar_nm = np.nan_to_num(nomelt_npz['h_nm_Ar'])
c_totals_Ar_nm = np.nan_to_num(nomelt_npz['c_nm_Ar'])
v_totals_Ar_nm = np.nan_to_num(nomelt_npz['v_nm_Ar'])
t_totals_Ar_nm = np.nan_to_num(nomelt_npz['t_nm_Ar'])
e_totals_Ar_nm = np.nan_to_num(nomelt_npz['e_nm_Ar'])
f_totals_Ar_nm = np.nan_to_num(nomelt_npz['f_nm_Ar'])


g_conf_m = 100*np.nan_to_num(melt_npz['cg_m'])
w_conf_m = 100*np.nan_to_num(melt_npz['ch_m'])
c_conf_m = 100*np.nan_to_num(melt_npz['cc_m'])
v_conf_m = 100*np.nan_to_num(melt_npz['cv_m'])
t_conf_m = 100*np.nan_to_num(melt_npz['ct_m'])
e_conf_m = 100*np.nan_to_num(melt_npz['ce_m'])
f_conf_m = 100*np.nan_to_num(melt_npz['cf_m'])
conf_m = 100*np.nan_to_num(melt_npz['call_m'])

g_conf_nm = 100*np.nan_to_num(nomelt_npz['cg_nm'])
w_conf_nm = 100*np.nan_to_num(nomelt_npz['ch_nm'])
c_conf_nm = 100*np.nan_to_num(nomelt_npz['cc_nm'])
v_conf_nm = 100*np.nan_to_num(nomelt_npz['cv_nm'])
t_conf_nm = 100*np.nan_to_num(nomelt_npz['ct_nm'])
e_conf_nm = 100*np.nan_to_num(nomelt_npz['ce_nm'])
f_conf_nm = 100*np.nan_to_num(nomelt_npz['cf_nm'])
conf_nm = 100*np.nan_to_num(nomelt_npz['call_nm'])



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

final_tstep = results[:,0]

tstart = results[:,1]
t_arr = np.unique(tstart)

f_ext = results[:,3]
fext_arr = np.unique(f_ext)

f_volc = results[:,4]
fvolc_arr = np.unique(f_volc)

z_erupt = results[:,6]

gwater = results[:,8]
gwat_arr = np.unique(gwater)

concCO2 = results[:,9]
CO2_arr = np.unique(concCO2)
CO2_plotarr = np.array([300,500,1000,2000])

concH2O = results[:,10]
concH2O[concH2O>0.6e-2] = np.round(concH2O[concH2O>0.6e-2],3)
H2O_arr = np.unique(concH2O)

FMQ = results[:,11]
FMQ_arr = np.unique(FMQ)


t_arr = np.array([0.5,1.5,3])
fext_arr = np.array([0.1,0.3,0.5,1])
fvolc_arr = np.array([0.1,0.5,0.9,1])
gwat_arr = np.array([10,50,100,300,500,700,1000])
CO2_arr = np.array([0.0003,0.0005,0.001,0.002])
H2O_arr = 1e-2*np.array([0.001,0.1,0.2,0.5,0.7,1.])
H2O_arr[4] = 7e-3
FMQ_arr = np.array([0,-1,-2,-3,-4])

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


capN = moles_CO*(-2) + moles_CH4*(-8) + moles_H2*(-2) + moles_O2*(4)




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
#success conditions: Pres = 93 bar, O2 < measured, H2O < measured

end = np.size(tstart)

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
print(success)

#%%

fig0,ax0 = plt.subplots()
scatter = ax0.scatter(FMQ[success==1]+(0.5*np.random.rand(np.size(FMQ[success==1]))),capN[success==1],c=np.log10(gwater[success==1]),cmap='viridis_r',alpha=0.5)
# ax0.set_yscale('log')
ax0.invert_xaxis()
cbar = fig0.colorbar(scatter)

fig0a,ax0a = plt.subplots()
scatter = ax0a.scatter(FMQ+(0.5*np.random.rand(np.size(FMQ))),capN,c=success,cmap='viridis_r',alpha=0.1)
# ax0.set_yscale('log')
ax0a.invert_xaxis()

#%%
#Generate individual contour plots for combinations of 2 parameters for corner plot:

#GEL and tstart
g_t = np.zeros([np.size(gwat_arr),np.size(t_arr)])

for i1 in range(0,np.size(gwat_arr)):
    for i2 in range(0,np.size(t_arr)):
        ind1 = np.where(gwater==gwat_arr[i1])
        ind2 = np.where(tstart == t_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        # print(ind)
        g_t[i1,i2] = 100*np.mean(success[ind])
        print(g_t[i1,i2])
        
        print('gwater',gwat_arr[i1],'t',t_arr[i2],'ind',np.size(ind))

fig_m,axs = plt.subplots(7,7,figsize=(85,80))
if melt == 0:
    fig_m.suptitle("a. No runaway greenhouse melting",fontsize=150,weight='bold')
elif melt == 1:
    fig_m.suptitle("b. With runaway greenhouse melting",fontsize=150,weight='bold')
    
fig_m.tight_layout(pad=15)
axs[0,1].axis('off')
axs[0,2].axis('off')
axs[0,3].axis('off')
axs[0,4].axis('off')
axs[0,5].axis('off')
axs[0,6].axis('off')
axs[1,2].axis('off')
axs[1,3].axis('off')
axs[1,4].axis('off')
axs[1,5].axis('off')
axs[1,6].axis('off')
axs[2,3].axis('off')
axs[2,4].axis('off')
axs[2,5].axis('off')
axs[2,6].axis('off')
axs[3,4].axis('off')
axs[3,5].axis('off')
axs[3,6].axis('off')
axs[4,5].axis('off')
axs[4,6].axis('off')
axs[5,6].axis('off')


fig,ax = plt.subplots(figsize=(15,10))

img = ax.imshow(g_t,vmin=0, vmax=20,  aspect='auto')
ax.set_yticklabels(np.append(0,gwat_arr))
ax.set_xticklabels(np.append(0,4.5-t_arr))

cbar = fig.colorbar(img,orientation='horizontal')
cbar.set_label('% successful runs',size=30)

ax.invert_yaxis()


vector= np.nan_to_num(np.concatenate(g_t))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax.contour(g_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax.contour(g_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

cbar.ax.tick_params(labelsize=40) 


img_g_t = axs[1,0].imshow(g_t, vmin=0, vmax=20, aspect='auto')
axs[1,0].set_xticks(np.linspace(-1,2,4))
axs[1,0].set_yticks(np.linspace(-1,6,8))
axs[1,0].set_yticklabels(axs[1,0].get_yticks(), rotation = 45)
axs[1,0].set_yticklabels(np.append(0,gwat_arr))
axs[1,0].set_xticklabels([])
axs[1,0].invert_yaxis()


CS = axs[1,0].contour(g_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[1,0].contour(g_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[1,0].set_ylabel('Habitable era' '\n' ' water (m GEL)', fontsize=80,weight='bold')
axs[1,0].tick_params(axis='both', which='major', labelsize=60)



#%%
#Melt water conc and tstart
w_t = np.zeros([np.size(H2O_arr),np.size(t_arr)])

for i1 in range(0,np.size(H2O_arr)):
    for i2 in range(0,np.size(t_arr)):
        ind1 = np.where(concH2O==H2O_arr[i1])
        ind2 = np.where(tstart == t_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        # print(ind)
        w_t[i1,i2] = 100*np.mean(success[ind])
        
fig2,ax2 = plt.subplots(figsize=(15,10))

img = ax2.imshow(w_t, vmin=0, vmax=20,  aspect='auto')
ax2.set_yticklabels(np.round(np.append(0,1e2*H2O_arr),3))
ax2.set_xticklabels(np.append(0,4.5-t_arr))
cbar = fig2.colorbar(img)

vector= np.nan_to_num(np.concatenate(w_t))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax2.contour(w_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))

CS = ax2.contour(w_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


ax2.set_ylabel('Melt H$\mathregular{_{2}}$O (wt%)', fontsize=80)
ax2.invert_yaxis()
ax2.set_xlabel('End of habitable era (Ga)', fontsize=80)
ax2.tick_params(axis='both', which='major', labelsize=60)
ax2.set_title('% successful runs', fontsize=80)
ax2.set_ylim([0,70])
cbar.ax.tick_params(labelsize=40) 
# ax2.invert_yaxis()

img_w_t = axs[2,0].imshow(w_t,vmin=0, vmax=20,  aspect='auto')
axs[2,0].set_xticks(np.linspace(-1,2,4))
axs[2,0].set_yticks(np.linspace(-1,5,7))
axs[2,0].set_xticklabels(axs[2,0].get_xticks(), rotation = 45)
axs[2,0].set_yticklabels(axs[2,0].get_yticks(), rotation = 45)
axs[2,0].set_yticklabels(np.append(0,np.round(1e2*H2O_arr,3)))
axs[2,0].invert_yaxis()
axs[2,0].set_xticklabels([])


CS = axs[2,0].contour(w_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[2,0].contour(w_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[2,0].set_ylabel('Melt H$\mathregular{_{2}}$O (wt%)', fontsize=80,weight='bold')
axs[2,0].tick_params(axis='both', which='major', labelsize=60)


#%%
#FMQ and tstart
f_t = np.zeros([np.size(FMQ_arr),np.size(t_arr)])

for i1 in range(0,np.size(FMQ_arr)):
    for i2 in range(0,np.size(t_arr)):
        ind1 = np.where(FMQ==FMQ_arr[i1])
        ind2 = np.where(tstart == t_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        # print(ind)
        f_t[i1,i2] = 100*np.mean(success[ind])
        
fig2a,ax2a = plt.subplots(figsize=(15,10))

img = ax2a.imshow(f_t, vmin=0, vmax=20,  aspect='auto')
ax2a.set_yticklabels(np.append(0,FMQ_arr))
ax2a.set_xticklabels(np.append(0,4.5-t_arr))
cbar = fig2a.colorbar(img)

vector= np.nan_to_num(np.concatenate(f_t))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax2.contour(f_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))

CS = ax2.contour(f_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


ax2a.set_ylabel('$\mathregular{\Delta}$FMQ', fontsize=80)
ax2a.invert_yaxis()
ax2a.set_xlabel('End of habitable era (Ga)', fontsize=80)
ax2a.tick_params(axis='both', which='major', labelsize=60)
ax2a.set_title('% successful runs', fontsize=80)
cbar.ax.tick_params(labelsize=60) 
# ax2.invert_yaxis()

img_f_t = axs[6,0].imshow(f_t,vmin=0, vmax=20,  aspect='auto')
axs[6,0].set_xticks(np.linspace(-1,2,4))
axs[6,0].set_yticks(np.linspace(-1,4,6)) #FMQ
axs[6,0].set_xticklabels(axs[6,0].get_xticks(), rotation = 45)
axs[6,0].set_yticklabels(axs[6,0].get_yticks(), rotation = 45)
axs[6,0].set_yticklabels(np.append(0,FMQ_arr))
axs[6,0].set_xticklabels(np.append(0,4.5-t_arr))
axs[6,0].invert_yaxis()


CS = axs[6,0].contour(f_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[6,0].contour(f_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[6,0].set_xlabel('End of habitable\nera (Ga)', fontsize=80,weight='bold')
axs[6,0].set_ylabel('$\mathregular{\Delta}$FMQ', fontsize=80,weight='bold')
axs[6,0].tick_params(axis='both', which='major', labelsize=60)


#%% f_ext f_volc

v_e = np.zeros([np.size(fvolc_arr),np.size(fext_arr)])

for i1 in range(0,np.size(fvolc_arr)):
    for i2 in range(0,np.size(fext_arr)):
        ind1 = np.where(f_volc==fvolc_arr[i1])
        ind2 = np.where(f_ext == fext_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        # print(ind)
        v_e[i1,i2] = 100*np.mean(success[ind])

fig3,ax3 = plt.subplots(figsize=(15,10))

img = ax3.imshow(v_e, vmin=0, vmax=20,  aspect='auto')
ax3.set_xticks(np.linspace(-1,3,5))
ax3.set_yticks(np.linspace(-1,3,5))
ax3.set_yticklabels(np.append(0,fvolc_arr))
ax3.set_xticklabels(np.append(0,fext_arr))
cbar = fig3.colorbar(img)

vector= np.nan_to_num(np.concatenate(v_e))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax3.contour(v_e,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax3.contour(v_e,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax3.set_ylabel('f$\mathregular{_{volc}}$', fontsize=80)
ax3.set_xlabel('f$\mathregular{_{ext}}$', fontsize=80)
ax3.tick_params(axis='both', which='major', labelsize=60)
ax3.set_title('% successful runs', fontsize=80)
ax3.invert_yaxis()


img_v_e = axs[5,4].imshow(v_e, vmin=0, vmax=20,  aspect='auto')
axs[5,4].set_xticks(np.linspace(-1,3,5))
axs[5,4].set_yticks(np.linspace(-1,3,5))
# axs[5,4].set_xticklabels(axs[5,4].get_xticks(), rotation = 45)
axs[5,4].set_yticklabels([])
axs[5,4].set_xticklabels([])
axs[5,4].invert_yaxis()


CS = axs[5,4].contour(v_e,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[5,4].contour(v_e,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

# axs[5,4].set_xlabel('f$\mathregular{_{ext}}$', fontsize=100,weight='bold')
axs[5,4].tick_params(axis='both', which='major', labelsize=60)

#%% FMQ f_ext 

f_e = np.zeros([np.size(FMQ_arr),np.size(fext_arr)])

for i1 in range(0,np.size(FMQ_arr)):
    for i2 in range(0,np.size(fext_arr)):
        ind1 = np.where(FMQ==FMQ_arr[i1])
        ind2 = np.where(f_ext == fext_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        # print(ind)
        f_e[i1,i2] = 100*np.mean(success[ind])

fig3a,ax3a = plt.subplots(figsize=(15,10))
cbar = fig3a.colorbar(img)

img = ax3a.imshow(f_e, vmin=0, vmax=20,  aspect='auto')
ax3a.set_xticks(np.linspace(-1,3,5))
ax3a.set_yticks(np.linspace(-1,4,6))
ax3a.set_yticklabels(np.append(0,FMQ_arr))
ax3a.set_xticklabels(np.append(0,fext_arr))

vector= np.nan_to_num(np.concatenate(f_e))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax3a.contour(f_e,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax3a.contour(f_e,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax3a.set_ylabel('$\mathregular{\Delta}$FMQ', fontsize=80)
ax3a.set_xlabel('f$\mathregular{_{ext}}$', fontsize=80)
ax3a.tick_params(axis='both', which='major', labelsize=60)
ax3a.set_title('% successful runs', fontsize=80)
ax3a.invert_yaxis()


img_f_e = axs[6,4].imshow(f_e, vmin=0, vmax=20,  aspect='auto')
axs[6,4].set_xticks(np.linspace(-1,3,5))
axs[6,4].set_yticks(np.linspace(-1,4,6))
axs[6,4].set_xticklabels(axs[6,4].get_xticks(), rotation = 45)
axs[6,4].set_yticklabels([])
axs[6,4].set_xticklabels(np.append(0,fext_arr))
axs[6,4].invert_yaxis()


CS = axs[6,4].contour(f_e,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[6,4].contour(f_e,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[6,4].set_xlabel('f$\mathregular{_{ext}}$', fontsize=100,weight='bold')
axs[6,4].tick_params(axis='both', which='major', labelsize=60)

#%% FMQ f_volc

f_v = np.zeros([np.size(FMQ_arr),np.size(fvolc_arr)])

for i1 in range(0,np.size(FMQ_arr)):
    for i2 in range(0,np.size(fvolc_arr)):
        ind1 = np.where(FMQ==FMQ_arr[i1])
        ind2 = np.where(f_volc == fvolc_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        # print(ind)
        f_v[i1,i2] = 100*np.mean(success[ind])

fig3b,ax3b = plt.subplots(figsize=(15,10))

img = ax3b.imshow(f_v, vmin=0, vmax=20,  aspect='auto')
ax3b.set_xticks(np.linspace(-1,3,5))
ax3b.set_yticks(np.linspace(-1,4,6))
ax3b.set_yticklabels(np.append(0,FMQ_arr))
ax3b.set_xticklabels(np.append(0,fvolc_arr))
cbar = fig3b.colorbar(img)

vector= np.nan_to_num(np.concatenate(f_v))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax3b.contour(f_v,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax3b.contour(f_v,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax3b.set_ylabel('$\mathregular{\Delta}$FMQ', fontsize=80)
ax3b.set_xlabel('f$\mathregular{_{volc}}$', fontsize=80)
ax3b.tick_params(axis='both', which='major', labelsize=60)
ax3b.set_title('% successful runs', fontsize=80)
ax3b.invert_yaxis()


img_f_v = axs[6,5].imshow(f_v, vmin=0, vmax=20,  aspect='auto')
axs[6,5].set_xticks(np.linspace(-1,3,5))
axs[6,5].set_yticks(np.linspace(-1,4,6))
axs[6,5].set_xticklabels(axs[6,5].get_xticks(), rotation = 45)
axs[6,5].set_yticklabels([])
axs[6,5].set_xticklabels(np.append(0,fvolc_arr))
axs[6,5].invert_yaxis()


CS = axs[6,5].contour(f_v,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[6,5].contour(f_v,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[6,5].set_xlabel('f$\mathregular{_{volc}}$', fontsize=100,weight='bold')
axs[6,5].tick_params(axis='both', which='major', labelsize=60)


#%%
#Melt water conc and GEL
w_g = np.zeros([np.size(H2O_arr),np.size(gwat_arr)])

for i1 in range(0,np.size(H2O_arr)):
    for i2 in range(0,np.size(gwat_arr)):
        ind1 = np.where(concH2O==H2O_arr[i1])
        ind2 = np.where(gwater == gwat_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        w_g[i1,i2] = 100*np.mean(success[ind])

fig4,ax4 = plt.subplots(figsize=(15,10))

img = ax4.imshow(w_g, vmin=0, vmax=20, aspect='auto')
ax4.set_xticks(np.linspace(-1,6,8))
ax4.set_yticks(np.linspace(-1,5,7))
ax4.set_yticklabels(np.round(np.append(0,1e2*H2O_arr),3))
ax4.set_xticklabels(np.append(0,gwat_arr))

cbar = fig4.colorbar(img)

vector= np.nan_to_num(np.concatenate(w_g))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1

CS = ax4.contour(w_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax4.contour(w_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


ax4.set_ylabel('Melt H$\mathregular{_{2}}$O (wt%)', fontsize=80)
ax4.set_xlabel('Initial water (m GEL)', fontsize=80)
ax4.tick_params(axis='both', which='major', labelsize=60)
ax4.set_title('% successful runs', fontsize=80)
ax4.invert_yaxis()


img_w_g = axs[2,1].imshow(w_g, vmin=0, vmax=20,  aspect='auto')
axs[2,1].set_xticks(np.linspace(-1,6,8))
axs[2,1].set_yticks(np.linspace(-1,5,7))
axs[2,1].set_xticklabels([])
axs[2,1].set_yticklabels([])
axs[2,1].invert_yaxis()

CS = axs[2,1].contour(w_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[2,1].contour(w_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[2,1].tick_params(axis='both', which='major', labelsize=60)


#%%
#Melt water conc and melt co2 conc
c_w = np.zeros([np.size(CO2_arr),np.size(H2O_arr)])

for i1 in range(0,np.size(CO2_arr)):
    for i2 in range(0,np.size(H2O_arr)):
        ind1 = np.where(concCO2 == CO2_arr[i1])
        ind2 = np.where(concH2O==H2O_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        c_w[i1,i2] = 100*np.mean(success[ind])

fig5,ax5 = plt.subplots(figsize=(15,10))

img = ax5.imshow(c_w, vmin=0, vmax=20,  aspect='auto')
ax5.set_yticks(np.linspace(-1,3,5))
ax5.set_yticklabels(np.append(0,CO2_plotarr))
ax5.set_xticklabels(np.round(np.append(0,1e2*H2O_arr),3))
cbar = fig5.colorbar(img)

vector= np.nan_to_num(np.concatenate(c_w))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1

CS = ax5.contour(c_w,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax5.contour(c_w,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax5.set_ylabel('Melt CO$\mathregular{_{2}}$ (ppm)', fontsize=80)
ax5.set_xlabel('Melt H$\mathregular{_{2}}$O (wt%)', fontsize=80)
ax5.tick_params(axis='both', which='major', labelsize=60)
ax5.set_title('% successful runs', fontsize=80)
ax5.invert_yaxis()

img_c_w = axs[3,2].imshow(c_w, vmin=0, vmax=20,  aspect='auto')
axs[3,2].set_yticks(np.linspace(-1,3,5))
axs[3,2].set_xticks(np.linspace(-1,5,7))
axs[3,2].set_xticklabels([])
axs[3,2].set_yticklabels([])
axs[3,2].invert_yaxis()

CS = axs[3,2].contour(c_w,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[3,2].contour(c_w,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[3,2].tick_params(axis='both', which='major', labelsize=60)

#%%
#Melt CO2 conc and tstart
c_t = np.zeros([np.size(CO2_arr),np.size(t_arr)])

for i1 in range(0,np.size(CO2_arr)):
    for i2 in range(0,np.size(t_arr)):
        ind2 = np.where(tstart==t_arr[i2])
        ind1 = np.where(concCO2 == CO2_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        # print(ind)
        c_t[i1,i2] = 100*np.mean(success[ind])

fig6,ax6 = plt.subplots(figsize=(15,10))

img = ax6.imshow(c_t, vmin=0, vmax=20,  aspect='auto')
ax6.set_yticks(np.linspace(-1,3,5))
ax6.set_xticks(np.linspace(-1,2,4))
ax6.set_xticklabels(np.append(0,4.5-t_arr))
ax6.set_yticklabels(np.append(0,CO2_plotarr))
cbar = fig6.colorbar(img)

vector= np.nan_to_num(np.concatenate(c_t))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
if np.size(np.where(pdf<0.05))<1:
    sig2 = -1
else:
    ii_2sig = np.max(np.where(pdf<0.05))
    sig2 = sort[ii_2sig]
sig1 = sort[ii_1sig]
if sig2 >= sig1:
    sig2 = -1

CS = ax6.contour(c_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax6.contour(c_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


ax6.set_xlabel('End of habitable era (Ga)', fontsize=80)
ax6.set_ylabel('Melt CO$\mathregular{_{2}}$ (ppm)', fontsize=80)
ax6.tick_params(axis='both', which='major', labelsize=60)
ax6.set_title('% successful runs', fontsize=80)
ax6.invert_yaxis()

img_c_t = axs[3,0].imshow(c_t, vmin=0, vmax=20,  aspect='auto')
axs[3,0].set_yticks(np.linspace(-1,3,5))
axs[3,0].set_xticks(np.linspace(-1,2,4))
axs[3,0].set_xticklabels(axs[3,0].get_xticks(), rotation = 45)
axs[3,0].set_yticklabels(axs[3,0].get_yticks(), rotation = 45)
axs[3,0].set_xticklabels([])
axs[3,0].set_yticklabels(np.append(0,CO2_plotarr))
axs[3,0].invert_yaxis()

CS = axs[3,0].contour(c_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[3,0].contour(c_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[3,0].set_ylabel('Melt CO$\mathregular{_{2}}$ (ppm)', fontsize=80,weight='bold')
axs[3,0].tick_params(axis='both', which='major', labelsize=60)


#%%
#Melt CO2 conc and GEL
c_g = np.zeros([np.size(CO2_arr),np.size(gwat_arr)])

for i1 in range(0,np.size(CO2_arr)):
    for i2 in range(0,np.size(gwat_arr)):
        ind1 = np.where(concCO2 == CO2_arr[i1])
        ind2 = np.where(gwater==gwat_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        c_g[i1,i2] = 100*np.mean(success[ind])

fig7,ax7 = plt.subplots(figsize=(15,10))

img = ax7.imshow(c_g, vmin=0, vmax=20,  aspect='auto')
ax7.set_yticks(np.linspace(-1,3,5))
ax7.set_yticklabels(np.append(0,CO2_plotarr))
ax7.set_xticklabels(np.append(0,gwat_arr))
cbar = fig7.colorbar(img)

vector= np.nan_to_num(np.concatenate(c_g))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax7.contour(c_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax7.contour(c_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax7.set_ylabel('Melt CO$\mathregular{_{2}}$ (ppm)', fontsize=80)
ax7.set_xlabel('Initial water (m GEL)', fontsize=80)
ax7.tick_params(axis='both', which='major', labelsize=60)
ax7.set_title('% successful runs', fontsize=80)
ax7.invert_yaxis()


img_c_g = axs[3,1].imshow(c_g, vmin=0, vmax=20,  aspect='auto')
axs[3,1].set_yticks(np.linspace(-1,3,5))
axs[3,1].set_xticks(np.linspace(-1,6,8))
axs[3,1].set_xticklabels([])
axs[3,1].set_yticklabels([])
axs[3,1].invert_yaxis()

CS = axs[3,1].contour(c_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[3,1].contour(c_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[3,1].tick_params(axis='both', which='major', labelsize=60)


#%%
#FMQ and GEL
f_g = np.zeros([np.size(FMQ_arr),np.size(gwat_arr)])

for i1 in range(0,np.size(FMQ_arr)):
    for i2 in range(0,np.size(gwat_arr)):
        ind1 = np.where(FMQ == FMQ_arr[i1])
        ind2 = np.where(gwater==gwat_arr[i2])
        ind = np.intersect1d(ind1,ind2)
        f_g[i1,i2] = 100*np.mean(success[ind])

fig7a,ax7a = plt.subplots(figsize=(15,10))

img = ax7.imshow(f_g, vmin=0, vmax=20,  aspect='auto')
ax7a.set_yticks(np.linspace(-1,3,5))
ax7a.set_yticklabels(np.append(0,CO2_plotarr))
ax7a.set_xticklabels(np.append(0,gwat_arr))
cbar = fig7a.colorbar(img)

vector= np.nan_to_num(np.concatenate(f_g))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax7.contour(f_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax7.contour(f_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax7a.set_ylabel('$\mathregular{\Delta}$FMQ', fontsize=80)
ax7a.set_xlabel('Initial water (m GEL)', fontsize=80)
ax7a.tick_params(axis='both', which='major', labelsize=60)
ax7a.set_title('% successful runs', fontsize=80)
ax7a.invert_yaxis()


img_f_g = axs[6,1].imshow(f_g, vmin=0, vmax=20,  aspect='auto')
axs[6,1].set_yticks(np.linspace(-1,4,6))
axs[6,1].set_xticks(np.linspace(-1,6,8))
axs[6,1].set_xticklabels([])
axs[6,1].set_yticklabels([])
axs[6,1].set_xticklabels(axs[6,1].get_xticks(), rotation = 45)
axs[6,1].set_xticklabels(np.append(0,gwat_arr))
axs[6,1].invert_yaxis()

CS = axs[6,1].contour(f_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[6,1].contour(f_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[6,1].tick_params(axis='both', which='major', labelsize=60)
axs[6,1].set_xlabel('Habitable era' '\n' ' water (m GEL)', fontsize=80,weight='bold')

#%%
#f-ext and tstart
e_t = np.zeros([np.size(fext_arr),np.size(t_arr)])

for i1 in range(0,np.size(fext_arr)):
    for i2 in range(0,np.size(t_arr)):
        ind2 = np.where(tstart==t_arr[i2])
        ind1 = np.where(f_ext == fext_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        e_t[i1,i2] = 100*np.mean(success[ind])

fig8,ax8 = plt.subplots(figsize=(15,10))

img = ax8.imshow(e_t, vmin=0, vmax=20,  aspect='auto')
ax8.set_yticks(np.linspace(-1,3,5))
ax8.set_xticks(np.linspace(-1,2,4))
ax8.set_xticklabels(np.append(0,4.5-t_arr))
ax8.set_yticklabels(np.append(0,fext_arr))
cbar = fig8.colorbar(img)

vector= np.nan_to_num(np.concatenate(e_t))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax8.contour(e_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax8.contour(e_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax8.set_xlabel('End of habitable era (Ga)', fontsize=80)
ax8.set_ylabel('f$\mathregular{_{ext}}$', fontsize=90)
ax8.tick_params(axis='both', which='major', labelsize=60)
ax8.set_title('% successful runs', fontsize=80)
ax8.invert_yaxis()

img_e_t = axs[4,0].imshow(e_t, vmin=0, vmax=20,  aspect='auto')
axs[4,0].set_yticks(np.linspace(-1,3,5))
axs[4,0].set_xticks(np.linspace(-1,2,4))
axs[4,0].set_xticklabels(axs[4,0].get_xticks(), rotation = 45)
axs[4,0].set_yticklabels(axs[4,0].get_yticks(), rotation = 45)
axs[4,0].set_xticklabels([])
axs[4,0].set_yticklabels(np.append(0,fext_arr))
axs[4,0].invert_yaxis()


CS = axs[4,0].contour(e_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[4,0].contour(e_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


axs[4,0].set_ylabel('f$\mathregular{_{ext}}$', fontsize=100,weight='bold')
axs[4,0].tick_params(axis='both', which='major', labelsize=60)


#%%
#f_volc and tstart
v_t = np.zeros([np.size(fvolc_arr),np.size(t_arr)])

for i1 in range(0,np.size(fvolc_arr)):
    for i2 in range(0,np.size(t_arr)):
        ind2 = np.where(tstart==t_arr[i2])
        ind1 = np.where(f_volc == fvolc_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        v_t[i1,i2] = 100*np.mean(success[ind])


fig9,ax9 = plt.subplots(figsize=(15,10))

img = ax9.imshow(v_t, vmin=0, vmax=20,  aspect='auto')
ax9.set_xticks(np.linspace(-1,2,4))
ax9.set_yticks(np.linspace(-1,3,5))
ax9.set_xticklabels(np.append(0,4.5-t_arr))
ax9.set_yticklabels(np.append(0,fvolc_arr))
cbar = fig9.colorbar(img)


vector= np.nan_to_num(np.concatenate(v_t))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax9.contour(v_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax9.contour(v_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax9.set_xlabel('End of habitable era (Ga)', fontsize=80)
ax9.set_ylabel('f$\mathregular{_{volc}}$', fontsize=80)
ax9.tick_params(axis='both', which='major', labelsize=60)
ax9.set_title('% successful runs', fontsize=80)
ax9.invert_yaxis()

img_v_t = axs[5,0].imshow(v_t, vmin=0, vmax=20,  aspect='auto')
axs[5,0].set_xticks(np.linspace(-1,2,4))
axs[5,0].set_yticks(np.linspace(-1,3,5))
axs[5,0].set_xticklabels(axs[5,0].get_xticks(), rotation = 45)
axs[5,0].set_yticklabels(axs[5,0].get_yticks(), rotation = 45)
# axs[5,0].set_yticklabels([])
axs[5,0].set_xticklabels([])
axs[5,0].set_yticklabels(np.append(0,fvolc_arr))
axs[5,0].invert_yaxis()

CS = axs[5,0].contour(v_t,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[5,0].contour(v_t,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


axs[5,0].set_ylabel('f$\mathregular{_{volc}}$', fontsize=100,weight='bold')
axs[5,0].tick_params(axis='both', which='major', labelsize=60)

#%%
#f_volc and GEL
v_g = np.zeros([np.size(fvolc_arr),np.size(gwat_arr)])

for i1 in range(0,np.size(fvolc_arr)):
    for i2 in range(0,np.size(gwat_arr)):
        ind2 = np.where(gwater==gwat_arr[i2])
        ind1 = np.where(f_volc == fvolc_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        v_g[i1,i2] = 100*np.mean(success[ind])


fig10,ax10 = plt.subplots(figsize=(15,10))

img = ax10.imshow(v_g, vmin=0, vmax=20,  aspect='auto')
ax10.set_xticks(np.linspace(-1,6,8))
ax10.set_yticks(np.linspace(-1,3,5))
ax10.set_xticklabels(np.append(0,gwat_arr))
ax10.set_yticklabels(np.append(0,fvolc_arr))
cbar = fig10.colorbar(img)

vector= np.nan_to_num(np.concatenate(v_g))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]

if sig2 >= sig1:
    sig2 = -1
    
CS = ax10.contour(v_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax10.contour(v_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax10.set_xlabel('Initial water (m GEL)', fontsize=80)
ax10.set_ylabel('f$\mathregular{_{volc}}$', fontsize=80)
ax10.tick_params(axis='both', which='major', labelsize=60)
ax10.set_title('% successful runs', fontsize=80)
ax10.invert_yaxis()

img_v_g = axs[5,1].imshow(v_g, vmin=0, vmax=20,  aspect='auto')
axs[5,1].set_xticks(np.linspace(-1,6,8))
axs[5,1].set_yticks(np.linspace(-1,3,5))
# axs[5,1].set_xticklabels(axs[5,1].get_xticks(), rotation = 45)
axs[5,1].set_yticklabels(axs[5,1].get_yticks(), rotation = 45)
axs[5,1].set_yticklabels([])
axs[5,1].set_xticklabels([])
axs[5,1].invert_yaxis()

CS = axs[5,1].contour(v_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[5,1].contour(v_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


axs[5,1].tick_params(axis='both', which='major', labelsize=60)


#%%
#f_volc and melt water conc
v_w = np.zeros([np.size(fvolc_arr),np.size(H2O_arr)])

for i1 in range(0,np.size(fvolc_arr)):
    for i2 in range(0,np.size(H2O_arr)):
        ind2 = np.where(concH2O==H2O_arr[i2])
        ind1 = np.where(f_volc == fvolc_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        v_w[i1,i2] = 100*np.mean(success[ind])

fig11,ax11 = plt.subplots(figsize=(15,20))

img = ax11.imshow(v_w, vmin=0, vmax=20,  aspect='auto')
ax11.set_xticks(np.linspace(-1,5,7))
ax11.set_yticks(np.linspace(-1,3,5))
ax11.set_xticklabels(np.round(np.append(0,H2O_arr*1e2),3))
ax11.set_yticklabels(np.append(0,fvolc_arr))
cbar = fig11.colorbar(img)

vector= np.nan_to_num(np.concatenate(v_w))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax11.contour(v_w,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax11.contour(v_w,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax11.set_xlabel('Melt H$\mathregular{_{2}}$O (wt%)', fontsize=80)
ax11.set_ylabel('f$\mathregular{_{volc}}$', fontsize=80)
ax11.tick_params(axis='both', which='major', labelsize=60)
ax11.set_title('% successful runs', fontsize=80)
ax11.invert_yaxis()


img_v_w = axs[5,2].imshow(v_w, vmin=0, vmax=20,  aspect='auto')
axs[5,2].set_yticks(np.linspace(-1,3,5))
axs[5,2].set_xticks(np.linspace(-1,5,7))
# axs[5,2].set_xticklabels(axs[5,2].get_xticks(), rotation = 45)
axs[5,2].set_yticklabels([])
axs[5,2].set_xticklabels([])
axs[5,2].invert_yaxis()

CS = axs[5,2].contour(v_w,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[5,2].contour(v_w,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


axs[5,2].tick_params(axis='both', which='major', labelsize=60)


#%%
#f_volc and melt CO2 conc
v_c = np.zeros([np.size(fvolc_arr),np.size(CO2_arr)])

for i1 in range(0,np.size(fvolc_arr)):
    for i2 in range(0,np.size(CO2_arr)):
        ind2 = np.where(concCO2==CO2_arr[i2])
        ind1 = np.where(f_volc == fvolc_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        v_c[i1,i2] = 100*np.mean(success[ind])

fig12,ax12 = plt.subplots(figsize=(15,30))

img = ax12.imshow(v_c, vmin=0, vmax=20,  aspect='auto')
ax12.set_xticks(np.linspace(-1,3,5))
ax12.set_yticks(np.linspace(-1,3,5))
ax12.set_xticklabels(np.append(0,CO2_plotarr))
ax12.set_yticklabels(np.append(0,fvolc_arr))
cbar = fig12.colorbar(img)

vector= np.nan_to_num(np.concatenate(v_c))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
if np.size(np.where(pdf<0.05))<1:
    sig2 = -1
else:
    ii_2sig = np.max(np.where(pdf<0.05))
    sig2 = sort[ii_2sig]
sig1 = sort[ii_1sig]

if sig2 >= sig1:
    sig2 = -1
CS = ax12.contour(v_c,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax12.contour(v_c,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax12.set_xlabel('Melt CO$\mathregular{_{2}}$ (ppm)', fontsize=80)
ax12.set_ylabel('f$\mathregular{_{volc}}$', fontsize=80)
ax12.tick_params(axis='both', which='major', labelsize=60)
ax12.set_title('% successful runs', fontsize=80)
ax12.invert_yaxis()

img_v_c = axs[5,3].imshow(v_c, vmin=0, vmax=20,  aspect='auto')
axs[5,3].set_xticks(np.linspace(-1,3,5))
axs[5,3].set_yticks(np.linspace(-1,3,5))
# axs[5,3].set_xticklabels(axs[5,3].get_xticks(), rotation = 45)
axs[5,3].set_yticklabels([])
axs[5,3].set_xticklabels([])
# axs[5,3].set_xticklabels(np.append(0,CO2_plotarr))
axs[5,3].invert_yaxis()

CS = axs[5,3].contour(v_c,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[5,3].contour(v_c,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))


axs[5,3].tick_params(axis='both', which='major', labelsize=60)


#%% f_ext and GEL

e_g = np.zeros([np.size(fext_arr),np.size(gwat_arr)])

for i1 in range(0,np.size(fext_arr)):
    for i2 in range(0,np.size(gwat_arr)):
        ind2 = np.where(gwater==gwat_arr[i2])
        ind1 = np.where(f_ext == fext_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        e_g[i1,i2] = 100*np.mean(success[ind])

fig13,ax13 = plt.subplots(figsize=(15,10))

img = ax13.imshow(e_g, vmin=0, vmax=20,  aspect='auto')
ax13.set_yticks(np.linspace(-1,3,5))
ax13.set_xticks(np.linspace(-1,6,8))
ax13.set_xticklabels(np.append(0,gwat_arr))
ax13.set_yticklabels(np.append(0,fext_arr))
cbar = fig13.colorbar(img)

vector= np.nan_to_num(np.concatenate(e_g))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax13.contour(e_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax13.contour(e_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax13.set_xlabel('Initial water (m GEL)', fontsize=80)
ax13.set_ylabel('f$\mathregular{_{ext}}$', fontsize=80)
ax13.tick_params(axis='both', which='major', labelsize=60)
ax13.set_title('% successful runs', fontsize=80)
ax13.invert_yaxis()

img_e_g = axs[4,1].imshow(e_g, vmin=0, vmax=20,  aspect='auto')
axs[4,1].set_yticks(np.linspace(-1,3,5))
axs[4,1].set_xticks(np.linspace(-1,6,8))
axs[4,1].set_xticklabels([])
axs[4,1].set_yticklabels([])
axs[4,1].invert_yaxis()

CS = axs[4,1].contour(e_g,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[4,1].contour(e_g,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[4,1].tick_params(axis='both', which='major', labelsize=60)

#%% f_ext mH2O

e_w = np.zeros([np.size(fext_arr),np.size(H2O_arr)])

for i1 in range(0,np.size(fext_arr)):
    for i2 in range(0,np.size(H2O_arr)):
        ind2 = np.where(concH2O==H2O_arr[i2])
        ind1 = np.where(f_ext == fext_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        e_w[i1,i2] = 100*np.mean(success[ind])
        # if melt ==1:
        #     e_w[e_w>10]=0.8

fig14,ax14 = plt.subplots(figsize=(15,10))

img = ax14.imshow(e_w, vmin=0, vmax=20,  aspect='auto')
ax14.set_yticks(np.linspace(-1,3,5))
ax14.set_xticks(np.linspace(-1,5,7))
ax14.set_xticklabels(np.round(np.append(0,1e2*H2O_arr),3))
ax14.set_yticklabels(np.append(0,fext_arr))
cbar = fig14.colorbar(img)

vector= np.nan_to_num(np.concatenate(e_w))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax14.contour(e_w,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax14.contour(e_w,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax14.set_xlabel('Melt H$\mathregular{_{2}}$O (wt%)', fontsize=80)
ax14.set_ylabel('f$\mathregular{_{ext}}$', fontsize=80)
ax14.tick_params(axis='both', which='major', labelsize=60)
ax14.set_title('% successful runs', fontsize=80)
ax14.invert_yaxis()

img_e_w = axs[4,2].imshow(e_w, vmin=0, vmax=20,  aspect='auto')
axs[4,2].set_yticks(np.linspace(-1,3,5))
axs[4,2].set_xticks(np.linspace(-1,5,7))

axs[4,2].set_xticklabels([])
axs[4,2].set_yticklabels([])
axs[4,2].invert_yaxis()

CS = axs[4,2].contour(e_w,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[4,2].contour(e_w,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[4,2].tick_params(axis='both', which='major', labelsize=60)

#%% FMQ mH2O

f_w = np.zeros([np.size(FMQ_arr),np.size(H2O_arr)])

for i1 in range(0,np.size(FMQ_arr)):
    for i2 in range(0,np.size(H2O_arr)):
        ind2 = np.where(concH2O==H2O_arr[i2])
        ind1 = np.where(FMQ == FMQ_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        f_w[i1,i2] = 100*np.mean(success[ind])
        # if melt ==1:
        #     f_w[f_w>10]=0.8

fig14a,ax14a = plt.subplots(figsize=(15,10))

img = ax14a.imshow(f_w, vmin=0, vmax=20,  aspect='auto')
ax14a.set_yticks(np.linspace(-1,4,6))
ax14a.set_xticks(np.linspace(-1,5,7))
ax14a.set_xticklabels(np.round(np.append(0,1e2*H2O_arr),3))
ax14a.set_yticklabels(np.append(0,FMQ_arr))
cbar = fig14a.colorbar(img)

vector= np.nan_to_num(np.concatenate(f_w))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
CS = ax14a.contour(f_w,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax14a.contour(f_w,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax14a.set_xlabel('Melt H$\mathregular{_{2}}$O (wt%)', fontsize=80)
ax14a.set_ylabel('$\mathregular{\Delta}$FMQ', fontsize=80)
ax14a.tick_params(axis='both', which='major', labelsize=60)
ax14a.set_title('% successful runs', fontsize=80)
ax14a.invert_yaxis()

img_f_w = axs[6,2].imshow(f_w, vmin=0, vmax=20,  aspect='auto')
axs[6,2].set_yticks(np.linspace(-1,4,6))
axs[6,2].set_xticks(np.linspace(-1,5,7))

axs[6,2].set_xticklabels([])
axs[6,2].set_yticklabels([])
axs[6,2].set_xticklabels(axs[6,2].get_xticks(), rotation = 45)
axs[6,2].set_xticklabels(np.append(0,np.round(1e2*H2O_arr,3)))
axs[6,2].invert_yaxis()

CS = axs[6,2].contour(f_w,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[6,2].contour(f_w,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[6,2].tick_params(axis='both', which='major', labelsize=60)
axs[6,2].set_xlabel('Melt H$\mathregular{_{2}}$O\n(wt%)', fontsize=80,weight='bold')

#%% f_ext melt CO2 conc

e_c = np.zeros([np.size(fext_arr),np.size(CO2_arr)])

for i1 in range(0,np.size(fext_arr)):
    for i2 in range(0,np.size(CO2_arr)):
        ind2 = np.where(concCO2==CO2_arr[i2])
        ind1 = np.where(f_ext == fext_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        e_c[i1,i2] = 100*np.mean(success[ind])

fig15,ax15 = plt.subplots(figsize=(15,10))

img = ax15.imshow(e_c, vmin=0, vmax=20,  aspect='auto')
ax15.set_xticks(np.linspace(-1,3,5))
ax15.set_yticks(np.linspace(-1,3,5))
ax15.set_xticklabels(np.append(0,CO2_plotarr))
ax15.set_yticklabels(np.append(0,fext_arr))
cbar = fig15.colorbar(img)

vector= np.nan_to_num(np.concatenate(e_c))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
if np.size(np.where(pdf<0.05))<1:
    sig2 = -1
else:
    ii_2sig = np.max(np.where(pdf<0.05))
    sig2 = sort[ii_2sig]
sig1 = sort[ii_1sig]
if sig2 >= sig1:
    sig2 = -1
    
CS = ax15.contour(e_c,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax15.contour(e_c,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax15.set_xlabel('Melt CO$\mathregular{_{2}}$ (ppm)', fontsize=80)
ax15.set_ylabel('f$\mathregular{_{ext}}$', fontsize=80)
ax15.tick_params(axis='both', which='major', labelsize=60)
ax15.set_title('% successful runs', fontsize=80)
ax15.invert_yaxis()

img_e_c = axs[4,3].imshow(e_c, vmin=0, vmax=20,  aspect='auto')
axs[4,3].set_xticks(np.linspace(-1,3,5))
axs[4,3].set_yticks(np.linspace(-1,3,5))
axs[4,3].set_xticklabels([])
axs[4,3].set_yticklabels([])
axs[4,3].invert_yaxis()

CS = axs[4,3].contour(e_c,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[4,3].contour(e_c,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[4,3].tick_params(axis='both', which='major', labelsize=60)
# axs[4,3].set_title('% successful runs', fontsize=80)
# c_label = fig_m.colorbar(img_e_c,ax=axs[4,3])

#%% FMQ CO2 conc

f_c = np.zeros([np.size(FMQ_arr),np.size(CO2_arr)])

for i1 in range(0,np.size(FMQ_arr)):
    for i2 in range(0,np.size(CO2_arr)):
        ind2 = np.where(concCO2==CO2_arr[i2])
        ind1 = np.where(FMQ == FMQ_arr[i1])
        ind = np.intersect1d(ind1,ind2)
        f_c[i1,i2] = 100*np.mean(success[ind])

fig15a,ax15a = plt.subplots(figsize=(15,10))

img = ax15a.imshow(f_c, vmin=0, vmax=20,  aspect='auto')
ax15a.set_xticks(np.linspace(-1,3,5))
ax15a.set_yticks(np.linspace(-1,4,6))
ax15a.set_xticklabels(np.append(0,CO2_plotarr))
ax15a.set_yticklabels(np.append(0,FMQ_arr))
cbar = fig15a.colorbar(img)

vector= np.nan_to_num(np.concatenate(f_c))
sort = np.sort(vector)
sum_sort = np.cumsum(sort)
pdf = sum_sort/ sum_sort[np.size(sum_sort)-1]
ii_1sig = np.max(np.where(pdf<0.32))
ii_2sig = np.max(np.where(pdf<0.05))
sig1 = sort[ii_1sig]
sig2 = sort[ii_2sig]
if sig2 >= sig1:
    sig2 = -1
    
CS = ax15a.contour(f_c,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = ax15a.contour(f_c,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

ax15a.set_xlabel('Melt CO$\mathregular{_{2}}$ (ppm)', fontsize=80)
ax15a.set_ylabel('f$\mathregular{_{ext}}$', fontsize=80)
ax15a.tick_params(axis='both', which='major', labelsize=60)
ax15a.set_title('% successful runs', fontsize=80)
ax15a.invert_yaxis()

img_f_c = axs[6,3].imshow(f_c, vmin=0, vmax=20,  aspect='auto')
axs[6,3].set_xticks(np.linspace(-1,3,5))
axs[6,3].set_yticks(np.linspace(-1,3,5))
axs[6,3].set_xticklabels([])
axs[6,3].set_yticklabels([])
axs[6,3].set_xticklabels(axs[5,3].get_xticks(), rotation = 45)
axs[6,3].set_xticklabels(np.append(0,CO2_plotarr))
axs[6,3].invert_yaxis()

CS = axs[6,3].contour(f_c,levels = [sig2],colors=('w',),linestyles=(':',),linewidths=(8,))
CS = axs[6,3].contour(f_c,levels = [sig1],colors=('w',),linestyles=('-',),linewidths=(8,))

axs[6,3].tick_params(axis='both', which='major', labelsize=60)
axs[6,3].set_xlabel('Melt CO$\mathregular{_{2}}$\n(ppm)', fontsize=80,weight='bold')
# axs[6,3].set_title('% successful runs', fontsize=80)
# c_label = fig_m.colorbar(img_f_c,ax=axs[4,3])


#%% Histograms for corner plot

if melt == 1:
    g_totals = g_totals_m
    w_totals = w_totals_m
    c_totals = c_totals_m
    v_totals = v_totals_m
    t_totals = t_totals_m
    e_totals = e_totals_m
    f_totals = f_totals_m

if melt == 0:
    g_totals = g_totals_nm
    w_totals = w_totals_nm
    c_totals = c_totals_nm
    v_totals = v_totals_nm
    t_totals = t_totals_nm
    e_totals = e_totals_nm
    f_totals = f_totals_nm
    

axs[0,0].bar(np.arange(np.size(t_totals)),t_totals,color='indigo')
# axs[0,0].set_ylim([0,3])
axs[0,0].set_xticklabels(axs[0,0].get_xticks(), rotation = 45)
# axs[0,0].set_yticklabels(axs[0,0].get_yticks(), rotation = 45)
axs[0,0].set_xticks(np.linspace(-1,2,4))
axs[0,0].set_xticklabels(np.append(" ",4.5-t_arr))
axs[0,0].tick_params(axis='both', which='major', labelsize=60)
axs[0,0].set_ylabel('% successful\nruns', fontsize=80,weight='bold')
# axs[0,0].set_xticklabels([])
axs[0,0].yaxis.set_major_formatter(FormatStrFormatter('%d%%'))

axs[1,1].bar(np.arange(np.size(g_totals)),g_totals,color='indigo')
# axs[1,1].set_ylim([0,5])
axs[1,1].set_xticklabels(axs[1,1].get_xticks(), rotation = 45)
# axs[1,1].set_yticklabels(axs[1,1].get_yticks(), rotation = 45)
axs[1,1].set_xticks(np.linspace(-1,6,8))
axs[1,1].set_xticklabels(np.append(" ",gwat_arr))
# axs[1,1].set_xticklabels([])
axs[1,1].tick_params(axis='both', which='major', labelsize=60)
axs[1,1].yaxis.set_major_formatter(FormatStrFormatter('%d%%'))


axs[2,2].bar(np.arange(np.size(w_totals)),w_totals,color='indigo')
# axs[2,2].set_ylim([0,5])
axs[2,2].set_xticklabels(axs[2,2].get_xticks(), rotation = 45)
# axs[2,2].set_yticklabels(axs[2,2].get_yticks(), rotation = 45)
axs[2,2].set_xticks(np.linspace(-1,5,7))
axs[2,2].set_xticklabels(np.append(" ",np.round(1e2*H2O_arr,3)))
# axs[2,2].set_xticklabels([])
axs[2,2].tick_params(axis='both', which='major', labelsize=60)
axs[2,2].yaxis.set_major_formatter(FormatStrFormatter('%d%%'))


axs[3,3].bar(np.arange(np.size(c_totals)),c_totals,color='indigo')
# axs[3,3].set_ylim([0,4])
axs[3,3].set_xticklabels(axs[3,3].get_xticks(), rotation = 45)
# axs[3,3].set_yticklabels(axs[3,3].get_yticks(), rotation = 45)
axs[3,3].set_xticks(np.linspace(-1,3,5))
axs[3,3].set_xticklabels(np.append(" ",CO2_plotarr))
# axs[3,3].set_xticklabels([])
axs[3,3].tick_params(axis='both', which='major', labelsize=60)
axs[3,3].yaxis.set_major_formatter(FormatStrFormatter('%d%%'))


axs[4,4].bar(np.arange(np.size(e_totals)),e_totals,color='indigo')
# axs[4,4].set_ylim([0,3])
axs[4,4].set_xticklabels(axs[4,4].get_xticks(), rotation = 45)
# axs[4,4].set_yticklabels(axs[4,4].get_yticks(), rotation = 45)
axs[4,4].set_xticks(np.linspace(-1,3,5))
axs[4,4].set_xticklabels(np.append(" ",fext_arr))
axs[4,4].tick_params(axis='both', which='major', labelsize=60)
axs[4,4].yaxis.set_major_formatter(FormatStrFormatter('%d%%'))
# axs[4,4].set_xticklabels([])


axs[5,5].bar(np.arange(np.size(v_totals)),v_totals,color='indigo')
# axs[5,5].set_ylim([0,4])
axs[5,5].set_xticklabels(axs[5,5].get_xticks(), rotation = 45)
# axs[5,5].set_yticklabels(axs[5,5].get_yticks(), rotation = 45)
axs[5,5].set_xticks(np.linspace(-1,3,5))
axs[5,5].set_xticklabels(np.append(" ",fvolc_arr))
axs[5,5].tick_params(axis='both', which='major', labelsize=60)
axs[5,5].yaxis.set_major_formatter(FormatStrFormatter('%d%%'))
# axs[5,5].set_xlabel('f$\mathregular{_{volc}}$', fontsize=100,weight='bold')
# axs[5,5].set_xticklabels([])

axs[6,6].bar(np.arange(np.size(f_totals)),f_totals,color='indigo')
# axs[6,6].set_ylim([0,4])
axs[6,6].set_xticklabels(axs[6,6].get_xticks(), rotation = 45)
# axs[6,6].set_yticklabels(axs[6,6].get_yticks(), rotation = 45)
axs[6,6].set_xticks(np.linspace(-1,4,6))
axs[6,6].set_xticklabels(np.append(" ",FMQ_arr))
axs[6,6].tick_params(axis='both', which='major', labelsize=60)
axs[6,6].set_xlabel('$\mathregular{\Delta}$FMQ', fontsize=100,weight='bold')
axs[6,6].yaxis.set_major_formatter(FormatStrFormatter('%d%%'))


#%%
fig,barax1 = plt.subplots(figsize=(15,10))

barax1.bar(np.arange(np.size(t_totals)),t_totals,color='indigo') #color=' indigo'
barax1.set_xticks(np.linspace(-1,2,4))
barax1.set_xticklabels(np.append(" ",(4.5-t_arr)))
barax1.tick_params(axis='both', which='major', labelsize=40)
barax1.set_xlabel('End of habitable era' '\n' '(Ga)', fontsize=40)
barax1.set_ylabel('% successful runs', fontsize=40)


fig,barax2 = plt.subplots(figsize=(15,10))

barax2.bar(np.arange(np.size(g_totals)),g_totals,color='indigo')
barax2.set_xticks(np.linspace(-1,6,8))
barax2.set_xticklabels(np.append(" ",gwat_arr))
barax2.tick_params(axis='both', which='major', labelsize=40)
barax2.set_xlabel('Habitable era' '\n' ' water (m GEL)', fontsize=40)
barax2.set_ylabel('% successful runs', fontsize=40)

#%%

fig,barax3 = plt.subplots(figsize=(15,10))

barax3.bar(np.arange(np.size(e_totals)),e_totals,color='indigo') #color=' indigo'
barax3.set_xticks(np.linspace(-1,3,5))
barax3.set_xticklabels(np.append(" ",(fext_arr)))
barax3.tick_params(axis='both', which='major', labelsize=40)
barax3.set_xlabel('f$_{ext}$', fontsize=40)
barax3.set_ylabel('% successful runs', fontsize=40)

#%%

#%%
fig_m,bxs = plt.subplots(4,2,figsize=(30,50))
fig_m.tight_layout(pad=15.0)

#t
n=np.size(t_totals)
r = np.arange(n)

width = 0.4
bxs[0,0].bar(r, t_totals_nm,  color = 'midnightblue',width = width, edgecolor = 'black')
bxs[0,0].bar(r + width, t_totals_m,  color = 'teal',width = width, edgecolor = 'black')
bxs[0,0].set_xticks(r + width/2)
# bxs[0,0].set_ylim([0,2.2])
bxs[0,0].set_xticklabels((4.5-t_arr))
bxs[0,0].tick_params(axis='both', which='major', labelsize=40)
bxs[0,0].set_ylabel('% successful runs', fontsize=50,weight='bold')
bxs[0,0].set_xlabel('End of habitable era (Ga)', fontsize=50,weight='bold')


n=np.size(g_totals)
r = np.arange(n)
width = 0.4
bxs[0,1].bar(r, g_totals_nm,  color = 'midnightblue',width = width, edgecolor = 'black',label='No RGH melting')
bxs[0,1].bar(r + width, g_totals_m,  color = 'teal',width = width, edgecolor = 'black',label='RGH melting')
bxs[0,1].set_xticks(r + width/2)
# bxs[0,1].set_ylim([0,5])
bxs[0,1].set_xticklabels((gwat_arr))
bxs[0,1].tick_params(axis='both', which='major', labelsize=40)
# bxs[0,1].set_ylabel('% successful runs', fontsize=50,weight='bold')
bxs[0,1].set_xlabel('Habitable Era Water (m GEL)', fontsize=50,weight='bold')
bxs[0,1].legend(fontsize=40)

#w
n=np.size(w_totals)
r = np.arange(n)
width = 0.4
bxs[1,0].bar(r, w_totals_nm,  color = 'midnightblue',width = width, edgecolor = 'black',label='No RGH melting')
bxs[1,0].bar(r + width, w_totals_m,  color = 'teal',width = width, edgecolor = 'black',label='RGH melting')
bxs[1,0].set_xticks(r + width/2)
bxs[1,0].set_xticklabels(np.round((1e2*H2O_arr),3))
bxs[1,0].tick_params(axis='both', which='major', labelsize=40)
bxs[1,0].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,0].legend(fontsize=30)
# bxs[1,0].set_ylim([0,2.7])
bxs[1,0].set_xlabel('Melt H$\mathregular{_{2}}$O\n(wt%)', fontsize=50,weight='bold')

#c
n=np.size(c_totals)
r = np.arange(n)
width = 0.4
bxs[1,1].bar(r, c_totals_nm,   color = 'midnightblue',width = width, edgecolor = 'black',label='No RGH melting')
bxs[1,1].bar(r + width, c_totals_m,  color = 'teal',width = width, edgecolor = 'black',label='RGH melting')
bxs[1,1].set_xticks(r + width/2)
bxs[1,1].set_xticklabels((CO2_plotarr))
bxs[1,1].tick_params(axis='both', which='major', labelsize=40)
# bxs[1,1].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,1].legend(fontsize=30)
# bxs[1,1].set_ylim([0,2.5])
bxs[1,1].set_xlabel('Melt CO$\mathregular{_{2}}$\n(ppm)', fontsize=50,weight='bold')


n=np.size(e_totals)

df_e = results_pd.groupby(['f_ext'])[['success']].agg('mean')

r = np.arange(n)
width = 0.4
bxs[2,0].bar(r, e_totals_nm,   color = 'midnightblue',width = width, edgecolor = 'black')
bxs[2,0].bar(r + width, e_totals_m,   color = 'teal',width = width, edgecolor = 'black')
bxs[2,0].set_xticks(r + width/2)
bxs[2,0].set_xticklabels((fext_arr))
bxs[2,0].tick_params(axis='both', which='major', labelsize=40)
bxs[2,0].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,1].legend(fontsize=30)
bxs[2,0].set_xlabel('f$\mathregular{_{ext}}$', fontsize=70,weight='bold')
bxs[2,0].scatter(r,e_totals)
# bxs[2,0].set_ylim([0,3.0])


#v
n=np.size(v_totals)
r = np.arange(n)
width = 0.4
bxs[2,1].bar(r, v_totals_nm,  color = 'midnightblue',width = width, edgecolor = 'black',label='No RGH melting')
bxs[2,1].bar(r + width, v_totals_m,   color = 'teal',width = width, edgecolor = 'black',label='RGH melting')
bxs[2,1].set_xticks(r + width/2)
bxs[2,1].set_xticklabels((fvolc_arr))
bxs[2,1].tick_params(axis='both', which='major', labelsize=40)
# bxs[2,1].set_ylim([0,2.5])
# bxs[1,2].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,2].legend(fontsize=30)
bxs[2,1].set_xlabel('f$\mathregular{_{volc}}$', fontsize=70,weight='bold')


#fmq
n=np.size(f_totals)

df_f = results_pd.groupby(['f_ext'])[['success']].agg('mean')

r = np.arange(n)
width = 0.4
bxs[3,0].bar(r, f_totals_nm,   color = 'midnightblue',width = width, edgecolor = 'black')
bxs[3,0].bar(r + width, f_totals_m,   color = 'teal',width = width, edgecolor = 'black')
bxs[3,0].set_xticks(r + width/2)
bxs[3,0].set_xticklabels((FMQ_arr))
bxs[3,0].tick_params(axis='both', which='major', labelsize=40)
bxs[3,0].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,1].legend(fontsize=30)
bxs[3,0].set_xlabel('f$\mathregular{_{ext}}$', fontsize=70,weight='bold')
bxs[3,0].scatter(r,f_totals)

# % Histograms with 40Ar results included
fig_m,cxs = plt.subplots(4,2,figsize=(30,50))
fig_m.tight_layout(pad=15.0)

#t
n=np.size(t_totals)
r = np.arange(n)

# print('mean',mean_succ,'95% confidence intervals',conf)
df_t = results_pd.groupby(['tstart'])[['success']].agg('mean')
t_totals_Ar = df_t.to_numpy()
t_totals_Ar = t_totals_Ar[:,0]

t_totals_Ar_m = t_totals_Ar * t_totals_m
t_totals_Ar_nm = t_totals_Ar * t_totals_nm

width = 0.4
bxs[0,0].bar(r, t_totals_Ar_nm,  color = 'cornflowerblue',width = width, edgecolor = 'black')
bxs[0,0].bar(r + width, t_totals_Ar_m,  color = 'lawngreen',width = width, edgecolor = 'black')
bxs[0,0].set_xticks(r + width/2)
bxs[0,0].set_ylim([0,6])
bxs[0,0].set_xticklabels((4.5-t_arr))
bxs[0,0].tick_params(axis='both', which='major', labelsize=40)
bxs[0,0].set_ylabel('% successful runs', fontsize=50,weight='bold')
bxs[0,0].set_xlabel('End of habitable era (Ga)', fontsize=50,weight='bold')
# bxs[0,0].errorbar(capsize = 10)


df_g = results_pd.groupby(['gwater'])[['success']].agg('mean')
g_totals_Ar = df_g.to_numpy()
g_totals_Ar = g_totals_Ar[:,0]

g_totals_Ar_m = g_totals_Ar * g_totals_m
g_totals_Ar_nm = g_totals_Ar * g_totals_nm

n=np.size(g_totals)
r = np.arange(n)
width = 0.4
bxs[0,1].bar(r, g_totals_Ar_nm,  color = 'cornflowerblue',width = width, edgecolor = 'black',label='No melting, $\mathregular{^{40}}$Ar')
bxs[0,1].bar(r + width, g_totals_Ar_m,  color = 'lawngreen',width = width, edgecolor = 'black',label='Melting, $\mathregular{^{40}}$Ar')
bxs[0,1].set_xticks(r + width/2)
bxs[0,1].set_ylim([0,6])
bxs[0,1].set_xticklabels((gwat_arr))
bxs[0,1].tick_params(axis='both', which='major', labelsize=40)
# bxs[0,1].set_ylabel('% successful runs', fontsize=50,weight='bold')
bxs[0,1].set_xlabel('Habitable Era Water (m GEL)', fontsize=50,weight='bold')
bxs[0,1].legend(fontsize=40)


df_w = results_pd.groupby(['H2O'])[['success']].agg('mean')
w_totals_Ar = df_w.to_numpy()
w_totals_Ar = w_totals_Ar[:,0]

w_totals_Ar_m = w_totals_Ar * w_totals_m
w_totals_Ar_nm = w_totals_Ar * w_totals_nm
#w
n=np.size(w_totals)
r = np.arange(n)
width = 0.4
bxs[1,0].bar(r, w_totals_Ar_nm,  color = 'cornflowerblue',width = width, edgecolor = 'black',label='No RGH melting')
bxs[1,0].bar(r + width, w_totals_Ar_m,  color = 'lawngreen',width = width, edgecolor = 'black',label='RGH melting')
bxs[1,0].set_xticks(r + width/2)
bxs[1,0].set_xticklabels(np.round(1e2*H2O_arr,3))
bxs[1,0].tick_params(axis='both', which='major', labelsize=40)
bxs[1,0].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,0].legend(fontsize=30)
bxs[1,0].set_ylim([0,6.5])
bxs[1,0].set_xlabel('Melt H$\mathregular{_{2}}$O (wt%)', fontsize=50,weight='bold')

#c

df_c = results_pd.groupby(['CO2'])[['success']].agg('mean')
c_totals_Ar = df_c.to_numpy()
c_totals_Ar = c_totals_Ar[:,0]

c_totals_Ar_m = c_totals_Ar * c_totals_m
c_totals_Ar_nm = c_totals_Ar * c_totals_nm

n=np.size(c_totals)
r = np.arange(n)
width = 0.4
bxs[1,1].bar(r, c_totals_Ar_nm,   color = 'cornflowerblue',width = width, edgecolor = 'black',label='No RGH melting')
bxs[1,1].bar(r + width, c_totals_Ar_m,  color = 'lawngreen',width = width, edgecolor = 'black',label='RGH melting')
bxs[1,1].set_xticks(r + width/2)
bxs[1,1].set_xticklabels((CO2_plotarr))
bxs[1,1].tick_params(axis='both', which='major', labelsize=40)
# bxs[1,1].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,1].legend(fontsize=30)
bxs[1,1].set_ylim([0,3])
bxs[1,1].set_xlabel('Melt CO$\mathregular{_{2}}$ (ppm)', fontsize=50,weight='bold')


#e
n=np.size(e_totals)

df_e = results_pd.groupby(['f_ext'])[['success']].agg('mean')
e_totals_Ar = df_e.to_numpy()
e_totals_Ar = e_totals_Ar[:,0]

e_totals_Ar_m = e_totals_Ar * e_totals_m
e_totals_Ar_nm = e_totals_Ar * e_totals_nm

r = np.arange(n)
width = 0.4
bxs[2,0].bar(r, e_totals_Ar_nm,   color = 'cornflowerblue',width = width, edgecolor = 'black')
bxs[2,0].bar(r + width, e_totals_Ar_m,   color = 'lawngreen',width = width, edgecolor = 'black')
bxs[2,0].set_xticks(r + width/2)
bxs[2,0].set_xticklabels((fext_arr))
bxs[2,0].tick_params(axis='both', which='major', labelsize=40)
bxs[2,0].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,1].legend(fontsize=30)
bxs[2,0].set_xlabel('f$\mathregular{_{ext}}$', fontsize=70,weight='bold')
# bxs[2,0].scatter(r,e_totals_Ar)
bxs[2,0].set_ylim([0,3])


df_v = results_pd.groupby(['f_volc'])[['success']].agg('mean')
v_totals_Ar = df_v.to_numpy()
v_totals_Ar = v_totals_Ar[:,0]

v_totals_Ar_m = v_totals_Ar * v_totals_m
v_totals_Ar_nm = v_totals_Ar * v_totals_nm

#v
n=np.size(v_totals)
r = np.arange(n)
width = 0.4
bxs[2,1].bar(r, v_totals_Ar_nm,  color = 'cornflowerblue',width = width, edgecolor = 'black',label='No RGH melting')
bxs[2,1].bar(r + width, v_totals_Ar_m,   color = 'lawngreen',width = width, edgecolor = 'black',label='RGH melting')
bxs[2,1].set_xticks(r + width/2)
bxs[2,1].set_xticklabels((fvolc_arr))
bxs[2,1].tick_params(axis='both', which='major', labelsize=40)
bxs[2,1].set_ylim([0,6])
# bxs[1,2].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,2].legend(fontsize=30)
bxs[2,1].set_xlabel('f$\mathregular{_{volc}}$', fontsize=70,weight='bold')


#FMQ
n=np.size(f_totals)

df_f = results_pd.groupby(['FMQ'])[['success']].agg('mean')
f_totals_Ar = df_f.to_numpy()
f_totals_Ar = f_totals_Ar[:,0]

f_totals_Ar_m = f_totals_Ar * f_totals_m
f_totals_Ar_nm = f_totals_Ar * f_totals_nm

r = np.arange(n)
width = 0.4
bxs[3,0].bar(r, f_totals_Ar_nm,   color = 'cornflowerblue',width = width, edgecolor = 'black')
bxs[3,0].bar(r + width, f_totals_Ar_m,   color = 'lawngreen',width = width, edgecolor = 'black')
bxs[3,0].set_xticks(r + width/2)
bxs[3,0].set_xticklabels((FMQ_arr))
bxs[3,0].tick_params(axis='both', which='major', labelsize=40)
bxs[3,0].set_ylabel('% successful runs', fontsize=50,weight='bold')
# bxs[1,1].legend(fontsize=30)
bxs[3,0].set_xlabel('$\mathregular{\Delta}$FMQ', fontsize=70,weight='bold')
# bxs[3,0].scatter(r,f_totals_Ar)
bxs[3,0].set_ylim([0,6.5])


bxs[0,0].text(-0.3,5.5,'a.', fontsize=50,weight='bold',color='black')
bxs[0,1].text(-0.3,5.5,'b.', fontsize=50,weight='bold',color='black')
bxs[1,0].text(-0.3,6,'c.', fontsize=50,weight='bold',color='black')
bxs[1,1].text(-0.3,2.75,'d.', fontsize=50,weight='bold',color='black')
bxs[2,0].text(-0.3,2.75,'e.', fontsize=50,weight='bold',color='black')
bxs[2,1].text(-0.3,5.5,'f.', fontsize=50,weight='bold',color='black')
bxs[3,0].text(-0.3,6,'g.', fontsize=50,weight='bold',color='black')