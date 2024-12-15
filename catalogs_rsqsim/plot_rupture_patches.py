#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 15:50:10 2024

@author: jackwilliams
"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker

import os
from osgeo import ogr

#load catalog rupture patch and event info 
all_rupture_patches=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_patches.npy')#event patch
all_rupture_events=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_events.npy')#event index
all_rupture_slip=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_slip.npy')#event patch slip


#read input fault info
#otago fault names only. Check is consistent with NSHM IFM list 
otago_fault_list=pd.read_csv('otago_rsqsim_catalog/otago_fault_list_20240428.csv', usecols=[0])
#geometry info about faults
otago_faults = pd.read_csv('otago_rsqsim_catalog/otago_faults_2500_tapered_slip.flt', sep=" ", header=None)

#load catalog
otago_rsqsim_catalog=pd.read_csv('otago_rsqsim_catalog/otago_1e6yr_nospinup_catalogue.csv')#overall catalog


#IMPORT AND CLEAN NZCFM AND COASTLINE DATA
os.chdir("..") #Go up one level in directory

nzcfm_ds=ogr.Open("gis_files/NZ_CFM_v1_0.shp", 1)
nzcfm_layer=nzcfm_ds.GetLayer()

nz_coast_ds=ogr.Open("gis_files/nz-coastline-polygon_nztm.shp", 0)
nz_coast_layer=nz_coast_ds.GetLayer()
paths = []

#Only read in South Island NZ Coast Feature
nz_coast_feature=nz_coast_layer.GetFeature(0)
geom = nz_coast_feature.geometry()
all_x = []
all_y = []
for i in range(geom.GetGeometryCount()):
    # Read ring geometry and create path
    r = geom.GetGeometryRef(i)
    x = [r.GetX(j) for j in range(r.GetPointCount())]
    y = [r.GetY(j) for j in range(r.GetPointCount())]
    all_x += x
    all_y += y
coastpoly = mpath.Path(np.column_stack((all_x,all_y)))
paths.append(coastpoly)

os.chdir('catalogs_rsqsim')


# %% Make Plot

cmap=plt.get_cmap('viridis')

#Select ruptures to plot

#rupture_list=[58852,59770,60923,62583,63724,64931] #Suggest list for Akatore Fault

rupture_list=[97773,115477,130551,145325,162209,176505] #Suggest list for Pisa Fault
#THE SUBPLOT MUST CORRESPOND TO THE LENGTH OF THE RUPTURE LIST
fig, axes = plt.subplots(3,2,sharex=True,layout='constrained',figsize=(6.5,11)) 

kk=0
fntsize=7
#LABEL LIST MUST CORRESPOND TO SUBPLOT SIZE
label_list=('(a) Mw ','(b) Mw ','(c) Mw ','(d) Mw ','(e) Mw ','(f) Mw ')

for ii in range(len(axes)):   
    for jj in range(len(axes[0])):
    # Add paths as patches to axes
        for coastpoly in paths:
            coast_patch = mpatches.PathPatch(coastpoly, \
                facecolor='grey', edgecolor='black', alpha=0.3)
            axes[ii][jj].add_patch(coast_patch)
            
        #index patches associated with each event
        rup_indx=np.where(all_rupture_events==rupture_list[kk]) 
        patches_indx=all_rupture_patches[rup_indx]
    
        for nn in range(len(patches_indx)):
        
            #index coordinates for each patch
            x_points=(otago_faults.loc[patches_indx[nn],0],otago_faults.loc[patches_indx[nn],3],otago_faults.loc[patches_indx[nn],6])
            y_points=(otago_faults.loc[patches_indx[nn],1],otago_faults.loc[patches_indx[nn],4],otago_faults.loc[patches_indx[nn],7])
        
            patch_disp=all_rupture_slip[rup_indx[0][nn]] #patch displacement
        
            #plot patch
            fault_patch=plt.Polygon([[x_points[0],y_points[0]],
                             [x_points[1],y_points[1]],
                             [x_points[2],y_points[2]],
                             [x_points[0],y_points[0]]],color=cmap(patch_disp/(8),0.3))
        
            #add patch to plot colored by displacement
            axes[ii][jj].add_patch(fault_patch)
        
        #Plot NZCFM
        for pp in range(nzcfm_layer.GetFeatureCount()):
            nzcfm_feature=nzcfm_layer.GetFeature(pp)
            nzcfm_geometry=nzcfm_feature.GetGeometryRef()
            nzcfm_x=np.empty([nzcfm_geometry.GetPointCount(),1])   
            nzcfm_y=np.empty([nzcfm_geometry.GetPointCount(),1])
            for rr in range(nzcfm_geometry.GetPointCount()):
                nzcfm_x[rr],nzcfm_y[rr],z=nzcfm_geometry.GetPoint(rr)
            axes[ii][jj].plot(nzcfm_x,nzcfm_y,'r-')

        #Set axes limits, corresponds to Akatore Fault in NZCM
        #axes[ii][jj].set_xlim([1.36*10e5, 1.415*10e5])
        #y_lim=([4.85*10e5, 4.92*10e5])
        
        #Set axes limits, corresponds to Pisa Fault in NZCM
        axes[ii][jj].set_xlim([1.26*10e5, 1.315*10e5])
        y_lim=([4.98*10e5, 5.043*10e5])
        
        axes[ii][jj].set_ylabel('NZTM_Y',fontsize=fntsize)
        axes[ii][jj].set_ylim(y_lim[0], y_lim[1]) 
        axes[ii][jj].set_aspect('equal', adjustable='box')  
        axes[ii][jj].set_xlabel('NZTM_X',fontsize=fntsize) 
        axes[ii][jj].set_ylabel('NZTM_Y',fontsize=fntsize)
        
        #set plot title
        event_label=label_list[kk]
        event_mw=str(round(otago_rsqsim_catalog.iloc[rupture_list[kk],3],2))
        event_yr=str(round(otago_rsqsim_catalog.iloc[rupture_list[kk],1]/(3600*365.25*24)))
        title_str=event_label+event_mw+', '+event_yr+' yrs'
        axes[ii][jj].set_title(title_str,fontsize=fntsize+1)
        
        #edit axis labels so numbers not so large
        x_tmp=axes[ii][jj].get_xticks().tolist()
        xx_tick=np.empty(len(x_tmp))

        for nn in range(len(x_tmp)):
            xx_tick[nn]=int(x_tmp[nn])/10e4

        axes[ii][jj].xaxis.set_major_locator(mticker.FixedLocator(x_tmp))
        axes[ii][jj].set_xticklabels(xx_tick, fontsize=fntsize)

        y_tmp=axes[ii][jj].get_yticks().tolist()
        yy_tick=np.empty(len(y_tmp))
        
        for nn in range(len(y_tmp)):
            yy_tick[nn]=int(y_tmp[nn])/10e5

        axes[ii][jj].yaxis.set_major_locator(mticker.FixedLocator(y_tmp))
        axes[ii][jj].set_yticklabels(yy_tick, fontsize=fntsize)
        axes[ii][jj].tick_params(axis='y', labelrotation = 90)
        
        kk=kk+1    

pcm = axes[ii][jj].pcolormesh(np.random.random((1, 1)),vmin=0, vmax=8)
fig.colorbar(pcm, ax=axes[ii][jj],location='bottom', shrink=0.4,label='rupture displacement (m)')


plt.show()  



