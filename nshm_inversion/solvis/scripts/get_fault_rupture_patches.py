#!/usr/bin/env python3
#Ensure necessary packages that are imported below are installed in python environment
#RUN FROM FOLDER nshm_inversion/solvis/scripts

# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:32:36 2024

@author: jackwilliams
"""

import numpy as np
import pandas as pd
import json
import pyproj #added by JW

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker

from shapely.geometry.polygon import Polygon 
from shapely.geometry import shape
from shapely.ops import transform 

import geopandas as gpd
import os

fault_select=['Akatore','Dunstan','Pisa']

fntsize=8

#define function to get fault geometry
def get_rup_polygons(fault_select):
    
    os.chdir("..")

    all_fault_sec_info=pd.read_csv('OtagoFaults/fault_sections_geo_w.csv')
    #select all ruptures that fault participates in
    fault_rups=all_fault_sec_info[all_fault_sec_info['ParentName']==fault_select]
        
    rupture_geom={}
    rup_rate=np.zeros((len(fault_rups),1))
    
    #Get rupture surfaces, as created in get_rupture_traces.py, and mean geologic rupture rates
    for ii in range(len(fault_rups)):
        rup_str='OtagoFaults/rupture_surfaces/rupture_surface_'+str(fault_rups.iloc[ii,2])+'.geojson'
        rupture_geom[ii]=json.load(open(rup_str)) #get rupture geometry
        rup_rate[ii]=fault_rups.iloc[ii,28] #get rupture rate (from mean geologic_dm)
     
    os.chdir("scripts") #Go up one level in directory    
    return(rupture_geom, rup_rate, fault_rups)

# define a projection from lat/lon to NZTM
projected_crs = pyproj.CRS.from_epsg(2193)
projector = pyproj.Transformer.from_crs('epsg:4326', projected_crs, always_xy=True).transform


os.chdir("..") #Go up one level in directory
os.chdir("..") #Go up one level in directory
os.chdir("..") #Go up one level in directory

nzcfm= gpd.read_file("gis_files/NZ_CFM_v1_0.shp")

nz_coast=gpd.read_file("gis_files/nz-coastline-polygon_nztm.shp")

os.chdir("nshm_inversion/solvis/scripts") #Go back to scripts

fig, ax1 = plt.subplots(1,1,sharex=True,layout='constrained',figsize=(6.5,11)) 
  
#Plot coast line 
nz_coast.plot(ax=ax1,facecolor="green",edgecolor="black",alpha=0.1)  
#Plot NZCFM       
nzcfm.plot(ax=ax1, linewidth=1,color='red')

rup_color=['orange','blue','purple']

#Plot each fault
for ii in range(len(fault_select)):

    #extract rupture geometry
    rupture_geom,rup_rate,tmp_file=get_rup_polygons(fault_select[ii])
    
    #for each rupture
    for jj in range(len(rupture_geom)):
        tmp_rup_geom=rupture_geom[jj]
    
        #For each section
        for kk in range(len(tmp_rup_geom['features'])):
        
            polygon: Polygon = shape(tmp_rup_geom['features'][kk]['geometry'])
        
            #project lat/lon coordinates of polygon to NZTM
            projected_polygon = transform(projector, polygon)
            x,y = projected_polygon.exterior.xy
        
            polygon_patch = mpatches.Polygon(np.stack([x,y],axis=1),edgecolor='black',
                   facecolor=rup_color[ii], alpha=((12+np.log10(rup_rate[jj]))/50)[0])
            #plot patch
            ax1.add_patch(polygon_patch)
 
ax1.set_xlim([1.245*10e5, 1.45*10e5])
y_lim=([4.85*10e5, 5.15*10e5])
ax1.set_ylim(y_lim[0], y_lim[1]) 

ax1.set_aspect('equal', adjustable='box')  

#edit axis lavels so numbers not so large
x_tmp=ax1.get_xticks().tolist()
xx_tick=np.empty(len(x_tmp))

for nn in range(len(x_tmp)):
    xx_tick[nn]=int(x_tmp[nn])/10e4

ax1.xaxis.set_major_locator(mticker.FixedLocator(x_tmp))
ax1.set_xticklabels(xx_tick, fontsize=fntsize)

y_tmp=ax1.get_yticks().tolist()
yy_tick=np.empty(len(y_tmp))

for nn in range(len(y_tmp)):
    yy_tick[nn]=int(y_tmp[nn])/10e5

ax1.yaxis.set_major_locator(mticker.FixedLocator(y_tmp))
ax1.set_yticklabels(yy_tick, fontsize=fntsize)
ax1.tick_params(axis='y', labelrotation = 90)
ax1.set_xlabel('NZTM X') 
ax1.set_ylabel('NZTM Y') 

os.chdir('..')
#Export figure, adjust filename as necessary
plt.savefig('OtagoFaults/pisa_akatore_dunstan_rup_maps.jpg')
os.chdir('scripts')




