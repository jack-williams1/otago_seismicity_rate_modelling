#!/usr/bin/env python3
#RUN FROM FOLDER: otago_recurrence_modelling/nshm_inversion/solvis
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
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker

from shapely.geometry.polygon import Polygon 
from shapely.geometry import shape
from shapely.ops import transform 

import os
from osgeo import ogr

fault_select=['Akatore','Dunstan','Pisa']


fntsize=8

#define function to get fault geometry
def get_rup_polygons(fault_select):
    os.chdir("..") #Go up one level in directory
    os.chdir("..") #Go up one level in directory
    os.chdir('mfd_analysis/by_fault')
    
    #Go to individual fault inversion results created in nshm_otago_inversion_results_by_fault.m
    fault_file=fault_select+'.csv'
    fault_info=pd.read_csv(fault_file)
    
    #Only select ruptures from weighted mean of explored geologic logic tree branches
    branch_indx=np.where(fault_info.iloc[:,8]>0)
    
    #go to directory with rupture surfaces in created from get_section_area.py
    os.chdir("..") #Go up one level in directory
    os.chdir("..") #Go up one level in directory
    os.chdir("solvis/WORK/rupture_surfaces") 
    
    rupture_geom={}
    rup_rate=[0] * len(branch_indx[0])
    
    for ii in range(len(branch_indx[0])):
        rup_str='rupture_surface_'+str(fault_info.iloc[branch_indx[0][ii],0])+'.geojson'
        rupture_geom[ii]=json.load(open(rup_str)) #get rupture geometry
        rup_rate[ii]=fault_info.iloc[branch_indx[0][ii],8] #get rupture rate (from mean geologic_dm)
       
    os.chdir("..") #Go up one level in directory
    os.chdir("..") #Go up one level in directory
    os.chdir("scripts") #Go up one level in directory
    
    return(rupture_geom, rup_rate, fault_info)

# define a projection from lat/lon to NZTM
projected_crs = pyproj.CRS.from_epsg(2193)
projector = pyproj.Transformer.from_crs('epsg:4326', projected_crs, always_xy=True).transform


os.chdir("..") #Go up one level in directory
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

os.chdir("nshm_inversion/solvis/scripts") #Go back to scripts


fig, ax1 = plt.subplots(1,1,sharex=True,layout='constrained',figsize=(6.5,11)) 

#Plot coastline
for coastpoly in paths:
    coast_patch = mpatches.PathPatch(coastpoly, \
            facecolor='green', edgecolor='black', alpha=0.1)
    ax1.add_patch(coast_patch)
  
#Plot NZCFM       
for pp in range(nzcfm_layer.GetFeatureCount()):
     nzcfm_feature=nzcfm_layer.GetFeature(pp)
     nzcfm_geometry=nzcfm_feature.GetGeometryRef()
     nzcfm_x=np.empty([nzcfm_geometry.GetPointCount(),1])   
     nzcfm_y=np.empty([nzcfm_geometry.GetPointCount(),1])
     for rr in range(nzcfm_geometry.GetPointCount()):
         nzcfm_x[rr],nzcfm_y[rr],z=nzcfm_geometry.GetPoint(rr)
     ax1.plot(nzcfm_x,nzcfm_y,'r-') 

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
                   facecolor=rup_color[ii], alpha=(12+np.log10(rup_rate[jj]))/50)
            #plot patch
            ax1.add_patch(polygon_patch)

#Adjust axis limits as necessary depending on selected faults
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





