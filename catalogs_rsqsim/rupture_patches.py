#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 08:53:22 2024

Search through all events in Otago RSQSim catalog and return only the 
seismic moment released on faults in Otago, and not on edge faults

@author: wilja48p
"""

import numpy as np
import pandas as pd


#load catalog rupture patch and event info 
all_rupture_patches=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_patches.npy')#event patch
all_rupture_slip=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_slip.npy')#event patch slip
all_rupture_events=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_events.npy')#event index
#
#list of numbers that is the same as the length of events
rupture_list=np.unique(all_rupture_events)

#read input fault info
#otago fault names only. Check is consistent with NSHM IFM list 
otago_fault_list=pd.read_csv('otago_rsqsim_catalog/otago_fault_list_20240428.csv', usecols=[0])
#geometry info about faults
otago_faults = pd.read_csv('otago_rsqsim_catalog/otago_faults_2500_tapered_slip.flt', sep=" ", header=None)

#define fuction for calcutaing area of patch
def calculate_triangle_area(coord1, coord2, coord3):
    # Convert coordinates to numpy arrays
    A = np.array(coord1)
    B = np.array(coord2)
    C = np.array(coord3)
    
    # Calculate two vectors representing edges of the triangle
    AB = B - A
    AC = C - A

    # Calculate the cross product
    cross_product = np.cross(AB, AC)

    # Calculate the magnitude of the cross product
    area = 0.5 * np.linalg.norm(cross_product)

    return area

#%% Downweight magnitude of ruptures that are not entirely on Otago faults

rup_mw=[0]*len(rupture_list)  
num_fault=[0]*len(rupture_list)  
num_otago_fault=[0]*len(rupture_list)  
num_patches=[0]*len(rupture_list)  
#load catalog
otago_rsqsim_catalog=pd.read_csv('otago_rsqsim_catalog/otago_1e6yr_nospinup_catalogue.csv')#overall catalog

simu_count=0
tmp_count=1
mu=3e10 #crustal rigidity as used in NSHM
c=9.05 #scaling between moment and magnitude
b=1.5#scaling between moment and magnitude

#for kk in range(2000): 
#for kk in [57]:
for kk in rupture_list: 
    
    simu_count=simu_count+1
    
    #Keep track of loop through all ruptures
    if simu_count==tmp_count*5000:
        print('rupture count is ' + str(simu_count) + ' out of '+ str(len(rupture_list)))
        tmp_count=tmp_count+1
    
    #index patches associated with each event
    rup_indx=np.where(all_rupture_events==rupture_list[kk]) 
    patches_indx=all_rupture_patches[rup_indx]

    #get individual faults that ruptured in that event
    rup_faults=list(set(otago_faults.loc[patches_indx,12]))
    
    check=[0]*len(rup_faults)
    
    for jj in range(len(rup_faults)):
        if ((otago_fault_list.eq(rup_faults[jj])).any()).bool()==True:
            check[jj]=1
            
    if sum(check) ==len(rup_faults): 
        #all ruptured faults are within otago, can maintain original magnitude
        #I DO CHANGE SCALING BETWEEN MOMENT AND MAGNITUDE FROM 9.1 TO 9.05
        rup_mw[kk]= (np.log10(otago_rsqsim_catalog.iloc[kk,2])-c)/b
        num_fault[kk]=len(rup_faults) #record number of faults involved in rupture
        num_otago_fault[kk]=len(rup_faults) #record number of Otago faults involved in rupture
        
    elif  sum(check)>0:
    #some of the ruptured faults are in Otago
        
     #print(kk) #uncomment to identify problematic events 
     
     tmp2=[0] * len(patches_indx) 
     patch_area=[0] * len(patches_indx)
     patch_mo=[0] * len(patches_indx)
    
     for ii in range(len(patches_indx)):
        
        if (otago_fault_list.eq(otago_faults.loc[patches_indx[ii],12]).any()).bool()==True:
            
            #if patch's fault name is an Otago fault, determine its area and seismic moment
            vertex1=(otago_faults.loc[patches_indx[ii],0],otago_faults.loc[patches_indx[ii],1],otago_faults.loc[patches_indx[ii],2])
            vertex2=(otago_faults.loc[patches_indx[ii],3],otago_faults.loc[patches_indx[ii],4],otago_faults.loc[patches_indx[ii],5])
            vertex3=(otago_faults.loc[patches_indx[ii],6],otago_faults.loc[patches_indx[ii],7],otago_faults.loc[patches_indx[ii],8])
            patch_area[ii]=calculate_triangle_area(vertex1,vertex2,vertex3)
            patch_mo[ii]=patch_area[ii]*all_rupture_slip[rup_indx[0][ii]]*mu
        
        #recalculate rupture's moment and magnitude based on combined moment of Otago fault patches
        rup_mo=sum(patch_mo)
        
        #note catalog produced assuming 9.1 scaling between moment and magnitude, but 9.05 used in 2022 NSHM
        rup_mw[kk]=(np.log10(rup_mo)-c)/b
        num_fault[kk]=len(rup_faults)
        num_otago_fault[kk]=sum(check)
 
#%% Insert revised magnitudes to dataframe and write to csv file     
     
otago_rsqsim_catalog.insert(9, "partial_mw", rup_mw, True)    
otago_rsqsim_catalog.insert(10, "num_fault", num_fault, True)  
otago_rsqsim_catalog.insert(11, "num_otago_fault", num_otago_fault, True)      
otago_rsqsim_catalog.to_csv("otago_rsqsim_catalog/eqs_otago_1e6yr_partial_mw.csv")
