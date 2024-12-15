#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 08:53:22 2024

Search through all events in Otago RSQSim catalog and return only the 
seismic moment released on faults in Otago, and not on edge faults

@author: wilja48p
\
"""
import numpy as np
import pandas as pd


#load catalog rupture patch and event info 
all_rupture_patches=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_patches.npy')#event patch
all_rupture_events=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_events.npy')#event index
all_rupture_times=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_slip_time.npy')#true event time
all_rupture_slip=np.load('otago_rsqsim_catalog/otago_1e6yr_nospinup_slip.npy')#true patch slip

#list of numbers that is the same as the length of events
rupture_list=np.unique(all_rupture_events)

#read input fault info
otago_fault_list=pd.read_csv('otago_rsqsim_catalog/otago_fault_list_20240428.csv', usecols=[0])#otago fault names only
#geometry info about faults
otago_faults = pd.read_csv('otago_rsqsim_catalog/otago_faults_2500_tapered_slip.flt', sep=" ", header=None)

otago_rsqsim_catalog=pd.read_csv('otago_rsqsim_catalog/otago_1e6yr_nospinup_catalogue.csv')#overall catalog

#create dictionary of  empty numpy arrays that can be called based on fault name
data_f={}

for ii in range(len(otago_fault_list)):
    tmp1=str(otago_fault_list.iloc[ii].tolist())
    tmp2=tmp1[2:len(tmp1)-2]
    
    data_f[tmp2]=np.empty((1,13))

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

#use indicies of rupture to check if each fault is surface rupturing

def surface_rupture_check(flt_geometry, flt_patches,rup_indx,all_rupture_slip):
         
    depth_cutoff=-10 #Cut off for minimum depth surface rupturing patch
    
    #use max function as patch depths are negative numbers
    min_depth=max([max(flt_geometry.loc[flt_patches,2].values),max(flt_geometry.loc[flt_patches,5].values),max(flt_geometry.loc[flt_patches,8].values)])
    
    #if a surface rupturing event, check average surface SED     
    if min_depth>=depth_cutoff:
        #index rupture ruptures within patches_indx list
        sur_rup_patches1=np.where(flt_geometry.loc[flt_patches,2].values>depth_cutoff)
        sur_rup_patches2=np.where(flt_geometry.loc[flt_patches,5].values>depth_cutoff)
        sur_rup_patches3=np.where(flt_geometry.loc[flt_patches,8].values>depth_cutoff)
        
        #unique array of each surface rupturing patch within patches_indx list
        sur_rup_patches=np.unique(np.concatenate((sur_rup_patches1,sur_rup_patches2,sur_rup_patches3),axis=1))
        #index of surface rupturing patches in rup_indx list
        sur_rup_indx=rup_indx[0][sur_rup_patches]
        mean_surface_sed=sum(all_rupture_slip[sur_rup_indx])/len(sur_rup_indx)
        
        if mean_surface_sed>1 and len(sur_rup_indx)>9:#1 m set as criteria for observing rupture in trench
        #1 m set as criteria, reasonably consistent with Otago SED, and 75% threshold in Biasi and Weldon (2013)
        #require 10 surface rupturing patches as there are many 1 patch events that rupture to the surface with >1 m displacement
            sur_rup=1 #surface rupturing earthquake
        else:
            sur_rup=0 #surface rupturing earthquake, but displacement too small to recognise in trench
    else:
        sur_rup=0 #blind rupture
        
    return sur_rup

#%% Downweight magnitude of ruptures that are not entirely on Otago faults
 
#load catalog
otago_rsqsim_catalog=pd.read_csv('otago_rsqsim_catalog/otago_1e6yr_nospinup_catalogue.csv')#overall catalog

simu_count=0
tmp_count=1
mu=3e10 #crustal rigidity as used in NSHM
c=9.05 #scaling between moment and magnitude
b=1.5#scaling between moment and magnitude
mflt_test=[0]*len(rupture_list) 
 

for kk in rupture_list: 
    
    simu_count=simu_count+1
    
    #Keep track of loop through all ruptures
    if simu_count==tmp_count*5000:
        print('rupture count is ' + str(simu_count) + ' out of '+ str(len(rupture_list)))
        tmp_count=tmp_count+1
    
    #index patches associated with each event in rupture_list
    rup_indx=np.where(all_rupture_events==rupture_list[kk]) 
    #index patches associated with each event in patches_list
    patches_indx=all_rupture_patches[rup_indx]
    
    #replace rupture time in otago_rsqsim_catalog df as the csv file has
    #overly rounded these event times. Instead obtain times directly from patches
    otago_rsqsim_catalog.loc[kk,"t0"]=all_rupture_times[rup_indx[0][0]]
    
    #get individual faults that ruptured in that event
    rup_faults=list(set(otago_faults.loc[patches_indx,12]))
    
    #if only 1 Otago fault ruptured, can maintain original magnitude
    #I DO CHANGE SCALING BETWEEN MOMENT AND MAGNITUDE FROM 9.1 TO 9.05
    if len(rup_faults)==1 and ((otago_fault_list.eq(rup_faults)).any()).bool()==True:  
        flt_mw=(np.log10(otago_rsqsim_catalog.iloc[kk,2])-c)/b
        
        #check if patches are surface rupturing
        sur_rup=surface_rupture_check(otago_faults,patches_indx,rup_indx,all_rupture_slip)
        
        tmp1=otago_rsqsim_catalog.loc[kk][0:9].values
        tmp2=np.append(tmp1,[flt_mw,sur_rup,1,1])
        
        data_f[rup_faults[0]]=np.vstack((data_f[rup_faults[0]],tmp2))
        
        mflt_test[kk]=1
    #if multi-fault rupture, find partial moment for each Otago fault
    else:
             
        check=[0]*len(rup_faults)
        
        #count number of Otago faults involved in rupture
        for jj in range(len(rup_faults)):
            if ((otago_fault_list.eq(rup_faults[jj])).any()).bool()==True:
                check[jj]=1
                
        flt_mo=[0]*len(rup_faults)
        flt_mw=[0]*len(rup_faults)
        
        #for each fault
        for jj in range(len(rup_faults)):
       
            #for each Otago fault
            if ((otago_fault_list.eq(rup_faults[jj])).any()).bool()==True:   
            
                #index rupture patches associated with Otago fault
                flt_patches_indx=np.where(otago_faults.loc[patches_indx,12]==rup_faults[jj])
                patch_area=[0] * len(flt_patches_indx[0])
                patch_mo=[0] * len(flt_patches_indx[0])
            
                #find seismic moment of Otago fault patch
                for ii in range(len(flt_patches_indx[0])):
                    tmp_indx=patches_indx[int(flt_patches_indx[0][ii])]
                    flt_rup_indx=np.where(tmp_indx==patches_indx) #index fault patch back to rup indx
                
                    vertex1=(otago_faults.loc[tmp_indx,0],otago_faults.loc[tmp_indx,1],otago_faults.loc[tmp_indx,2])
                    vertex2=(otago_faults.loc[tmp_indx,3],otago_faults.loc[tmp_indx,4],otago_faults.loc[tmp_indx,5])
                    vertex3=(otago_faults.loc[tmp_indx,6],otago_faults.loc[tmp_indx,7],otago_faults.loc[tmp_indx,8])
                    patch_area[ii]=calculate_triangle_area(vertex1,vertex2,vertex3)
                    patch_mo[ii]=patch_area[ii]*all_rupture_slip[rup_indx[0][flt_rup_indx[0][0]]]*mu
                
                #recalculate rupture's moment and magnitude based on the seismic moment of each fault's patch
                flt_mo[jj]=sum(patch_mo)
                #note catalog produced assuming 9.1 scaling between moment and magnitude, but 9.05 used in 2022 NSHM
                flt_mw[jj]=(np.log10(flt_mo[jj])-c)/b  
                
                #check if patches are surface rupturing
                sur_rup=surface_rupture_check(otago_faults, flt_patches_indx,rup_indx,all_rupture_slip)
                
                tmp1=otago_rsqsim_catalog.loc[kk].values 
                tmp2=np.append(tmp1,[flt_mw[jj],sur_rup,len(rup_faults),sum(check)]) #add partial mag, sur rup info, and multifault info
                
                data_f[rup_faults[jj]]=np.vstack((data_f[rup_faults[jj]],tmp2)) 
        
        #test whether the partial moment of each fault in a multifault rupture 
        #is equivalent to Mw 6.7 (i.e., the same as one section in the NSHM IFM)
            
        if sum(check)==len(rup_faults):
            mflt_test[kk]=sum(np.where(np.array(flt_mw)>6.7,1,0))
        
        #if rupture is crossing Otago to edge faults, require multifault rupture
        #to have a partial moment release of Mw6.7 outside Otago
        else:
            non_otago_mo=otago_rsqsim_catalog.iloc[kk,2]-sum(flt_mo)
            non_otago_indx=np.where(np.array(flt_mw)==0)[0][0]
            flt_mw[non_otago_indx]=(np.log10(non_otago_mo)-c)/b
            mflt_test[kk]=sum(np.where(np.array(flt_mw)>6.7,1,0))
    
#%%  Save output to fault specific csv files

mflt_test_df=pd.DataFrame(mflt_test)
mflt_test_df.to_csv('otago_rsqsim_catalog/mflt_test.csv')

#Create empty dataframe with correct column names
tmp_df1=pd.DataFrame(columns=otago_rsqsim_catalog.columns)
tmp_df1.insert(9, "partial_mw", True) 
tmp_df1.insert(10, "surface_rupture", True) 
tmp_df1.insert(11, "num_fault", True) 
tmp_df1.insert(12, "num_otago_fault", True) 

for ii in range(len(otago_fault_list)):
#for ii in [34]:
    tmp1=str(otago_fault_list.iloc[ii].tolist())
    tmp2=tmp1[2:len(tmp1)-2]
    
    #delete first row of each array where empty values
    data_f[tmp2]=np.delete(data_f[tmp2],0,0)

    tmp2_df = pd.DataFrame(data_f[tmp2],columns=tmp_df1.columns)
    
    #Save dataframe to csv file
    tmp_str=''.join(['otago_rsqsim_catalog/fault_catalog/', tmp2,'.csv'])
    tmp2_df.to_csv(tmp_str)     
     