#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 12:48:21 2025

-Combine rates of 3 geologic and 3x geodetic logic tree branches
-Select_ruptures_by_polygon_unscaled.py and get_section_area.py needs to be run first
-Run in scripts/OtagoFaults folder

returns rates in csv file

@author: jackwilliams
"""

import pandas as pd
import numpy as np
import os

os.chdir('..')

#Get list of all rupture ids in RSQSimFault polygon from get_section_area.py
all_ruptures=pd.read_csv('OtagoFaults/filtered_rupture_ids.csv',header=None)
all_ruptures.columns=["rup_id"]
num_branches=3

count=0

tmp_col=np.zeros(len(all_ruptures))

nb_weights=[0.24,0.53,0.23] #weights for each logic tree branch

tmp_series1={'rup_id':all_ruptures.rup_id, 'mag':tmp_col, 'w_mag':tmp_col, 'geo_rate_nb1':tmp_col, 'geo_rate_nb2':tmp_col, 'geo_rate_nb3':tmp_col, 
            'geol_w': tmp_col, 'ged_rate_nb1':tmp_col, 'ged_rate_nb2':tmp_col, 'ged_rate_nb3':tmp_col, 'ged_w': tmp_col} 
all_branches_rates=pd.DataFrame(data=tmp_series1)


for ii in range(2): #loop through geologic_dm and geodetic_dm logic tree branches

    if ii==0:

        branches=["SW52ZXJzaW9uU29sdXRpb246MTEzMDUz", #dm: geologic, n-b pair: 0.823-2.7
        "SW52ZXJzaW9uU29sdXRpb246MTEzMDMy", #dm: geologic, n-b pair: 0.959-3.4
        "SW52ZXJzaW9uU29sdXRpb246MTEzMDM5"] #%dm: geologic, n-b pair: 1.089-4.6
        count=0
        
    else:

        branches=["SW52ZXJzaW9uU29sdXRpb246MTEzMDY3", # %dm geodetic, n-b pair: 0.823-2.7
        "SW52ZXJzaW9uU29sdXRpb246MTEzMDgw", #dm geodetic, n-b pair: 0.959-3.4
        "SW52ZXJzaW9uU29sdXRpb246MTEzMDYz"] #dm geodetic, n-b pair: 1.089-4.6
        count=4
    

    
    for jj in range(len(branches)):#loop through each branch
        
        branch_rate_info = pd.read_csv('OtagoFaults/'+branches[jj]+'_Mmod_Rtrim.csv',delimiter=",",header=0,
                                       usecols=[1,2,8,9],names=["rup_id","w_rate","mag","w_mag"])
        
        #for each branch, find if all_ruptures[kk] has a rate
        for kk in range(len(all_ruptures)):
            
            tmp=branch_rate_info[branch_rate_info['rup_id'] == all_ruptures.rup_id[kk]]
            
            if tmp.empty==True: #rup_id[kk] not included in branch's IFM solution
                all_branches_rates.iloc[kk,jj+3+count]=0
            else:  #take branch's IFM solution mag and rupture rate for rup_idk
                all_branches_rates.iloc[kk,1]=tmp.mag     
                all_branches_rates.iloc[kk,2]=tmp.w_mag 
                all_branches_rates.iloc[kk,jj+3+count]=tmp.w_rate
               
    for kk in range(len(all_ruptures)):  
        
        if ii==0:
            all_branches_rates.iloc[kk,jj+4]=sum((all_branches_rates.geo_rate_nb1[kk]*nb_weights[0],
                                         all_branches_rates.geo_rate_nb2[kk]*nb_weights[1],
                                         all_branches_rates.geo_rate_nb3[kk]*nb_weights[2]))
            
            tmp_w_rates=all_branches_rates.iloc[:,jj+4]
            file_name='geo_rates_w'
            
        else:
            all_branches_rates.iloc[kk,jj+count+4]=sum((all_branches_rates.ged_rate_nb1[kk]*nb_weights[0],
                                         all_branches_rates.ged_rate_nb2[kk]*nb_weights[1],
                                         all_branches_rates.ged_rate_nb3[kk]*nb_weights[2]))
            tmp_w_rates=all_branches_rates.iloc[:,jj+count+4]
            file_name='ged_rates_w'
         
    tmp_series2={'rup_id':all_branches_rates.iloc[:,0], 'mag':all_branches_rates.iloc[:,1], 'w_mag':all_branches_rates.iloc[:,2], 'w_rate':tmp_w_rates}   
    tmp_df2=pd.DataFrame(data=tmp_series2) #collate just ruptures ids, weighted magnitudes, and weighted rates
    tmp_df3=tmp_df2.loc[~(tmp_df2.w_rate==0),:]  #Remove rupture ids with zero rates 
    tmp_df3.to_csv('OtagoFaults/'+file_name+'.csv',index=False) 
  
        
os.chdir('scripts')        