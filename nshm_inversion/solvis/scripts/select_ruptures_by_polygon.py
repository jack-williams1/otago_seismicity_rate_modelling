##### SELECT 2022 NSHM RUPTURE SETS THAT INTERSECT AND WITHIN SEARCH POLYGON ######
### THEN RETURN THE RUPTURES'S MOMENT MAGNITUDE THAT IS RELEASED WITHIN THE POLOGON ###

import nzshm_model as nm
import json
import pathlib

from shapely.geometry import shape
from shapely.geometry.polygon import Polygon 

from shapely.ops import transform #added by JW
import pyproj #added by JW
import numpy as np #added by JW
import os #added by JW

from solvis import CompositeSolution, FaultSystemSolution, export_geojson

#Set logic tree exploration option
# 1 = weighted mean rates from all branches
# 2 = select certain branch

logic_tree_opt=1

#select branch. 
if logic_tree_opt==2:
    #uncomment to print branch options
    #ref_table=open("nshm_logictreebranch_lookuptable.txt",'r')
    #branch_ref_table=ref_table.read()
    #print(branch_ref_table) 
    branch="U2NhbGVkSW52ZXJzaW9uU29sdXRpb246MTIwNzUz"


# load the search polygon to select ruptures
os.chdir("..") #Go up one level in directory
os.chdir("..") #Go up one level in directory
geojson = json.load(open('gis_files/orb_area_polygon.geojson'))
polygon: Polygon = shape(geojson['features'][0]['geometry'])

# load the model
os.chdir('nshm_inversion')
slt = nm.get_model_version("NSHM_v1.0.4").source_logic_tree
comp = CompositeSolution.from_archive("NSHM_v1.0.4_CompositeSolution.zip", slt)

#Go back to solvis directory
os.chdir('solvis')

# get the crustal fault solutions
fss: FaultSystemSolution = comp._solutions['CRU']  # NB this API call will change in the future

# get the rupture_ids that intersect and within the selected polygon
rupture_ids: list = fss.get_ruptures_intersecting(polygon)
rupture_ids_list=rupture_ids.tolist()

#get rupture rates that are within search polygon based on logic tree opt

if logic_tree_opt==1: #Use weighted rupture rates from all logic tree branches
    rr = fss.ruptures_with_rupture_rates
    rup_info: list = rr[rr['Rupture Index'].isin(rupture_ids)]
    rup_info_list=rup_info['Rupture Index'].tolist()

    
elif logic_tree_opt==2: #Use rupture rates from specific logic tree branch
    cr = fss.composite_rates 
    #select rupture rates by logic tree branch
    rup_info=cr[(cr["Rupture Index"].isin(rupture_ids)) & (cr.solution_id==branch)]
    rup_info_list=rup_info['Rupture Index'].tolist()  
    #select only rupture ids that are in logic tree branch solution
    rupture_ids=rupture_ids[rupture_ids.isin(rup_info_list)]
    rupture_ids_list=rupture_ids.tolist()


#get fault sections of ruptures
#from line 285 at: https://github.com/GNS-Science/solvis-graphql-api/blob/181a42ee119f664a7fba16ed55c067073512935a/solvis_graphql_api/composite_solution/cached.py#L246
fsr = fss.fault_sections_with_rupture_rates

#get fault sections of intersecting ruptures
rup_fs = fsr[fsr['Rupture Index'].isin(rupture_ids)]
#rup_fs.to_csv("WORK/otago_fault_sections.csv") #Uncomment to write fault section csv file

#set dummy variables 
rupture_mw=[0] * len(rupture_ids) 
weighted_rupture_area=[0] * len(rupture_ids) 
weighted_rupture_mw=[0] * len(rupture_ids) 

def section_intersecting(fss: FaultSystemSolution, rupture_id: int):
    
    #get index of intersecting rupture
    #Rupture Index in rupture_ids and rup_info are NOT in the same order
    #so need seperate indicies
    tmp_idx1=rupture_ids_list.index(rupture_id)
    tmp_idx2=rup_info_list.index(rupture_id)
    
    #get surfaces of each rupture section
    rupture_surface_gdf=fss.rupture_surface(rupture_ids[tmp_idx1])
    rupture_sections_info=rupture_surface_gdf[['key_0','section','geometry','DipDeg','Magnitude','Area (m^2)']]
    rupture_polygon=rupture_sections_info['geometry'].tolist()
    section_dip=rupture_sections_info['DipDeg'].tolist()
    
    contained_indx=[0] * len(rupture_polygon) 
    weighted_section_area=[0] * len(rupture_polygon)
    
    #index the ruture sections individual geometry and area
    for ss in range(len(rupture_polygon)):
        
        section_polygon=rupture_polygon[ss] #define the geometry of each rupture section
        
        # define a projection from lat/lon to NZTM
        projected_crs = pyproj.CRS.from_epsg(2193)
        projector = pyproj.Transformer.from_crs('epsg:4326', projected_crs, always_xy=True).transform
        
        #project lat/lon coordinates of polygon to NZTM
        projected_polygon = transform(projector, section_polygon)
        
        #derive area of polgon in m^2
        tmp_area=projected_polygon.area #area of section polygon with no info about dip
        section_area=tmp_area/np.cos(np.radians(section_dip[ss])) #section polygon area once projected through fault dip
        
        #return value of 1 if rupture section entirely within search polygon
        if polygon.contains(section_polygon):
            contained_indx[ss]=1
        else:
            #return value of 0.5 if rupture section partly within polygon
            if section_polygon.intersects(polygon): 
                contained_indx[ss]=0.5
            
            #return value of 0 if rupture section not within polygon
            else:
                contained_indx[ss]=0
        
        #return section area weighted by whether it crosses polygon
        weighted_section_area[ss]=contained_indx[ss]*section_area
    
    #if rupture only partly within polygon, derive new area and magnitude
    if  min(contained_indx)<1:   
        weighted_rupture_area[tmp_idx2]=sum(weighted_section_area) #Total area of rupture within polygon
        rupture_area=(rupture_sections_info['Area (m^2)'].tolist())[0] #Total area of rupture
        
        rupture_mw[tmp_idx2]=(rupture_sections_info['Magnitude'].tolist())[0]#get total magnitude of rupture
        rupture_mo=10**(rupture_mw[tmp_idx2]*1.5+9.05) #Total seismic moment of rupture
        weighted_rupture_mo=rupture_mo*weighted_rupture_area[tmp_idx2]/rupture_area #seismic moment released within Otago polygon
        weighted_rupture_mw[tmp_idx2]=(np.log10(weighted_rupture_mo)-9.05)/1.5 #moment magnitude of polygon released within Otago
    
    #if rupture entirely within polygon no need to change area and magnitude
    else: 
        rupture_mw[tmp_idx2]=(rupture_sections_info['Magnitude'].tolist())[0]
        weighted_rupture_area[tmp_idx2]=rupture_sections_info['Area (m^2)'].tolist()[0]    
        weighted_rupture_mw[tmp_idx2]=rupture_mw[tmp_idx2]
    
    return weighted_rupture_area, rupture_mw, weighted_rupture_mw 


for rid in rupture_ids[:len(rupture_ids)]:  
    contain=section_intersecting(fss, rid)


#Reshape rows in contain to columns and add to new rup_info_w dataframe
rup_info_w=rup_info.assign(weighted_area=np.reshape(contain[0:1],(len(rupture_ids),1)))
rup_info_w=rup_info_w.assign(magnitude=np.reshape(contain[1:2],(len(rupture_ids),1)))
rup_info_w=rup_info_w.assign(weighted_magnitude=np.reshape(contain[2:],(len(rupture_ids),1)))

#Return rupture info in .csv file
if logic_tree_opt ==1:
    rup_info_w.to_csv("WORK/otago_ruptures_weighted.csv") #Uncomment to write csv file
elif logic_tree_opt ==2:
    filename_list=["WORK/otago_ruptures_" + branch + ".csv"]
    tmp_str=""
    filename=tmp_str.join(filename_list)
    rup_info_w.to_csv(filename)
    
#print(rup_info_w)
