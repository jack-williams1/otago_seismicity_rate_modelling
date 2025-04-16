##### SELECT 2022 NSHM RUPTURE SETS THAT INTERSECT AND WITHIN SEARCH POLYGON ######
### THEN RETURN A GEOJSON FILE FOR RUPTURES THAT INTERSECT THE POLYGON  #######

#A paticipation rate is also calculated for each section. For rupture sections within the 
#polygon, this should be equal to the full partipation rate of that section. For rupture sections
#outside the polygon, it is the rate for ruptures that intersect the polygon only

#Can be run for inversion solution rates from individual logic tree branches (opt=1) or for 
#combined geologic (opt=2) or geodetic weighted average (opt =3) rates. For opt=2 or 3, then
#need to run combine_branches.py first to get weighed average rates. 

import json
import pathlib
import geopandas
import pandas as pd

from shapely.geometry import shape
from shapely.geometry.polygon import Polygon 

import numpy as np #added by JW
import os

from solvis import InversionSolution, CompositeSolution, FaultSystemSolution
from solvis.filter import FilterRuptureIds

def export_geojson(gdf, filename, **kwargs):
    print(f"Exporting to {filename}")
    f = open(filename, 'w')
    f.write(gdf.to_json(**kwargs))
    f.close()
 

logic_tree_branches = {
    "geol_b0.823": "SW52ZXJzaW9uU29sdXRpb246MTEzMDUz",
    "geol_b0.959": "SW52ZXJzaW9uU29sdXRpb246MTEzMDMy",
    "geol_b1.089": "SW52ZXJzaW9uU29sdXRpb246MTEzMDM5",
    "geod_b0.823": "SW52ZXJzaW9uU29sdXRpb246MTEzMDY3",
    "geod_b0.959": "SW52ZXJzaW9uU29sdXRpb246MTEzMDgw",
    "geod_b1.089": "SW52ZXJzaW9uU29sdXRpb246MTEzMDYz",
}



# ================ MODIFY THESE LINES ================ 

#Select if want rates from individual or combined branches
#Combined branches requires running combine_branches.py first
#1 = individual logic tree branch rates
#2 = combined weighted mean geologic rates
#3 = combined weighted mean geodetic rates

opt=1

#if opt = 1 select branch
branch = "geol_b0.823"

# set your path
work_dir = pathlib.Path("WORK/")
output_dir= pathlib.Path("OtagoFaults/")

polygon_filename = "orb_area_polygon.geojson"
    
# ================ MODIFY THESE LINES ================ 

#Load model
if opt==1:
    
    soln_id = logic_tree_branches[branch]
    # load individual logic tree inversion solution
    soln_filepath = work_dir / f"{soln_id}.zip"
    soln = InversionSolution.from_archive(soln_filepath)
    filepath=output_dir / f"{soln_id}.geojson"
    
else:
    #Need to use full logic tree solution
    import nzshm_model as nm
    slt = nm.get_model_version("NSHM_v1.0.4").source_logic_tree
    comp = CompositeSolution.from_archive("WORK/NSHM_v1.0.4_CompositeSolution.zip", slt)

    # get the crustal fault solutions
    soln: FaultSystemSolution = comp._solutions['CRU'] 
    if opt==2: #import files created in combine_patches.py
        rates_w=pd.read_csv('OtagoFaults/geo_rates_w.csv')
        filepath=output_dir / "geo_rates.geojson"
    elif opt==3:
        rates_w=pd.read_csv('OtagoFaults/ged_rates_w.csv') 
        filepath=output_dir / "ged_rates.geojson"

# load the search polygon to select ruptures

os.chdir("..") #Go up one level in directory
os.chdir("..") #Go up one level in directory
geojson = json.load(open('gis_files/orb_area_polygon.geojson'))
polygon: Polygon = shape(geojson['features'][0]['geometry'])

os.chdir("nshm_inversion/solvis") 


# get the rupture_ids that intersect and within the selected polygon
if opt==1:
    rupture_ids = FilterRuptureIds(soln).for_polygon(polygon)
else:
    rupture_ids = rates_w.rup_id
    
rupture_ids_list=rupture_ids.tolist()

# extract selected ruptures from solution model
rr = soln.model.fault_sections_with_rupture_rates
rup_info: list = rr[rr['Rupture Index'].isin(rupture_ids)]
rup_info_list=rup_info['Rupture Index'].tolist()


#get fault sections of intersecting ruptures
rup_fs = rr[rr['Rupture Index'].isin(rupture_ids)]
rup_fs_list=rup_fs['Rupture Index'].tolist()
rup_sec=rup_fs['section'].tolist() 
unique_rup_sec=np.unique(rup_fs['section']) #all unique fault sections in rup_fs

#If opt is 2 or 3, add combined weighted rates and save to csv file for ALL sections
if opt!=1:
    tmp_w_rates=np.zeros(len(rup_fs))
    
    for kk in range(len(rupture_ids)):
        tmp_indx1=np.where(np.isin(rup_fs['Rupture Index'],rates_w.rup_id[kk],))[0].tolist()
        tmp_w_rates[tmp_indx1]=rates_w.w_rate[kk]
    
    #Add weighted average rates to rup_fs and export to csv
    rup_fs["w_rates"] = tmp_w_rates

    
    if opt==2:
        rup_fs.to_csv("OtagoFaults/fault_sections_geo_w.csv") 
    elif opt==3:
        rup_fs.to_csv("OtagoFaults/fault_sections_ged_w.csv") 

tmp_indx2=[0]*len(unique_rup_sec)
sec_rate=[0]*len(unique_rup_sec) 

#loop through each unique section
for rid in range(unique_rup_sec.size):

    tmp_indx2[rid]=rup_sec.index(unique_rup_sec[rid])

    #for each section, find all the ruptures its associated with in rup_fs
    tmp_indx3 = np.where(np.isin(rup_fs['section'], unique_rup_sec[rid]))[0].tolist()
    rate=[0] * len(tmp_indx3) 

    for jj in range(len(tmp_indx3)):
        if opt==1:
            rate[jj]=rup_fs.iloc[tmp_indx3[jj],2]#Use rates from branch's inversion solution
        else:
            rate[jj]=rup_fs.iloc[tmp_indx3[jj],26]#Use weighted average rates calculated from combined_branches.py 
    
    #combine the section's rate for all ruptures it particpated in    
    sec_rate[rid]=sum(rate) 

   
#remove rows in rup_fs with duplicate sections
rup_fs_unique=rup_fs.iloc[np.array(tmp_indx2)]

#add section prticipation rates to dataframe
rup_fs_unique.insert(len(rup_fs_unique.columns),"Participation_Rate",sec_rate)

#convert section particpation rates to geojson file
gdf = geopandas.GeoDataFrame(rup_fs_unique, crs="EPSG:4326")

export_geojson(gdf,filepath) 

#Uncomment to write rupture surfaces to geojson files (only for geologic weighted mean)
#This is needed to run get_fault_rupture_patches.py

if opt==2:
    for rid in range(len(rupture_ids_list)):
        rupture_surface_gdf=soln.rupture_surface(rupture_ids_list[rid])
        rupture_surface_gdf=rupture_surface_gdf.assign(weighted_geologic_rate=rates_w.w_rate[rid]) 
        export_geojson(rupture_surface_gdf, filename=pathlib.Path("OtagoFaults/rupture_surfaces/rupture_surface_" + str(rupture_ids[rid]) + ".geojson"))

