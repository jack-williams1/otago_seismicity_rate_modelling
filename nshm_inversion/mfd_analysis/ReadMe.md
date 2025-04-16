## ANALYSE MAGNITUDE FREQUENCY DISTRIBUTION OF OTAGO FAULTS FORECAST BY THE NZ NSHM 2022 INVERSION FAULT MODEL ## [RUNS IN MATLAB]

Updated by JW 14/04/2025

## PREAMBLE ##

These MATLAB scripts require the following outputs first, which are generated from selecting Otago fault ruptures in the NZ NSHM 2022 Inversion Fault Model (IFM) using solvis:

-nshm_inversion/solvis/OtagoFaults/xxxxx_Mmod_Rtrim.csv: Table of all ruptures that are within or intersect the Otago study area for a given logic tree branch xxxxx [created in select_ruptures_by_polygon_unscaled.py]. In total, we include tables for the 6x logic tree branches that we explore in this work.
-nshm_inversion/solvis/OtagoFaults/geo_rates.csv: List and weighted mean of magnitudes and rates for ruptures from the 3 geologic_dm logic tree branches assessed in this study [created in combine_branches.py]. 
-nshm_inversion/solvis/OtagoFaults/ged_rates.csv: List and weighted mean of magnitudes and rates for ruptures from the 3 geodetic_dm logic tree branches assessed in this study [created in combine_branches.py]. 
-nshm_inversion/solvis/OtagoFaults/fault_sections_geo_w.csv: Information about the section of each rupture that is the inversion solution of the geologic dm branches, and intersects or is located within the Otago section area [created in get_rupture_traces_unscaled.py]. 
-nshm_inversion/solvis/OtagoFaults/fault_sections_ged_w.csv: Information about the section of each rupture that is the inversion solution of the geodetic dm branches, and intersects or is located within the Otago section area [created in get_rupture_traces_unscaled.py]. 
-nshm_inversion/solvis/OtagoFaults/otago_section_area; List and area of each Otago fault section [created in get_section_area.py].

The following scripts also needs to be run first

-nshm_urz/nshm_urz_analysis.m: Outputs of URZ results

This folder also contains:

-otago_fault_list_20240428: List and NZ CFM ID of all Otago faults. Needs to be consistent with catalogs_rsqsim/otago_rsqsim_catalog/otago_fault_list_20240428.csv except that 'titricombined' is divided into 'Titri North' and 'Titri Central'

A description of the NZ NSHM 2022 IFM can be found out Gerstenberger et al (2024).

## SCRIPTS ##

-nshm_inversion/mfd_analysis/nshm_otago_inversion_results_analysis.m: Analysis of all Otago seismicity forecast by the NZ NSHM 2022 IFM. Includes comparing: (1) MFDs of all ruptures, and ruptures once scaled down to area in Otago, (2) MFD for different logic tree branches, and (3) the URZ negative binomial model forecast. Also quantifies proportion of multifault ruptures in the IFM solution, and the solution's moment rate. Saves analysis to nshm_otago_inversion_results.mat for use in other analysis

-nshm_inversion/mfd_analysis/nshm_otago_inversion_results_analysis_by_fault.m: Extracts information about ruptures that each Otago fault participates in for either the geologic or geodetic deformation models. Collates in a fault specific .csv file that is stored in the folder '/by_fault_geo' or 'by_fault_ged.' Also stores MATLAB variable nshm_fault_stats.mat that compiles for each fault: its moment rate, and weighed mean and range of rupture magnitudes for these two deformation models. This is necessary for plotting Figure 4 in the manuscript in 'catalogs_rsqsim/rsqsim_byfaultanalysis.m' Option at end of script to plot MFD for a specific fault, or for multiple faults. Can also make plot that compares a fault's NZ CFM slip rate (Seebeck et al 2024) with its participation rate in the IFM, and where available, a paleoseismic constraint for its recurrence interval (Coffey et al 2024) (this is equivalent to Figure S20 in the manuscript).


## REFERENCES ## 

-Coffey, G. L., Rollins, C., Van Dissen, R. J., Rhoades, D. A., Gerstenberger, M. C., Litchfield, N. J., & Thingbaijam, K. K. (2024). Paleoseismic earthquake recurrence interval derivation for the 2022 revision of the New Zealand National Seismic Hazard Model. Seismological Research Letters, 95(1), 78-94.
-Seebeck, H., Van Dissen, R., Litchfield, N., Barnes, P. M., Nicol, A., Langridge, R., ... & Lee, J. (2024). The New Zealand Community Fault Model–version 1.0: An improved geological foundation for seismic hazard modelling. New Zealand Journal of Geology and Geophysics, 67(2), 209-229.
-Gerstenberger, M. C., Van Dissen, R., Rollins, C., DiCaprio, C., Thingbaijim, K. K., Bora, S., ... & Williams, C. (2024). The seismicity rate model for the 2022 Aotearoa New Zealand national seismic hazard model. Bulletin of the Seismological Society of America, 114(1), 182-216.
