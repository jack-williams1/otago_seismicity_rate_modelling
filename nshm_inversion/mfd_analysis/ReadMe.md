## ANALYSE MAGNITUDE FREQUENCY DISTRIBUTION OF OTAGO FAULTS FORECAST BY THE NZ NSHM 2022 INVERSION FAULT MODEL ## [RUNS IN MATLAB]

## PREAMBLE ##

These MATLAB scripts require the following outputs, which are generated from selecting Otago fault ruptures in the NZ NSHM 2022 Inversion Fault Model (IFM) using solvis:

-nshm_inversion/solvis/WORK/otago_ruptures_weighted.csv: Table of all ruptures that are within or intersect the Otago study area, and includes a down weighted rupture magnitude based on how much of the rupture is contained within the study area [created in select_ruptures_by_polygon.py]
-nshm_inversion/solvis/WORK/otago_ruptures_xxxxx.csv: Table of all ruptures that are within or intersect the Otago study area for a given logic tree branch xxxxx [created in select_ruptures_by_polygon.py]. In total, we include tables for the 6x logic tree branches that we explore in this work.
-nshm_inversion/solvis/WORK/otago_fault_sections.csv: List and info of all ruptures that a given fault section participates [created in selected_ruptures_by_polygon.py]. 
-nshm_inversion/solvis/WORK/otago_section_area; List and area of each Otago fault section [created in get_section_area.py].
-nshm_urz/nshm_urz_analysis.m: Outputs of URZ results
-otago_fault_list_20240428: List and NZ CFM ID of all Otago faults. Needs to be consistent with catalogs_rsqsim/otago_rsqsim_catalog/otago_fault_list_20240428.csv except that 'titricombined' is divided into 'Titri North' and 'Titri Central'

A description of the NZ NSHM 2022 IFM can be found out Gerstenberger et al (2024).

## SCRIPTS ##

-nshm_inversion/mfd_analysis/nshm_otago_inversion_results_analysis.m: Analysis of all Otago seismicity forecast by the NZ NSHM 2022 IFM. Includes comparing: (1) MFDs of all ruptures, and ruptures once scaled down to area in Otago, (2) MFD for different logic tree branches, and (3) the URZ negative binomial model forecast. Also quantifies proportion of multifault ruptures in the IFM solution, and the solution's moment rate. Saves analysis to nshm_otago_inversion_results.mat for use in other analysis
-nshm_inversion/mfd_analysis/nshm_otago_inversion_results_analysis_by_fault.m: Extracts information about ruptures that each Otago fault participates in. Collates in a fault specific .csv file that is stored in the folder '/by_fault.' Also stores MATLAB variable nshm_fault_stats.mat that compiles for each fault its moment rate, and weighed mean and range of rupture magnitudes. This is necessary for plotting Figure 4 in the manuscript in 'catalogs_rsqsim/rsqsim_byfaultanalysis.m' Option at end of script to plot MFD for a specific fault, or for multiple faults. Can also make plot that compares a fault's NZ CFM slip rate (Seebeck et al 2024) with its participation rate in the IFM, and where available, a paleoseismic constraint for its recurrence interval (Coffey et al 2024) (this is equivalent to Figure S19 in the manuscript).


## REFERENCES ## 

-Coffey, G. L., Rollins, C., Van Dissen, R. J., Rhoades, D. A., Gerstenberger, M. C., Litchfield, N. J., & Thingbaijam, K. K. (2024). Paleoseismic earthquake recurrence interval derivation for the 2022 revision of the New Zealand National Seismic Hazard Model. Seismological Research Letters, 95(1), 78-94.
-Seebeck, H., Van Dissen, R., Litchfield, N., Barnes, P. M., Nicol, A., Langridge, R., ... & Lee, J. (2024). The New Zealand Community Fault Model–version 1.0: An improved geological foundation for seismic hazard modelling. New Zealand Journal of Geology and Geophysics, 67(2), 209-229.
-Gerstenberger, M. C., Van Dissen, R., Rollins, C., DiCaprio, C., Thingbaijim, K. K., Bora, S., ... & Williams, C. (2024). The seismicity rate model for the 2022 Aotearoa New Zealand national seismic hazard model. Bulletin of the Seismological Society of America, 114(1), 182-216.