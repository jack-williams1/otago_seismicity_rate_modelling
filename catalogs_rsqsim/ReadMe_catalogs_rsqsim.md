#### ANALYSIS OF OTAGO RSQSIM CATALOG #### [RUNS IN PYTHON AND MATLAB]

Analysis of 1 myr long catalog of earthquakes on Otago faults and surrounding 'edge' faults that was generated through the physics-based earthquake simulator RSQSim (Dietrich and Richards-Dinger 2010, Richards-Dinger and Dieterich 2012). Plots also allow comparison to fault seismicity forecast by the NZ NSHM 2022 Inversion Fault Model (IFM; Gerstenberger et al 2024).

!IMPORTANT! The Otago RSQSim Catalog must be downloaded first from: https://doi.org/10.5281/zenodo.13943280 and stored in the folder 'otago_rsqsim_catalog' to run these scripts. 

### DEPENDENCIES ###

NumPy: https://numpy.org
Pandas: https://pandas.pydata.org/docs/index.html
Matplotlib: https://matplotlib.org
Os: https://github.com/python/cpython/blob/3.13/Lib/os.py
gdal: https://pypi.org/project/GDAL/

### SCRIPTS ###

!RUN THESE FIRST!

-catalogs_rsqsim/rupture_patches.py: Analyses complete Otago RSQSim Catalog ('otago_rsqsim_catalog/otago_1e6yr_nospinup_catalogue.csv') and assess: (1) how much of each event's seismic moment is released on an Otago Fault (as specified in 'otago_rsqsim_catalog/otago_fault_list_20240428.csv'), and (2) how many faults participated in each event. This information is then store in a new catalog: otago_rsqsim_catalog/eqs_otago_1e6yr_partial_mw.csv. NOTE, in RSQSim, relationship between seismic moment and magnitude log10(M_0)=bMw+c is calculated using c=9.1. For consistency with the NZ NSHM 2022, where c=9.05 (Gerstenberger et al 2024), we re-calculate the magnitude of ALL events in the Otago RSQSim Catalog using c=9.05.
-catalogs_rsqsim/rupture_patches_by_fault.py: Downweight magnitude of all Otago RSQSim Catalog events as with rupture_patches.py. In addition: (1) store each event in a fault specific catalog that is returned in a .csv file in the folder 'otago_rsqsim_catalog/fault_catalog', (2) assess if event is 'surface rupturing' as indicated by having 10 surface rupture patches where the average single event displacement is >=1 m, and (3) test how many events are multi-fault, in which >=2 faults contribute a partial moment release of Mw>=6.7 (returned as file 'mflt_test.csv').

!CATALOG ANALYSIS!

-catalogs_rsqsim/plot_rupture_patches.py: Plot spatial distribution and slip of rupture patches selected in the 'rupture_list' attribute. Fault traces from the NZ CFM also shown (Seebeck et al 2024). Used for plots Figure S12 and S13 in manuscript, which show the rupture patches for events on the Akatore and Pisa faults.
-catalogs_rsqsim/rsqsim_allcatalog_analysis.m: Plots Otago RSQSim Catalog MFD, for both full event magnitudes and magnitudes after down weighting in rupture_patches.py. Derive catalog moment rate and proportion of multifault events. Take 1000x catalog samples of 50-, 70-, or 10,000-year duration, and plot subsamples with overall catalog mfd. Store analysis in 'catalog_rsqsim_statistics.' This is necessary for running 'all_catalog_analysis.m' and 'otago_rsqsim_catalog/rsqsim_byfaultanalysis.m.' Note this requires that the 50-, 70-, and 10,000-yr catalog sampling is run.
-catalogs_rsqsim/rsqsim_byfaultanalysis. Returns info about each fault's events in the Otago RSQSim Catalog (generated by 'rupture_patches_by_fault.py') to derive fault specific moment rates, variation of interevent times, proportion of multifault events, and average magnitude. Plots: (1) RSQSim-NZ NSHM 2022 IFM comparison of fault magnitudes (equivalent to Figure 4 in the manuscript), (2) comparisons of fault moment rate with the NZ CFM and coefficient of variation of events in the Otago RSQSim Catalog (equivalent to Figure S11 in manuscript), (3) fault specific RSQSim-NZ NSHM 2022 IFM MFD comparisons, and time window of fault events in the Otago RSQSim Catalog (equivalent to Figure 5 in the manuscript), and (4) MFD and time window of events simultaneously for multiple faults. These require that 'nshm_inversion/mfd_analysis/nshm_otago_inversion_results_by_fault.m' is run first. 

### REFERENCES ###

-Dieterich, J. H., & Richards-Dinger, K. B. (2010). Earthquake recurrence in simulated fault systems. In Seismogenesis and earthquake forecasting: The Frank Evison (Vol. 2, pp. 233–250). Springer.
-Gerstenberger, M. C., Van Dissen, R., Rollins, C., DiCaprio, C., Thingbaijim, K. K., Bora, S., ... & Williams, C. (2024). The seismicity rate model for the 2022 Aotearoa New Zealand national seismic hazard model. Bulletin of the Seismological Society of America, 114(1), 182-216.
-Richards‐Dinger, K., & Dieterich, J. H. (2012). RSQSim earthquake simulator. Seismological Research Letters, 83(6), 983-990.
-Seebeck, H., Van Dissen, R., Litchfield, N., Barnes, P. M., Nicol, A., Langridge, R., ... & Lee, J. (2024). The New Zealand Community Fault Model–version 1.0: An improved geological foundation for seismic hazard modelling. New Zealand Journal of Geology and Geophysics, 67(2), 209-229.