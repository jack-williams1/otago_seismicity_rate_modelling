## GENERATION AND ANALYSIS OF STOCHASTIC CATALOGS [RUNS IN MATLAB]

Generate and analyse 1 million year long stochastic event catalogs for Otago fault seismicity using the NZ CFM.

## RUN THIS FIRST

catalog_stochastic/fault_recurrence_parameters.m: Extracts NZ CFM attributes on Otago fault from 'OtagoRangeBasinFaults.xlsx' to derive their moment rate and recurrence model (i.e., fault specific annual rate of events >=M_min and magnitude probability distribution function). Here, the area of fault has been calculated from the function orb_fault_geometiees/faultgeometries.m and saved in 'OtagoRangeBasinFaults.xlsx.' If recurrence model is set to 4, uses slip rate data from the NZ NSHM 2022 Geodetic Model (Johnson et al 2022,2024). If recurrence model is set to 1-3, uses slip rate data from the NZ CFM (Seebeck et al 2022,2024, see also folder 'orb_fault_geometries'). Recurrence model options 1-3 represent different aperiodicites [0.5, 2, 4] for the Brownian Passage Time and Weibull catalogs. Recurrence parameters are stored in folders model1, model2, or model3 depending on the aperiodicity that is selected (or folder model4 if geodetic slip rates are selected). If creating recurrence models for geodetic rates, run nshm_geodeticrates/extract_otago_sections.py first.

All recurrence models use NZ CFM geometry model. Recurrence models are then developed assuming Gutenberg-Richter and characteristic on-fault magnitude frequency distributions (Youngs and Coppersmith 1985, Convertito 2006) using the function "characteristic_magnitude_YC1985.m," originally written by Katsu Goda (e.g., Goda and Sharipov 2021). Code generates inputs 'orb_fault_combined,' 'orb_fault_parameters,' and 'orb_faults_segmented,' which are *REQUIRED* to generate the various stochastic event catalogs

## THEN GENERATE CATALOGS!

catalog_stochastic/catalog_bpt.m: Generates stochastic event catalogs for segmented-char, combined-char, and combined-GR recurrence models and a Brownian Passage Time (BPT) renewal process. Uses a Brownian Relaxation Oscillator (Matthews et al 2002) to realize the BPT model. Need to select aperiodicity of BPT model that is being simulated, and whether using NZ CFM or geodetic model slip rates. Returns catalog as 'catalog_BPT' in relevant folders model1-model4. Also option to save the evolution of the load state variable for each fault within the catalogs.

catalog_stochastic/catalog_poisson.m: Generates stochastic event catalog for Poisson earthquake interevent times. Note there is no strict need to run for different aperiodicities, as it is implicit that this = 1 in a Poisson process. Hence, can generate catalog once, and then copy .mat file to other model folders. Does need to be run separately for Model 4 as that uses a different slip rate input. Catalogs run for all on-fault MFDs.

catalog_stochastic/catalog_weibull.m: Generates stochastic event catalog for Weibull earthquake interevent times and given aperiodicity. If aperiodicity set to 0.5, Weibull parameters derived from methods of moment estimator. If aperiodicity 2 or 4, derived from fminunc function. Catalogs run for all on-fault MFDs.

## ONCE CATALOGS RUN, CAN THEN ANALYSE

catalog_stochastic/catalog_50year_samp.m: Extract 50-year period samples from the stochastic event catalog, which can be used to compare to 50-year NSHM URZ forecasts in all_catalog_comparison.m

catalog_stochastic/catalog_analysis.m: Analyses catalog selected by aperiodicity (set by 'model_opt') and renewal process (set by 'cat_option'). Analyses include: (1) evaluating moment rate of entire catalog, and by fault, (2) separately plotting catalog MFD for different on-fault MFDs, (3) plotting catalog interoccurrence times, and (4) subsampling 70-year periods from catalog, and plotting overall and sample's MFD. !Important! saves subsamples and moment rate calcs, which are needed for running all_catalog_comparison.m and plot_emptycounts.m

catalog_stochastic/catalog_analysis_byfault.m: Analyses catalog events for individual fault source (set by 'ff' and 'select_opt') and aperiodicity. Can simultaneously analyse catalog events for multiple renewal processes. Analyses include: (1) plot fault's MFD in catalogs with different renewal processes, (2) plot fault's MFD against recurrence model's median theoretical MFD, (3) plot all of the fault's theoretical MFDs (from exploring different b-value and Mmax combinations), and (4) plot theoretical and catalog interevent times and hazard functions (following Yakovlev et a 2006) for fault. 

catalog_stochastic/plot_bpt_loadstate.m: Plot evolution of load state variable of the BPT catalog for a given fault and catalog time window. Single plot shows load state variable for both G-R and characteristic catalogs. Option to show multiple plots for catalogs with different fault aperiodicities (is equivalent to Figure S6 in manuscript). 

catalog_stochastic/plot_catalog_comparison.m: Provides capability to plot multiple stochastic event catalogs together. Here is configured to compare the overall catalog's MFD for difference renewal plots, moment rate for individual faults, MFD for single fault, and mean and standard deviation of fault recurrence interval in Weibull catalogs. Equivalent to Figure S7 in manuscript.

catalog_stochastic/plot_emptycounts.m: For each stochastic event catalog, plot number of 70-year catalog samples with zero events. Figure arranged, so each plot shows the same aperiodicity value, but different renewal processes and on-fault MFDs. Equivalent to Figure 8 in manuscript. Requires catlaog_analysis.m to have been run for all catalogs prior to plotting.

catalog_stochastic/plot_intertimes.m: Plots interevent times and hazard functions for fault ff in different stochastic catalogs. Used for Figure S8 in manuscript.

catalogs_stochastic/plot_recurrence_model.m: Compare median MFDs for multiple faults from Young and Coppersmith (1985) recurrence models (e.g., Figure 2c in manuscript), and compare all MFD's for a single fault from exploring all b-value M-max combinations (e.g., Figure S5b).

## REFERENCE

-Convertito, V., Emolo, A., & Zollo, A. (2006). Seismic-hazard assessment for a characteristic earthquake scenario: an integrated probabilistic–deterministic method. Bulletin of the Seismological Society of America, 96(2), 377-391.

-Goda, K., & Sharipov, A. (2021). Fault-source-based probabilistic seismic hazard and risk analysis for Victoria, British Columbia, Canada: A case of the Leech River Valley Fault and Devil’s Mountain Fault System. Sustainability, 13(3), 1440.

-Johnson, K., L. Wallace, J. Maurer, I. Hamling, C. Williams, C. Rollins, M. Gerstenberger, and R. Dissen (2022). Geodetic Deformation Model for the 2022 Update of the New Zealand National Seismic Hazard Model, GNS Science Report, GNS Science, Lower Hutt, New Zealand

-Johnson, K. M., Wallace, L. M., Maurer, J., Hamling, I., Williams, C., Rollins, C., ... & Van Dissen, R. (2024). Inverting geodetic strain rates for slip deficit rate in complex deforming zones: An application to the New Zealand plate boundary. Journal of Geophysical Research: Solid Earth, 129(3), e2023JB027565.
-Matthews, M. V., Ellsworth, W. L., & Reasenberg, P. A. (2002). A Brownian model for recurrent earthquakes. Bulletin of the Seismological Society of America, 92(6), 2233-2250.

-Seebeck, H., R. J. Van Dissen, N. J. Litchfield, P. M. Barnes, A. Nicol, R. M. Langridge, D. J. A. Barrell, P. Villamor, S. M. Ellis, M. S. Rattenbury, et al. (2022). New Zealand Community Fault Model–version 1.0, GNS Science report 2021/57, GNS Science, Lower Hutt, New Zealand, 97 pp., doi: 10.21420/GA7S-BS61.

-Seebeck, H., Van Dissen, R., Litchfield, N., Barnes, P. M., Nicol, A., Langridge, R., ... & Lee, J. (2024). The New Zealand Community Fault Model–version 1.0: An improved geological foundation for seismic hazard modelling. New Zealand Journal of Geology and Geophysics, 67(2), 209-229.

-Tinti, S., & Mulargia, F. (1987). Confidence intervals of b values for grouped magnitudes. Bulletin of the Seismological Society of America, 77(6), 2125-2134.

-Yakovlev, G., Turcotte, D. L., Rundle, J. B., & Rundle, P. B. (2006). Simulation-based distributions of earthquake recurrence times on the San Andreas fault system. Bulletin of the Seismological Society of America, 96(6), 1995-2007.

-Youngs, R. R., & Coppersmith, K. J. (1985). Implications of fault slip rates and earthquake recurrence models to probabilistic seismic hazard estimates. Bulletin of the Seismological society of America, 75(4), 939-964.
