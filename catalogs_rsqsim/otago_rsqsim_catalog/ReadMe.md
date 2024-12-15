### Place the Otago RSQSim Catalog (https://doi.org/10.5281/zenodo.13943280) HERE ####

The following files are required to run the codes stored in "catalogs_rsqsim:"

-otago_1e6yr_nospinup_catalogue.csv: catalog of events in the Otago RSQSim Catalog
-otago_1e6yr_nospinup_patches.npy: index for which fault each rupture patch is associated with 
-otago_1e6yr_nospinup_slip.npy: amount of slip for each event's rupture patches
-otago_1e6yr_nospinup_events.npy: index for which event each rupture patch is associated with
-otago_1e6yr_nospinup_slip_time.npy: index for the timings of rupture on each event's rupture patches
-otago_faults_2500_tapered_slip.flt: look-up table for coordinates of each fault's rupture patches 

[Already Included]
-otago_fault_list_20240428.csv: look-up table for the names of Otago faults included in the catalog (and not 'edge' faults'). This should be self consistent with fault names in nshm_inversion/mfd_analysis/otago_fault_list_20240428.csv except that 'Titri North' and 'Titri Central' are almagamated into 'titricombined'
