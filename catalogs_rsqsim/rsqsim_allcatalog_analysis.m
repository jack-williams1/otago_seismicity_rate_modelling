%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ANALYSIS OF RSQSIM CATALOG FOR OTAGO R&B FAULTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

c=9.05; b=1.5; %NZ NSHM 2022 scaling between moment and magnitude

%1 myr-long catalog produced by RSQSIM for Otago and surrounding faults (https://doi.org/10.5281/zenodo.13943280)

%This needs to be filtered using rupture_patches.py to find moment release on Otago faults first!
addpath('otago_rsqsim_catalog')
otago_rsqsimcatalog=readtable('eqs_otago_1e6yr_partial_mw.csv');
mflt_test=readtable('mflt_test.csv');
num_simu_rsqsim=1e6;   

%recalc magnitude in rsqsim catalog using NSHM moment-mag scaling
otago_rsqsimcatalog.mw=(log10(otago_rsqsimcatalog.m0)-c)/b;    
%moment for partial magnitude
otago_rsqsimcatalog.mo=10.^(otago_rsqsimcatalog.partial_mw.*b+c);

mw_max=max(otago_rsqsimcatalog.mw);
mw_max_w=max(otago_rsqsimcatalog.partial_mw);

dm=0.005; %should be consistent with dm in catalogs_stochastic/fault_recurrence_parameters
rsqsim_mag_range=5:dm:7.8;

% Seismic moment rate for whole catalog
rsqsim_MoRate_catalog = sum(otago_rsqsimcatalog.mo)/(num_simu_rsqsim);

%% Freq mag plot comparison for entire catalog

rsqsim_mfd_rate=zeros(length(rsqsim_mag_range),1); rsqsim_partial_mfd_rate=zeros(length(rsqsim_mag_range),1);

for ll = 1:length(rsqsim_mag_range)
     rsqsim_mfd_rate (ll) = length(find(otago_rsqsimcatalog.mw >= rsqsim_mag_range(ll)))/(num_simu_rsqsim);  
     rsqsim_partial_mfd_rate(ll)=length(find(otago_rsqsimcatalog.partial_mw >= rsqsim_mag_range(ll)))/(num_simu_rsqsim); 
end
  
%Make MFD Plot comparison

figure(1);

%Comparing all different fault-based catalogs
semilogy(rsqsim_mag_range,rsqsim_mfd_rate,'r-','LineWidth',1.5); hold on;
multi_fault_prop=zeros(2,1);

semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'k-','LineWidth',1.5); hold on;
legend('all catalog','Scaled by area within Otago');
   
 num_otago=length(find(otago_rsqsimcatalog.num_fault>0)); %total number of Otago events
    
 % Multifault stats
 num_multi_fault=length(find(otago_rsqsimcatalog.num_fault>1)); %total number of multifault ruptures including Otago faults
 num_multi_fault_otago=length(find(otago_rsqsimcatalog.num_otago_fault>1)); %total number of multifault ruptures exclusively in Otago faults
 num_multi_fault_edge=num_multi_fault-num_multi_fault_otago;  %total number of multifault ruptures on Otago and edge faults  
    
 %proportion of multi-fault events in rsqsim catalog to compare to
 %inversion (use 6.9, which is equivalent to one section patch in IFM
 tmp_indx=find(otago_rsqsimcatalog.mw>6.9 & otago_rsqsimcatalog.num_otago_fault>0);
 multi_fault_count1=length(find(otago_rsqsimcatalog.num_fault(tmp_indx,:)>1));
 multi_fault_prop(1)=multi_fault_count1/length(tmp_indx); %propotion of multifault events for all M>6.9 events

 multi_fault_count2=length(find(mflt_test.Var2>1));
 multi_fault_prop(2)=multi_fault_count2/length(tmp_indx); %propotion of multifault events for all M>6.9 events
    
 %create list of events m>6.9
 otago_rsqsimcatalog_mw69=[otago_rsqsimcatalog.mw(tmp_indx),otago_rsqsimcatalog.partial_mw(tmp_indx)];
 set(gca,'fontsize',12); axis([rsqsim_mag_range(1) rsqsim_mag_range(end) 10^-5 5*10^-1]); 
 xlabel('Magnitude'); ylabel('Annual frequency'); grid on; axis square;

%% Sampling of RSQSim catalog

%Select level of sampling
%1) 70 years for comparison to instrumental record
%2) 50 years for comparison to DSM forecast
%3) 10,000 years for comparison to paleoseismic record
samp_opt=3;

if samp_opt==1
    catalog_duration=70;
elseif samp_opt==2
    catalog_duration=50; 
elseif samp_opt==3
    catalog_duration=1e4; 
end

nsim=1000; %number of subsamples
empty_count=0; 

sampleAnnualRate_rsqsim=zeros(1000,length(rsqsim_mag_range)); 
clock_old = cputime; simu_count=1; sampleMoRate_rsqsim=zeros(nsim,1); 

rand_t_indx=rand(nsim,1).*(num_simu_rsqsim-catalog_duration);%nsmim random numbers between 0 and num_simu-cat_duration
rsqsim_catalog_t=otago_rsqsimcatalog.t0/(365.25*24*3600);%convert catalog time to years

for gg=1:nsim
    
    if gg == 200*simu_count
        clock_new = cputime;
        disp(['Simulation cycle is: ',num2str(gg),' & Required time is: ',num2str(clock_new-clock_old)]);
        clock_old = clock_new;
            simu_count = simu_count + 1;
    end 
    
     %index random catalog increment from rsqsim catalog
     tmp_catalog_indx=find(rsqsim_catalog_t>=rand_t_indx(gg) & rsqsim_catalog_t<(rand_t_indx(gg)+catalog_duration));

     if isempty(tmp_catalog_indx)==1
         empty_count=empty_count+1;
     end

     tmp_catalog_rsqsim_mw=otago_rsqsimcatalog.partial_mw(tmp_catalog_indx);%sample catalog partial magnitudes
     
     %record and remove events with 0 magnitude
     if max(tmp_catalog_rsqsim_mw)==0
           empty_count=empty_count+1; sampleMoRate_rsqsim(gg)=0;
     else
           tmp_catalog_rsqsim_mw=tmp_catalog_rsqsim_mw(tmp_catalog_rsqsim_mw>0);
           sampleMoRate_rsqsim(gg)=sum(10.^(tmp_catalog_rsqsim_mw.*b+c))/catalog_duration;
     end
     
     for ll=1:length(rsqsim_mag_range)
        sampleAnnualRate_rsqsim(gg,ll) = length(find(tmp_catalog_rsqsim_mw >= rsqsim_mag_range(ll)))/catalog_duration;
     end
     
end

if samp_opt==1
    sampleAnnualRate_rsqsim_70=sampleAnnualRate_rsqsim; sampleMoRate_rsqsim_70=sampleMoRate_rsqsim; empty_count1=empty_count;
elseif samp_opt==2
    sampleAnnualRate_rsqsim_50=sampleAnnualRate_rsqsim; sampleMoRate_rsqsim_50=sampleMoRate_rsqsim; empty_count2=empty_count;
elseif samp_opt==3
    sampleAnnualRate_rsqsim_10kyr=sampleAnnualRate_rsqsim; sampleMoRate_rsqsim_10kyr=sampleMoRate_rsqsim; empty_count3=empty_count;
end

otago_indx=find(otago_rsqsimcatalog.partial_mw>0);

%Get average rate of M>5 events in the catalogs, plus std of this rate from catalog subsamples
m5_rate=[length(otago_indx)/num_simu_rsqsim std(sampleAnnualRate_rsqsim(:,1))];

%% Plot overall catalog MFD with sub samples

samp_opt=3;

if samp_opt==1; sampleAnnualRate_rsqsim=sampleAnnualRate_rsqsim_70;
elseif samp_opt==2; sampleAnnualRate_rsqsim=sampleAnnualRate_rsqsim_50;
elseif samp_opt==3; sampleAnnualRate_rsqsim=sampleAnnualRate_rsqsim_10kyr;
end

figure(5);
cat_col=vertcat([0.6 0.6 0.6 0.5]);

semilogy(rsqsim_mag_range,squeeze(sampleAnnualRate_rsqsim),'Color',cat_col(1,:),'LineWidth',0.5); hold on; %MFD from catalog samples
p3=semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'k-','LineWidth',1.5); hold on; %MFD from overall catalog

set(gca,'fontsize',12);  axis([rsqsim_mag_range(1) rsqsim_mag_range(end) 10^-5 5*10^0]); 

xlabel('Magnitude'); ylabel('Annual frequency'); grid on; axis square;

%% Save outputs

save('catalog_rsqsim_statistics','sampleAnnualRate_rsqsim_50','sampleMoRate_rsqsim_50','sampleAnnualRate_rsqsim_70','sampleMoRate_rsqsim_70', ...
    'sampleAnnualRate_rsqsim_10kyr','sampleMoRate_rsqsim_10kyr','rsqsim_mfd_rate','rsqsim_partial_mfd_rate','otago_rsqsimcatalog_mw69',...
   'rsqsim_mag_range','rsqsim_MoRate_catalog');


