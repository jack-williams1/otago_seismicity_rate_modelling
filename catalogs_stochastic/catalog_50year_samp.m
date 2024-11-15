%% Get 50 year catalog samples to compare to the URZ %%

%Select recurrence model for Otago Faults to assess:
%1) NZCFM slip rates - Aperiodicity = 0.5
%2) NZCFM slip rates - Aperiodicity = 2
%3) NZCFM slip rates - Aperiodicity = 4
%4) Geodetic slip rates - Aperiodicity = 2

model_opt=2;

model_path=strcat('model',num2str(model_opt));
addpath(model_path)

load orb_fault_parameters
load('orb_fault_segmented','num_fault','orb_faults');
load('orb_fault_combined','num_comb','orb_faults_comb','c_area','w_slip_rate_m');
load('orb_catalog','catalog_duration');

%select catalog to analyse: (1) Poisson, (2) BPT, (3) Weibull
cat_option = 2;

if cat_option ==1
    load catalog_poisson
    StochasticEventCatalog=StochasticEventCatalog_poisson;
elseif cat_option ==2
    load('catalog_bpt','StochasticEventCatalog_bpt');
    StochasticEventCatalog=StochasticEventCatalog_bpt;
elseif cat_option ==3
    load('catalog_weibull');
    StochasticEventCatalog=StochasticEventCatalog_weibull;    
    
end

all_mag_range_GR=Mmin:dm:Mmax;

% StochasticEventCatalog
% 1    ) Event number
% 2    ) Simulation cycle
% 3    ) Occurrence time
% 4    ) Earthquake magnitude
% 5    ) NCZFM ID
% 6    ) Recurrence type (1 = characteristic versus 2 = exponential)
% 7    ) Reccurrece parameter
% 8    ) SED (for BPT only)

%% Sampling of stochastic catalogs

nsim=1000; %number of subsamples
catalog_duration=50;%set to equal NZ NSHM 2022 DSM forecast length
empty_count=zeros(length(StochasticEventCatalog),1);

sampleAnnualRate=zeros(nsim,length(all_mag_range_GR),2); sampleMoRate=zeros(nsim,2);
clock_old = cputime; simu_count=1;
    
rand_t_indx=rand(nsim,1).*(num_simu-catalog_duration);%nsmim random numbers between 0 and num_simu-cat_duration
    
for gg=1:nsim
    
    if gg == 2000*simu_count
        clock_new = cputime;
        disp(['Simulation cycle is: ',num2str(gg),' & Required time is: ',num2str(clock_new-clock_old)]);
        clock_old = clock_new;
            simu_count = simu_count + 1;
    end 
    
    %loop through seg-char and combined-char catalogs only
    for ss=1:2
    
        %index random 50 year catalog increment for random fault based catalog
        tmp_catalog_indx=find(StochasticEventCatalog{ss}(:,2)>=rand_t_indx(gg) & StochasticEventCatalog{ss}(:,2)<(rand_t_indx(gg)+catalog_duration));
     
        if isempty(tmp_catalog_indx)==1
            empty_count(ss)=empty_count(ss)+1;
        end
    
        tmp_catalog=StochasticEventCatalog{ss}(tmp_catalog_indx,:);%sample catalog
        sampleMoRate(gg,ss)=sum(10.^(1.5*StochasticEventCatalog{ss}(tmp_catalog_indx,4)+9.05))/catalog_duration; %sample catalog moment rate
    
        for ll=1:length(all_mag_range_GR)
            sampleAnnualRate(gg,ll,ss) = length(find(tmp_catalog(:,4) >= all_mag_range_GR(ll)))/catalog_duration;
        end
     
    clear tmp_catalog_indx
    
    end %end ss loop for each catalog
end%end gg loop for each sample

save('sampleAnnualRate_50','sampleAnnualRate','sampleMoRate')