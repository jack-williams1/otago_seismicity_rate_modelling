%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   PLOT AND ANALYSE POISSON, BPT, AND %%%%%
%%%%    WEIBULL STOCHASTIC CATALOGS    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1),'/catalogs_instrumental']);

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
cat_option = 3;

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

%% Result check: Seismic moment for each fault and MFD and intervent case:

% Seismic moment rate for whole catalog and individual faults
Mo_StochasticEventCatalog=cell(3,1); MoRate_StochasticEventCatalog_all=zeros(3,1);

for ss=1:length(StochasticEventCatalog)

    Mo_StochasticEventCatalog_all = 10.^(1.5*StochasticEventCatalog{ss}(:,4)+9.05);
    MoRate_StochasticEventCatalog_all(ss) = sum(Mo_StochasticEventCatalog_all)/(num_simu*t_limit);
    
    %Added by JW March 2025
    m5_rate(ss)=length(StochasticEventCatalog{ss})/num_simu;

    if ss==1
        num_sources=num_fault; source_id=orb_faults.OBJECTID;
        load('orb_fault_segmented','YC85SeismicMoment_char');
    else
        num_sources=num_comb; source_id=orb_faults_comb.Combined_ID;
        load('orb_fault_combined','YC85SeismicMoment_char','YC85SeismicMoment_exp');
    end
    
    MoRate_StochasticEventCatalog{ss}=zeros(num_sources,2); f_indx={}; 
    
    for ff=1:num_sources %loop through each fault
        
        %Theoritical moment rates
        if ss==1 || ss==2
            MoRate_StochasticEventCatalog{ss}(ff,1)=median(YC85SeismicMoment_char{ff});
        else
            MoRate_StochasticEventCatalog{ss}(ff,1)=median(YC85SeismicMoment_exp{ff});
        end
        
        %Moment rate for individual faults in catalog
        f_indx{ff}=find(StochasticEventCatalog{ss}(:,5) == source_id(ff)); %Index all events for fault
        Mo_StochasticEventCatalog = 10.^(1.5*StochasticEventCatalog{ss}(f_indx{ff},4)+9.05);
        MoRate_StochasticEventCatalog{ss}(ff,2) = sum(Mo_StochasticEventCatalog)/(num_simu*t_limit);
        clear f_indx  
      
    end %end ff loop
     clear YC85SeismicMoment_char num_source
    
end
    
%% Moment rate comparison plots %%

if model_opt==4; axislimits=[15.5 17]; else; axislimits=[14 16.5]; end
col_opt=['r','b']; 
gr=[0.2 0.5 0.2];%dark green colour for GR can't be called by letter

figure(1);

tiledlayout(1,2,"TileSpacing","compact"); nexttile

%Compare catalog moment rates to median input Y&C85 reccurence model moment rate
plot([axislimits(1) axislimits(2)],[axislimits(1) axislimits(2)],'k--'); hold on

for ss=1:length(StochasticEventCatalog)
    
    if ss==3
        plot(log10(MoRate_StochasticEventCatalog{ss}(:,1)),log10(MoRate_StochasticEventCatalog{ss}(:,2)),...
        'x','Color',gr,'MarkerSize',12,'LineWidth',1.5);     
    else
        plot(log10(MoRate_StochasticEventCatalog{ss}(:,1)),log10(MoRate_StochasticEventCatalog{ss}(:,2)),...
        'x','Color',col_opt(ss),'MarkerSize',12,'LineWidth',1.5);
    end
end

axis([axislimits(1) axislimits(2) axislimits(1)  axislimits(2)]); axis square; hold on;

xlabel('Y&C 85 Reccurrence Model Moment Rate (log Nm/yr)'),ylabel('Stochastic Event Catalog Moment Rate (log Nm/yr)'); hold on;
legend('y=x','segmented-char','combined-char','combined-GR','Location','southeast');
grid on; set(gca,'fontsize',11)

t1=title('(a)','fontsize',13,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile

%Compare catalog moment rates to theoritical moment rate
%(calc already incorporates 20% reduction of fault area) but not 10% fault reduction in slip rate

%moment rate from geodetic model slip rates
if model_opt==4
    load('model4/otago_geodetic_slip_rates');
    input_mo_rate=orb_faults.Area_m2.*orb_sliprates*3*10^10*10^-3;
    xlabel_txt='Geodetic Moment Rate (log Nm/yr)';
%moment rate from NZCFM slip rates
else
    input_mo_rate=orb_faults.SR_pref.*10^-3.*orb_faults.Area_m2*rigidity;
    xlabel_txt='NZCFM Moment Rate (log Nm/yr)';
end
    comb_mo_rate=c_area.*w_slip_rate_m.*rigidity;

plot([axislimits(1) axislimits(2)],[axislimits(1) axislimits(2)],'k--'); hold on

for ss=1:length(StochasticEventCatalog)
    
    if ss==1
        tmp_mo_rate=input_mo_rate;
    else
        tmp_mo_rate=comb_mo_rate;
    end
    
    if ss==3
         plot(log10(tmp_mo_rate),log10(MoRate_StochasticEventCatalog{ss}(:,2)),...
        'x','Color',gr,'MarkerSize',12,'LineWidth',1.5); hold on       
    else
        plot(log10(tmp_mo_rate),log10(MoRate_StochasticEventCatalog{ss}(:,2)),...
        'x','Color',col_opt(ss),'MarkerSize',12,'LineWidth',1.5); hold on
    end
end

axis([axislimits(1) axislimits(2) axislimits(1)  axislimits(2)]); axis square; hold on;
xlabel(xlabel_txt); ylabel('Stochastic Event Catalog Moment Rate (log Nm/yr)'); hold on;
legend('y=x','segmented-char','combined-char','combined-GR','Location','southeast');
grid on; set(gca,'fontsize',11);

t1=title('(b)','fontsize',13,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);


set(gcf,'position',[440 368 817 429]);

%% Plot target and actual mean recurrence intervals and standard deviation for Weibull distribution

if cat_option==3

    ss=2; %select catalog doing analysis for

    if ss==1
        num_sources=num_fault; source_id=orb_faults.OBJECTID;
        load('orb_fault_segmented','T','alpha_bpt'); indx_opt=1;
    elseif ss==2
        num_sources=num_comb; source_id=orb_faults_comb.Combined_ID;
        load('orb_fault_combined','T','alpha_bpt'); indx_opt=1;
    elseif ss==3
        num_sources=num_comb; source_id=orb_faults_comb.Combined_ID;
        load('orb_fault_combined','T','alpha_bpt'); indx_opt=2;
    end

    mean_RI=zeros(num_sources,1); std_RI=zeros(num_sources,1);
    median_T=zeros(num_sources,1); std_T=zeros(num_sources,1);
 
    for kk=1:num_sources
        
        flt_ruptures_weibull=find(StochasticEventCatalog_weibull{ss}(:,5)==source_id(kk));
        
        intertimes_wbl=zeros(length(flt_ruptures_weibull),1);
        
        for mm=1:length(flt_ruptures_weibull)-1
            intertimes_wbl(mm)=StochasticEventCatalog_weibull{ss}(flt_ruptures_weibull(mm+1),2)-StochasticEventCatalog_weibull{ss}(flt_ruptures_weibull(mm),2);
        end
        
        mean_RI(kk)=mean(intertimes_wbl); std_RI(kk)=std(intertimes_wbl);
        median_T(kk)=median(T{kk}(:,indx_opt)); alpha_tmp=median(alpha_bpt{kk}(:,indx_opt));
        std_T(kk)=median_T(kk)*alpha_tmp;
        
    end
    
    figure(101);
    
    plot(median_T,mean_RI,'rv','LineWidth',1);hold on; plot(std_T,std_RI,'b^','LineWidth',1);hold on
    plot([0 10000],[0 10000],'k--');hold on
    legend('RI mean','RI standard deviation','x=y','location','northwest');
    grid on; set(gca,'fontsize',11);
    xlabel('Target (years)'); ylabel('Weibull catalog record (Years)');
    xlim([0 10000]); ylim([0 10000]);
    axis square
end

%% Freq mag plot comparison for all earthquake records

fault_syncatRate=zeros(length(all_mag_range_GR),4);

for ll = 1:length(all_mag_range_GR)
     fault_syncatRate (ll,1) = length(find(StochasticEventCatalog{1}(:,4) >= all_mag_range_GR(ll)))/(num_simu*t_limit);%Rate from segmented-char catalog
     fault_syncatRate (ll,2) = length(find(StochasticEventCatalog{2}(:,4) >= all_mag_range_GR(ll)))/(num_simu*t_limit);%Rate from combined-GR catalog
     fault_syncatRate (ll,3) = length(find(StochasticEventCatalog{3}(:,4) >= all_mag_range_GR(ll)))/(num_simu*t_limit);%Rate from combined-GR catalog
end
  
%Make MFD Plot comparison

figure(3);

%Comparing all different fault-based catalogs
semilogy(all_mag_range_GR,fault_syncatRate(:,1),'r-','LineWidth',1.5); hold on; %MFD from segmented-char event catalog
semilogy(all_mag_range_GR,fault_syncatRate (:,2),'b-','LineWidth',1.5); hold on; %MFD from combined-char event catalog
semilogy(all_mag_range_GR,fault_syncatRate (:,3),'k-','LineWidth',1.5); hold on; %MFD from combined-GR event catalog
set(gca,'fontsize',12);  axis([Mmin Mmax 10^-5 2*10^-1]); 

if cat_option==1
    legend('segmented-char-Poisson','combined-char-Poisson','combined-GR-Poisson','Location','southwest');hold on;
elseif cat_option==2
    legend('segmented-char-BPT','combined-char-BPT','combined-GR-BPT','Location','southwest');hold on;
elseif cat_option==3
    legend('segmented-char-wbl','combined-char-wbl','combined-GR-wbl','Location','southwest');hold on;
end

xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;

%% Plot distribution of interoccurrence times

inter_occ_times=cell(3,1);

figure(4);

for ss=1:length(StochasticEventCatalog)
    
    inter_occ_times{ss}=zeros(length(StochasticEventCatalog{ss})-1,1);
    
    for kk=2:length(StochasticEventCatalog{ss})
        inter_occ_times{ss}(kk)=StochasticEventCatalog{ss}(kk,2)-StochasticEventCatalog{ss}(kk-1,2);
    end
   
    histogram(inter_occ_times{ss},'BinWidth',5,'facealpha',0.5); hold on

end

xlabel('interoccurrence time (yrs)'); xlim([0 50]);
legend(['segmented-char,' newline 'n = ' num2str(length(StochasticEventCatalog{1}))],...
        ['combined-char,' newline 'n = ' num2str(length(StochasticEventCatalog{2}))],...
        ['combined-gr,' newline 'n = ' num2str(length(StochasticEventCatalog{3}))]);


%% Sampling of stochastic catalogs

nsim=1000; %number of subsamples
empty_count=zeros(length(StochasticEventCatalog),1);

sampleAnnualRate=zeros(nsim,length(all_mag_range_GR),length(StochasticEventCatalog)); sampleMoRate=zeros(nsim,3);
clock_old = cputime; simu_count=1;
    
rand_t_indx=rand(nsim,1).*(num_simu-catalog_duration);%nsmim random numbers between 0 and num_simu-cat_duration
    
for gg=1:nsim
    
    if gg == 2000*simu_count
        clock_new = cputime;
        disp(['Simulation cycle is: ',num2str(gg),' & Required time is: ',num2str(clock_new-clock_old)]);
        clock_old = clock_new;
            simu_count = simu_count + 1;
    end 
    
    %loop through each catalog for each subsample
    for ss=1:length(StochasticEventCatalog)
    
        %index random 70 year catalog increment for random fault based catalog
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

%% For Poisson catalogs, run seperate subsampling for 10,000 periods

if model_opt==2 & cat_option == 1
sampleAnnualRate_10kyr=zeros(nsim,length(all_mag_range_GR),length(StochasticEventCatalog)); 

rand_t_indx=rand(nsim,1).*(num_simu-10000);%nsmim random numbers between 0 and num_simu-cat_duration    
    for gg=1:nsim
    
    %loop through each catalog for each subsample
    for ss=1:length(StochasticEventCatalog)
    
        %index random 10,000 year catalog increment for random fault based catalog
        tmp_catalog_indx=find(StochasticEventCatalog{ss}(:,2)>=rand_t_indx(gg) & StochasticEventCatalog{ss}(:,2)<(rand_t_indx(gg)+10000));
     
        if isempty(tmp_catalog_indx)==1
            empty_count(ss)=empty_count(ss)+1;
        end
    
        tmp_catalog=StochasticEventCatalog{ss}(tmp_catalog_indx,:);%sample catalog
    
        for ll=1:length(all_mag_range_GR)
            sampleAnnualRate_10kyr(gg,ll,ss) = length(find(tmp_catalog(:,4) >= all_mag_range_GR(ll)))/10000;
        end
     
    clear tmp_catalog_indx
    
    end %end ss loop for each catalog
end%end gg loop for each sample


save(strcat('model',num2str(model_opt),'/catalog_poi_10kyr_samples'),'sampleAnnualRate_10kyr')
end


%% Plot stochastic catalog MFD with sub samples

figure(5);
cat_col=vertcat([0.8 0.6 0.6 0.4],[0.6 0.6 0.8 0.4],[0.6 0.6 0.6 0.4]);
%Plot sub sampled catalogs
semilogy(all_mag_range_GR,squeeze(sampleAnnualRate(:,:,1)),'Color',cat_col(1,:),'LineWidth',0.5); hold on; 
semilogy(all_mag_range_GR,squeeze(sampleAnnualRate(:,:,2)),'Color',cat_col(2,:),'LineWidth',0.5); hold on; 
semilogy(all_mag_range_GR,squeeze(sampleAnnualRate(:,:,3)),'Color',cat_col(3,:),'LineWidth',0.5); hold on; 
p3=semilogy(all_mag_range_GR,fault_syncatRate (:,1),'r-','LineWidth',1.5); hold on; %MFD from segmented-char
p4=semilogy(all_mag_range_GR,fault_syncatRate(:,2),'b-','LineWidth',1.5); hold on; %MFD from combined-char
p5=semilogy(all_mag_range_GR,fault_syncatRate (:,3),'k-','LineWidth',1.5); hold on; %MFD from combined-GR

set(gca,'fontsize',12);  axis([Mmin Mmax 10^-5 5*10^-1]); legend([p3 p4 p5],'segmented-char','combined-char','combined-GR',...
    'Location','southwest');hold on;
xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;

m5_rate=zeros(3,2);
%Get average rate of M>5 events in the catalogs, plus std of this rate from catalog subsamples
for mm=1:3
    m5_rate(mm,:)=[length(StochasticEventCatalog{mm})/num_simu std(sampleAnnualRate(:,1,mm))];
end

%% Save analysis
    
if cat_option ==1 %Poisson
    variablename=strcat('model',num2str(model_opt),'/catalog_poi_statistics');
    sampleAnnualRate_poi=sampleAnnualRate; 
    sampleMoRate_poi=sampleMoRate;
    fault_syncatRate_poi=fault_syncatRate;
    MoRate_StochasticEventCatalog_poi=MoRate_StochasticEventCatalog;
    MoRate_StochasticEventCatalog_all_poi=MoRate_StochasticEventCatalog_all;
    
  
     save(variablename,'sampleAnnualRate_poi','sampleMoRate_poi','fault_syncatRate_poi',...
            'MoRate_StochasticEventCatalog_poi','MoRate_StochasticEventCatalog_all_poi','input_mo_rate');
    

elseif cat_option==2 %catalog is bpt
    variablename=strcat('model',num2str(model_opt),'/catalog_bpt_statistics');
    sampleAnnualRate_bpt=sampleAnnualRate; 
    sampleMoRate_bpt=sampleMoRate;
    fault_syncatRate_bpt=fault_syncatRate;
    MoRate_StochasticEventCatalog_bpt=MoRate_StochasticEventCatalog;
    MoRate_StochasticEventCatalog_all_bpt=MoRate_StochasticEventCatalog_all;
    
    save(variablename,'sampleAnnualRate_bpt','sampleMoRate_bpt','fault_syncatRate_bpt',...
        'MoRate_StochasticEventCatalog_bpt','MoRate_StochasticEventCatalog_all_bpt','input_mo_rate');

elseif cat_option==3 %catalog is weibull
    variablename=strcat('model',num2str(model_opt),'/catalog_wbl_statistics');
    sampleAnnualRate_wbl=sampleAnnualRate; 
    sampleMoRate_wbl=sampleMoRate;
    fault_syncatRate_wbl=fault_syncatRate;
    MoRate_StochasticEventCatalog_wbl=MoRate_StochasticEventCatalog;
    MoRate_StochasticEventCatalog_all_wbl=MoRate_StochasticEventCatalog_all;
    
    save(variablename,'sampleAnnualRate_wbl','sampleMoRate_wbl','fault_syncatRate_wbl',...
        'MoRate_StochasticEventCatalog_wbl','MoRate_StochasticEventCatalog_all_wbl','input_mo_rate');    
end

