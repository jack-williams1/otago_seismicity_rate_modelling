%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PLOT AND ANALYSE POISSON, BPT, AND WEIBULL  %%%
%%%    STOCHASTIC CATALOG FOR A SINGLE FAULT     %%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%select Poisson, BPT and Weibull catalogs

load catalog_poisson;
load catalog_bpt;
load catalog_weibull

% StochasticEventCatalog
% 1    ) Event number
% 2    ) Simulation cycle
% 3    ) Occurrence time
% 4    ) Earthquake magnitude
% 5    ) NCZFM ID
% 6    ) Recurrence type (1 = characteristic versus 2 = exponential)
% 7    ) Reccurrece parameter
% 8    ) SED (for BPT only)

%% MFD for particular fault or multifault

%Select fault by (1) NZCFM ID or (2) Combined ID
select_opt=1;

ff=184; %NZCFM ID for Akatore
%ff=30; %Combined ID for Dunstan

if select_opt==1
    f_indx=find(orb_faults.OBJECTID==ff); syncat_opt=1; indx_opt=0;
    load('orb_fault_segmented','mag_range_GR','fm_char_var1','fm_exp_var1',...
        'RecurrenceVar','T','lambda');
elseif select_opt==2
    f_indx=find(orb_faults_comb.Combined_ID==ff); syncat_opt=[1,2]; indx_opt=1;
     load('orb_fault_combined','mag_range_GR','fm_char_var1','fm_exp_var1',...
         'RecurrenceVar','T','lambda');
end

cat_rate_poi=zeros(length(mag_range_GR{f_indx}),length(syncat_opt));
cat_rate_bpt=zeros(length(mag_range_GR{f_indx}),length(syncat_opt));
cat_rate_weibull=zeros(length(mag_range_GR{f_indx}),length(syncat_opt));

for ss=1:length(syncat_opt)
    
    %index fault events in catalog
    flt_ruptures_poi{ss}=find(StochasticEventCatalog_poisson{ss+indx_opt}(:,5)==ff);
    flt_ruptures_bpt{ss}=find(StochasticEventCatalog_bpt{ss+indx_opt}(:,5)==ff);
    flt_ruptures_weibull{ss}=find(StochasticEventCatalog_weibull{ss+indx_opt}(:,5)==ff);

    for gg=1:length(mag_range_GR{f_indx})
        cat_rate_poi(gg,ss)=length(find(StochasticEventCatalog_poisson{ss+indx_opt}(flt_ruptures_poi{ss},4)>=mag_range_GR{f_indx}(gg)))/(t_limit*num_simu);
        cat_rate_bpt(gg,ss)=length(find(StochasticEventCatalog_bpt{ss+indx_opt}(flt_ruptures_bpt{ss},4)>=mag_range_GR{f_indx}(gg)))/(t_limit*num_simu);
        cat_rate_weibull(gg,ss)=length(find(StochasticEventCatalog_weibull{ss+indx_opt}(flt_ruptures_weibull{ss},4)>=mag_range_GR{f_indx}(gg)))/(t_limit*num_simu);
    end
    
end

col_opt=["k","r","m"];%set colours for plotting poisson, bpt, or weibull results

%% Make MFD plots

%Plot MFD from catalogs only

figure(200);

legendfntsize=9; fntsize=10; ylimits=[10^-5 10^-1];%these parameters may need to change depending on fault
marker_style=["*","v"]; 

for ss=1:length(syncat_opt)
    
    marker_indicies=[1:30:length(cat_rate_poi(:,ss))];%plot markers at increments of 30 indicies
    
    semilogy(mag_range_GR{f_indx},cat_rate_poi(:,ss),col_opt(1),'LineWidth',1.2,'Marker',marker_style(ss),'MarkerIndices',marker_indicies); hold on;
    semilogy(mag_range_GR{f_indx},cat_rate_bpt(:,ss),col_opt(2),'LineWidth',1.2,'Marker',marker_style(ss),'MarkerIndices',abs(marker_indicies-10)); hold on;
    semilogy(mag_range_GR{f_indx},cat_rate_weibull(:,ss),col_opt(3),'LineWidth',1.2,'Marker',marker_style(ss),'MarkerIndices',abs(marker_indicies-20)); hold on;
end

if select_opt==1
    legend({'poi-char','bpt-char','wbl-char'},...
        'FontSize',legendfntsize,'Location','southwest'); hold on;
else
    legend({'poi-char','bpt-char','wbl-char','poi-GR','bpt-GR','wbl-GR'},...
        'FontSize',legendfntsize,'Location','southwest'); hold on;    
end

set(gca,'fontsize',fntsize); axis([Mmin Mmax ylimits(1) ylimits(2)]); hold on;
    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;

%Plot MFD from catalogs and YC85 Recurrence model

figure(201);

for ss=1:length(syncat_opt)
    
    marker_indicies=[1:30:length(cat_rate_poi(:,ss))];%plot markers at increments of 30 indicies
    
    semilogy(mag_range_GR{f_indx},cat_rate_poi(:,ss),col_opt(1),'LineWidth',1.2,'LineStyle','--','Marker',marker_style(ss),'MarkerIndices',marker_indicies); hold on
    semilogy(mag_range_GR{f_indx},cat_rate_bpt(:,ss),col_opt(2),'LineWidth',1.2,'LineStyle','--','Marker',marker_style(ss),'MarkerIndices',marker_indicies+10); hold on
    semilogy(mag_range_GR{f_indx},cat_rate_weibull(:,ss),col_opt(3),'LineWidth',1.2,'LineStyle','--','Marker',marker_style(ss),'MarkerIndices',marker_indicies+20); hold on
end

%select median (i.e. intermediate) recurrence curve
median_ri=(length(RecurrenceVar{f_indx})/2)+0.5;

if select_opt==1
    semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,median_ri*3),'b-','LineWidth',1.5); hold on;
    legend({'poi-char','bpt-char','weibull-char','YC85-char'},...
      'FontSize',legendfntsize,'Location','southwest'); hold on;
else
    semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,median_ri*3),'b-','LineWidth',1.5); hold on;
    semilogy(mag_range_GR{f_indx},fm_exp_var1{f_indx}(:,median_ri*3),'Color',[0.2 0.5 0.2],'LineStyle','-','LineWidth',1.5); hold on;
    legend({'poi-char','bpt-char','weibull-char','poi-GR','bpt-GR','weibull-GR','YC85-char','YC85-GR'},...
        'FontSize',legendfntsize,'Location','southwest'); hold on;
end

set(gca,'fontsize',fntsize); axis([Mmin Mmax-0.2 ylimits(1) ylimits(2)]); hold on;
    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;
   
%Plot all YC85 Recurrence model recurrence models
figure(202);

for ss=1:length(RecurrenceVar{f_indx})
    
    tmp_width=(RecurrenceVar{f_indx}(ss,4))^0.5*5;
    tmp_weight=(RecurrenceVar{f_indx}(ss,4))^0.3;
    
    semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,3*(ss-1)+3),'Color',[0 0 1 tmp_weight],'LineWidth',tmp_width); hold on;
    
    if select_opt==2
        semilogy(mag_range_GR{f_indx},fm_exp_var1{f_indx}(:,3*(ss-1)+3),'Color',[0.2 0.5 0.2 tmp_weight],'LineWidth',tmp_width); hold on;
    end
  
end

p1=semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,median_ri*3),'Color',[0 0 1 0.5]); hold on;

if select_opt==1 
     legend(p1,'YC85-char','FontSize',legendfntsize,'Location','northeast'); hold on;
else
    p2=semilogy(mag_range_GR{f_indx},fm_exp_var1{f_indx}(:,median_ri*3),'Color',[0.2 0.5 0.5]); hold on;
    legend([p1,p2],{'YC85-char','YC85-GR'},'FontSize',legendfntsize,'Location','northeast'); hold on;    
end

set(gca,'fontsize',fntsize); axis([Mmin Mmax ylimits(1) ylimits(2)]); hold on;
    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;

%% Empirical and theoritical ecdf's of intervent times and hazard functions

figure(203);
tiledlayout(length(syncat_opt),2,"TileSpacing","compact");

tmp_limit=[1500,150]; bylim=[0.01 0.07]; %set based on indivdual fault

for ss=1:length(syncat_opt)
    
    intertimes_bpt{ss}=zeros(length(flt_ruptures_bpt{ss}),1); intertimes_poi{ss}=zeros(length(flt_ruptures_poi{ss}),1);
    intertimes_wbl{ss}=zeros(length(flt_ruptures_weibull{ss}),1);
    
    %obtain empiricial interevent times from stochastic catalogs
    
    for kk=1:length(flt_ruptures_bpt{ss})-1   
        intertimes_bpt{ss}(kk)=StochasticEventCatalog_bpt{ss+indx_opt}(flt_ruptures_bpt{ss}(kk+1),2)-StochasticEventCatalog_bpt{ss+indx_opt}(flt_ruptures_bpt{ss}(kk),2);   
    end

    for kk=1:length(flt_ruptures_poi{ss})-1   
        intertimes_poi{ss}(kk)=StochasticEventCatalog_poisson{ss+indx_opt}(flt_ruptures_poi{ss}(kk+1),2)-StochasticEventCatalog_poisson{ss+indx_opt}(flt_ruptures_poi{ss}(kk),2);   
    end    

    for kk=1:length(flt_ruptures_weibull{ss})-1   
        intertimes_wbl{ss}(kk)=StochasticEventCatalog_weibull{ss+indx_opt}(flt_ruptures_weibull{ss}(kk+1),2)-StochasticEventCatalog_weibull{ss+indx_opt}(flt_ruptures_weibull{ss}(kk),2);   
    end 
    
        
    intertimes={}; 
    intertimes={intertimes_poi{ss},intertimes_bpt{ss},intertimes_wbl{ss}};
    
    %obtain interevent time distributions from input probability distributions using median recurrence model
    
    exp_cdf=cdf('Exponential',1:tmp_limit(ss),median(T{f_indx}(:,ss)));
    exp_pdf=pdf('Exponential',1:tmp_limit(ss),median(T{f_indx}(:,ss)));
    
    %cdf of interevent times follows Inverse Gaussian Distribution for BPT, see Eq. 12 of Matthews et al (2002)
    bpt_cdf=cdf('InverseGaussian',1:tmp_limit(ss),median(T{f_indx}(:,ss)),median(lambda{f_indx}(:,ss)));
    bpt_pdf=pdf('InverseGaussian',1:tmp_limit(ss),median(T{f_indx}(:,ss)),median(lambda{f_indx}(:,ss)));    
    
    wbl_cdf=cdf('Weibull',1:tmp_limit(ss),tau_wbl{ss+indx_opt}(f_indx),beta_wbl{ss+indx_opt}(f_indx));
    wbl_pdf=pdf('Weibull',1:tmp_limit(ss),tau_wbl{ss+indx_opt}(f_indx),beta_wbl{ss+indx_opt}(f_indx));
    
    nexttile %first column, plot empirical and theorictical interevent time distributions
    
    plot(1:tmp_limit(ss),exp_cdf,col_opt(1),'Linewidth',1.2); hold on  
    plot(1:tmp_limit(ss),bpt_cdf,col_opt(2),'Linewidth',1.2); hold on
    plot(1:tmp_limit(ss),wbl_cdf,col_opt(3),'Linewidth',1.2); hold on 
    
    for dd=1:3
        [f,x]=ecdf(intertimes{dd});plot(x,f,'Color',col_opt(dd),'Linewidth',1.2,'LineStyle','--'); hold on
    end
    
    xlim([0 tmp_limit(ss)]); axis square
    xlabel ('interevent time (yrs)');
    
    if ss==2 
        legend('GR-target poi','GR-target bpt','GR-target wbl',...
            'GR-catalog poi','GR-catalog bpt','GR-catalog wbl','Location','Southeast');
    else
         legend('char-target poi','char-target bpt','char-target wbl',...
            'char-catalog poi','char-catalog bpt','char-catalog wbl','Location','Southeast');
    end    
        
    nexttile %second column, plot empirical and theoritical hazard functions
    
    h_exp=exp_pdf./(1-exp_cdf); %theoritical hazard function following generic eq. 3 of Yakovlev et al (2006) 
    h_bpt=bpt_pdf./(1-bpt_cdf); 
    h_wbl=wbl_pdf./(1-wbl_cdf); 
    
    %plot theoritical hazard functions
    plot(1:tmp_limit(ss),h_exp,col_opt(1),'Linewidth',1.4);hold on
    plot(1:tmp_limit(ss),h_bpt,col_opt(2),'Linewidth',1.4);hold on
    plot(1:tmp_limit(ss),h_wbl,col_opt(3),'Linewidth',1.4);hold on
    
    %derive cumulative hazard function (i.e. integral of hazard function)

    bin=[50,10]; %time period over which averaging hazard rate    
    edges=[1:bin(ss):tmp_limit(ss)+bin(ss)];
    
    
    %calculate time average hazard rate from cumulative hazard function over
    %time period of bin following: https://www.itl.nist.gov/div898/handbook/apr/section1/apr123.htm
     
    for dd=1:3
            [f,x]=ecdf(intertimes{dd},'function','cumhazard'); 
            afr=zeros(length(edges)-1,1);
            
            for ii=1:(length(edges)-1)
                %find closest index to empirical cumlualtive hazard(t) and t
                [~,tmpindx1]=min(abs(edges(ii)-x));
                [~,tmpindx2]=min(abs(edges(ii+1)-x));
                afr(ii)=(f(tmpindx2)-f(tmpindx1))/bin(ss); %time average hazard rate
            end  
            
            plot(edges(1:end-1),afr,'Color',col_opt(dd),'Linewidth',1.4,'LineStyle','--'); hold on
    end

    xlim([0 tmp_limit(ss)]); axis square; ylim([0 bylim(ss)])
    xlabel ('time since last event (yrs)'); ylabel('hazard rate');    

    if ss==2 
        legend('GR-target poi','GR-target bpt','GR-target wbl',...
            'GR-catalog poi','GR-catalog bpt','GR-catalog wbl','Location','Northeast');
    else
         legend('char-target poi','char-target bpt','char-target wbl',...
            'char-catalog poi','char-catalog bpt','char-catalog wbl','Location','Northeast');
    end  
    
    
end

if select_opt==2
    set(gcf,'Position',[680 425 551 552]);
else
    set(gcf,'Position',[562 457 691 364]);
end
