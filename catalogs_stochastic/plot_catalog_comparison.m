%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CATALOG AND MO RATE COMPARISON PLOT FOR PAPER (FIGURE S7)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%Select recurrence model for Otago Faults to assess:
%1) NZCFM slip rates - Aperiodicity = 0.5
%2) NZCFM slip rates - Aperiodicity = 2
%3) NZCFM slip rates - Aperiodicity = 4
%4) Geodetic slip rates - Aperiodicity = 2


marker_style=["^","v","*"]; line_cols=["k","r","m"];%set colours for plotting poisson, bpt, or weibull results
leg_col=vertcat([0 0 0],[1 0 0],[1 0 1]);

title_opt=["(a)","(b)"]; fntsize=7;

model_opt=[1,2];

catalog_opt=["poisson","bpt","weibull"];

tiledlayout(3,2,"TileSpacing","compact");

for mm=1:length(model_opt)
    
    nexttile

    model_path=strcat('model',num2str(model_opt(mm)));
    addpath(model_path)

    load catalog_poisson; load('catalog_bpt','StochasticEventCatalog_bpt'); load catalog_weibull

    load orb_fault_parameters
    load('orb_fault_segmented','num_fault','orb_faults');
    load('orb_fault_combined','num_comb','orb_faults_comb','c_area','w_slip_rate_m');
    
    all_mag_range_GR=Mmin:dm:Mmax;
    marker_indicies=[1:30:length(all_mag_range_GR)];%plot markers at increments of 30 indicies
    
    for ss=1:3
            
        StochasticEventCatalog=eval(strcat('StochasticEventCatalog_',catalog_opt(ss)));%loop through catalog for each recurrence model
        fault_syncatRate=zeros(length(all_mag_range_GR),4);

        %find event rate for each model
        for ll = 1:length(all_mag_range_GR)
            fault_syncatRate (ll,1) = length(find(StochasticEventCatalog{1}(:,4) >= all_mag_range_GR(ll)))/(num_simu*t_limit);%Rate from segmented-char catalog
            fault_syncatRate (ll,2) = length(find(StochasticEventCatalog{2}(:,4) >= all_mag_range_GR(ll)))/(num_simu*t_limit);%Rate from combined-GR catalog
            fault_syncatRate (ll,3) = length(find(StochasticEventCatalog{3}(:,4) >= all_mag_range_GR(ll)))/(num_simu*t_limit);%Rate from combined-GR catalog
        end
        
    
    mfd_seg_char(ss)=semilogy(all_mag_range_GR,fault_syncatRate(:,1),'Color',line_cols(1),'LineWidth',1.2,...
        'Marker',marker_style(ss),'MarkerIndices',marker_indicies); hold on; 
    mfd_comb_char(ss)=semilogy(all_mag_range_GR,fault_syncatRate(:,2),'Color',line_cols(2),'LineWidth',1.2,...
        'Marker',marker_style(ss),'MarkerIndices',abs(marker_indicies-10)); hold on;
    mfd_comb_gr(ss)=semilogy(all_mag_range_GR,fault_syncatRate(:,3),'Color',line_cols(3),'LineWidth',1.2,...
        'Marker',marker_style(ss),'MarkerIndices',abs(marker_indicies-20)); hold on; 
    end
    

    [legh,objh] =legend([line([1,1],[1,1],'LineWidth',1.2),...
    mfd_seg_char(ss),mfd_comb_char(ss),mfd_comb_gr(ss),mfd_seg_char(1),mfd_seg_char(2),mfd_seg_char(3)],...
    {'segmented-char','combined-char','combined-gr','poisson','BPT','weibull',},...
    'FontSize',fntsize,'Location','southwest');

    %Shrink markers so not visible on legend    
    lineh1 = findobj(objh(7:8),'type','line'); set(lineh1,'MarkerSize',0.1); set(lineh1,'Color',leg_col(1,:));
    lineh2 = findobj(objh(9:10),'type','line'); set(lineh2,'MarkerSize',0.1); set(lineh2,'Color',leg_col(2,:));
    lineh3 = findobj(objh(11:12),'type','line'); set(lineh3,'MarkerSize',0.1); set(lineh3,'Color',leg_col(3,:));

    %Change some line colors to grey
    linei = findobj(objh(13:18),'type','line'); set(linei,'Color',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);

    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;
    xlim([Mmin Mmax]); ylim([10^-5 5*10^-2]); set(gca,'fontsize',fntsize)
    
    t2=title(title_opt(mm),'fontsize',fntsize+2,'fontweight','normal');
    set(t2, 'horizontalAlignment', 'left');
    set(t2, 'units', 'normalized');
    h2 = get(t2, 'position');
    set(t2, 'position', [-0.15 h2(2) h2(3)]);


    rmpath(model_path)
end

nexttile

% Moment rate comparison plots %%
addpath(model_path)

axislimits=[14 16.5];
marker_cols=['k','r','m']; 
%gr=[0.2 0.5 0.2];%dark green colour for GR can't be called by letter

%Compare catalog moment rates to median input Y&C85 reccurence model moment rate
plot([axislimits(1) axislimits(2)],[axislimits(1) axislimits(2)],'k--'); hold on

for ss=1:length(StochasticEventCatalog)
    
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


     plot(log10(MoRate_StochasticEventCatalog{ss}(:,1)),log10(MoRate_StochasticEventCatalog{ss}(:,2)),...
        'x','Color',marker_cols(ss),'MarkerSize',12,'LineWidth',1.5);     

end

axis([axislimits(1) axislimits(2) axislimits(1)  axislimits(2)]); axis square; hold on;

xlabel('Y&C 85 Reccurrence Model Moment Rate (log Nm/yr)'),ylabel(['Stochastic Event Catalog',char(10),'Moment Rate (log Nm/yr)']); hold on;
legend('y=x','segmented-char','combined-char','combined-GR','Location','southeast');
grid on; set(gca,'fontsize',fntsize)

t1=title('(c)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile

%Compare catalog moment rates to theoritical moment rate
%(calc already incorporates 20% reduction of fault area) but not 10% fault reduction in slip rate

input_mo_rate=orb_faults.SR_pref.*10^-3.*orb_faults.Area_m2*rigidity;

comb_mo_rate=c_area.*w_slip_rate_m.*rigidity;

plot([axislimits(1) axislimits(2)],[axislimits(1) axislimits(2)],'k--'); hold on

for ss=1:length(StochasticEventCatalog)
    
    if ss==1
        tmp_mo_rate=input_mo_rate;
    else
        tmp_mo_rate=comb_mo_rate;
    end 

   plot(log10(tmp_mo_rate),log10(MoRate_StochasticEventCatalog{ss}(:,2)),...
        'x','Color',marker_cols(ss),'MarkerSize',12,'LineWidth',1.5); hold on

end

axis([axislimits(1) axislimits(2) axislimits(1)  axislimits(2)]); axis square; hold on;
xlabel('NZCFM Moment Rate (log Nm/yr)'); ylabel(['Stochastic Event Catalog', char(10), 'Moment Rate (log Nm/yr)']); hold on;
legend('y=x','segmented-char','combined-char','combined-GR','Location','southeast');
grid on; set(gca,'fontsize',fntsize);

t1=title('(d)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile

%Plot MFD for single fault against its median recurrence model

%Select fault by (1) NZCFM ID or (2) Combined ID
select_opt=2;

%ff=184; %NZCFM ID for Akatore
ff=30; %Combined ID for Dunstan

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
    
     marker_indicies=[1:30:length(cat_rate_poi(:,ss))];
    
    %index fault events in catalog
    flt_ruptures_poi{ss}=find(StochasticEventCatalog_poisson{ss+indx_opt}(:,5)==ff);
    flt_ruptures_bpt{ss}=find(StochasticEventCatalog_bpt{ss+indx_opt}(:,5)==ff);
    flt_ruptures_weibull{ss}=find(StochasticEventCatalog_weibull{ss+indx_opt}(:,5)==ff);

    for gg=1:length(mag_range_GR{f_indx})
        cat_rate_poi(gg,ss)=length(find(StochasticEventCatalog_poisson{ss+indx_opt}(flt_ruptures_poi{ss},4)>=mag_range_GR{f_indx}(gg)))/(t_limit*num_simu);
        cat_rate_bpt(gg,ss)=length(find(StochasticEventCatalog_bpt{ss+indx_opt}(flt_ruptures_bpt{ss},4)>=mag_range_GR{f_indx}(gg)))/(t_limit*num_simu);
        cat_rate_weibull(gg,ss)=length(find(StochasticEventCatalog_weibull{ss+indx_opt}(flt_ruptures_weibull{ss},4)>=mag_range_GR{f_indx}(gg)))/(t_limit*num_simu);
    end
    
    semilogy(mag_range_GR{f_indx},cat_rate_poi(:,ss),line_cols(ss+indx_opt),'LineWidth',1.2,'LineStyle','-','Marker',marker_style(1),'MarkerIndices',marker_indicies); hold on
    semilogy(mag_range_GR{f_indx},cat_rate_bpt(:,ss),line_cols(ss+indx_opt),'LineWidth',1.2,'LineStyle','-','Marker',marker_style(2),'MarkerIndices',marker_indicies+10); hold on
    semilogy(mag_range_GR{f_indx},cat_rate_weibull(:,ss),line_cols(ss+indx_opt),'LineWidth',1.2,'LineStyle','-','Marker',marker_style(3),'MarkerIndices',marker_indicies+20); hold on
    
    
end

%select median (i.e. intermediate) recurrence curve
median_ri=(length(RecurrenceVar{f_indx})/2)+0.5;

if select_opt==1
    semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,median_ri*3),'b-','LineWidth',1.7); hold on;
    legend({'poi-char','bpt-char','weibull-char','YC85-char'},...
      'FontSize',legendfntsize,'Location','northeast'); hold on;
else
    semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,median_ri*3),'b-','LineWidth',1.7); hold on;
    semilogy(mag_range_GR{f_indx},fm_exp_var1{f_indx}(:,median_ri*3),'Color',[0.3 0.3 0.3],'LineStyle','-','LineWidth',1.7); hold on;
    legend({'poi-char','bpt-char','weibull-char','poi-GR','bpt-GR','weibull-GR','YC85-char','YC85-GR'},...
        'Location','northeast'); hold on;
end

set(gca,'fontsize',fntsize); axis([Mmin Mmax-0.2 10^-5 10^-1]); hold on;
    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;
    
t1=title('(e)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);
    
nexttile

% Plot target and actual mean recurrence intervals and standard deviation for Weibull distribution

ss=1 ; %select catalog doing analysis for

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
    
plot(median_T,mean_RI,'rv','LineWidth',1);hold on; plot(std_T,std_RI,'b^','LineWidth',1);hold on
plot([0 10000],[0 10000],'k--');hold on
legend('RI mean','RI standard deviation','y=x','location','northwest');
grid on; set(gca,'fontsize',fntsize);
xlabel('Target (years)'); ylabel('Weibull catalog record (Years)');
xlim([0 10000]); ylim([0 10000]); axis square

t1=title('(f)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);


set(gcf,'position',[603 86 561 711]);