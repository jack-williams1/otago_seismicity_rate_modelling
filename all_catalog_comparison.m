%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   COMPARE ALL CATALOGS: INSTRUMENTAL, STOCHASTICS, RSQSIM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

addpath('catalogs_instrumental'); addpath('catalogs_stochastic'); addpath('paleoseismic_constraints')
addpath('catalogs_rsqsim'); addpath('nshm_inversion/mfd_analysis'); addpath('nshm_urz'); addpath('nshm_inversion/solvis/work')

%Select stochastic earthquake catalog recurrence model for Otago Faults to assess:
%1) NZ CFM slip rates - Aperiodicity = 0.5
%2) NZ CFM slip rates - Aperiodicity = 2
%3) NZ CFM slip rates - Aperiodicity = 4
%4) Geodetic slip rates - Aperiodicity = 2

model_opt=2;
model_path=strcat('catalogs_stochastic/model',num2str(model_opt)); addpath(model_path)

load catalog_poisson; load('catalog_bpt','StochasticEventCatalog_bpt'); load catalog_weibull
load catalog_poi_statistics; load catalog_bpt_statistics; load catalog_wbl_statistics
sampleMoRate_tmp={sampleMoRate_poi,sampleMoRate_bpt,sampleMoRate_wbl};  
load orb_catalog; load catalog_rsqsim_statistics
load('orb_fault_parameters','Mmin','Mmax','dm'); load('orb_fault_segmented','orb_faults');
load('nshm_otago_inversion_results'); load nshm_urz_analysis; load otago_paleoseismic_constraint;
load sampleAnnualRate_50; 

runs1=0; runs2=0; runs3=0; %so only resets empty sample finder once
fntsize=9;  all_mag_range_GR=Mmin:dm:Mmax; 

if model_opt==4
    tmp_label_opt='Geodetic Model';
    load('model4/otago_geodetic_slip_rates');
    input_mo_rate=orb_faults.Area_m2.*orb_sliprates.*3*10^10*10^-3;
else
    tmp_label_opt='NZ CFM';
    input_mo_rate=orb_faults.SR_pref.*10^-3.*orb_faults.Area_m2*3*10^10;
end

b=1.5; c=9.05;

% MFD for instrumental catalog

ins_mag_range = [ins_Mmin:0.05:ins_Mmax]; ins_AnnualRate=zeros(length(ins_mag_range),1);

for mm=1:length(ins_mag_range)
    ins_AnnualRate(mm) = length(find(otago_aug_catalog_mat2(:,8) >= ins_mag_range(mm)))/(catalog_duration); %annual event rate 
end

ins_mo_rate=sum(10.^(1.5*otago_aug_catalog_mat2(:,8)+9.05))/catalog_duration;

%% Compare stochastic catalogs MFD shape and catalog subsamples

figure(3);

tiledlayout(3,4,"TileSpacing","compact");


legend_text=vertcat(["segmented-char-poi","segmented-char-BPT","segmented-char-wbl"],...
                    ["combined-char-poi","combined-char-BPT","combined-char-wbl"],...
                    ["combined-GR-poi","combined-GR-BPT","combined-GR-wbl"]);
                                        
line_cols=['k','r','m']; cat_col=vertcat([0.6 0.6 0.6 1.0]);

%kk loop is for on-fault MFD. %mm loop is for interevent distribution
for kk=1:3
    
    if runs1==0
        %fix samples with 0 moment rate to 1 so can plot in pdf
         for mm=1:3
            empty_sample=find(sampleMoRate_tmp{mm}(:,kk)==0); 
            empty_count(kk,mm)=length(empty_sample); sampleMoRate_tmp{mm}(empty_sample,kk)=1;
         end
    end
    
    nexttile
    
    for mm=1:3 %for each of the interevent time distributions
        
        if mm==1
        
         p1=semilogy(all_mag_range_GR,squeeze(sampleAnnualRate_poi(:,:,kk)),'Color',cat_col,'LineWidth',0.5); hold on; %MFD from poi event catalog sample
         p2=semilogy(all_mag_range_GR,fault_syncatRate_poi(:,kk),'k-','LineWidth',1.5); hold on; %MFD from poi catalog
         %p3=semilogy(ins_mag_range,ins_AnnualRate,'b--','LineWidth',1.5); hold on; %Uncomment to also plot MFD from event catalog
         
        elseif mm ==2
        
         p1=semilogy(all_mag_range_GR,squeeze(sampleAnnualRate_bpt(:,:,kk)),'Color',cat_col,'LineWidth',0.5); hold on; %MFD from poi catalog
         p2=semilogy(all_mag_range_GR,fault_syncatRate_bpt(:,kk),'k-','LineWidth',1.5); hold on; %MFD from BPT event catalog
         %p3=semilogy(ins_mag_range,ins_AnnualRate,'b--','LineWidth',1.5); hold on; %Uncomment to also plot MFD from event catalog
        elseif mm==3
            
         p1=semilogy(all_mag_range_GR,squeeze(sampleAnnualRate_wbl(:,:,kk)),'Color',cat_col,'LineWidth',0.5); hold on; %MFD from poi catalog
         p2=semilogy(all_mag_range_GR,fault_syncatRate_wbl(:,kk),'k-','LineWidth',1.5); hold on; %MFD from BPT event catalog  
         %p3=semilogy(ins_mag_range,ins_AnnualRate,'b--','LineWidth',1.5); hold on; %Uncomment to also plot MFD from event catalog
         
        end
        
        set(gca,'fontsize',fntsize);  axis([Mmin Mmax 10^-5 10^0]);
        
        legend_txt=[{strjoin({char(legend_text(kk,mm)), ['empty samples = ', num2str(empty_count(kk,mm))]},'\n')}];
        legend([p2],legend_txt,'Location','southwest','fontsize',fntsize-1);hold on;
        
        xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;
   
        nexttile
    end
   
    %Plot kernel distribution for moment rate
    
    for mm=1:3
        
        pd_i = fitdist(log10(sampleMoRate_tmp{mm}(:,kk)),'kernel');
        x_i = 0:0.05:max(log10(sampleMoRate_tmp{mm}(:,kk)));
        y_i = pdf(pd_i,x_i);
        
        plot(x_i,y_i,'Color',line_cols(mm),'linewidth',1.3);hold on
    end
    
    plot([log10(ins_mo_rate) log10(ins_mo_rate)],[0 1],'Color', [1 0.6 0],'linestyle','--','linewidth',1.3);
    plot([log10(sum(input_mo_rate)) log10(sum(input_mo_rate))],[0 1],'b--','linewidth',1.3);
    legend([legend_text(kk,:),'Instrumental',tmp_label_opt],'Location','northwest','fontsize',fntsize-1);hold on;

    set(gca,'fontsize',fntsize-1); axis square; xlabel('log10 (moment rate)'); ylim([0 1]); xlim([1 19]);
   
end

set(gcf,'Position',[680 130 965 847]);
runs1=1;


%% Compare all stochastic catalogs MFD (for a given aperiodicity)

figure(4);

line_cols=['k','r','m']; leg_col=vertcat([0 0 0],[1 0 0],[1 0 1]);
marker_style=["^","v","*"]; 
                
for kk=1:3
    
    marker_indicies=[1:30:length(all_mag_range_GR)];%plot markers at increments of 30 indicies
    
    mfd_poi(kk)=semilogy(all_mag_range_GR,fault_syncatRate_poi(:,kk),'Color',line_cols(1),'LineWidth',1.2,...
        'Marker',marker_style(kk),'MarkerIndices',marker_indicies); hold on; %MFD from poi catalog
    mfd_bpt(kk)=semilogy(all_mag_range_GR,fault_syncatRate_bpt(:,kk),'Color',line_cols(2),'LineWidth',1.2,...
        'Marker',marker_style(kk),'MarkerIndices',abs(marker_indicies-10)); hold on; %MFD from bpt catalog
    mfd_wbl(kk)=semilogy(all_mag_range_GR,fault_syncatRate_wbl(:,kk),'Color',line_cols(3),'LineWidth',1.2,...
        'Marker',marker_style(kk),'MarkerIndices',abs(marker_indicies-20)); hold on; %MFD from wbl catalog
    
end


[legh,objh] =legend([line([1,1],[1,1],'LineWidth',1.2),...
    mfd_poi(kk),mfd_bpt(kk),mfd_wbl(kk),mfd_poi(1),mfd_poi(2),mfd_poi(3)],...
    {'poisson','BPT','weibull','combined-gr','segmented-char','combined-char'},...
    'FontSize',9,'Location','northeast');

%Shrink markers so not visible on legend    
lineh1 = findobj(objh(7:8),'type','line'); set(lineh1,'MarkerSize',0.1); set(lineh1,'Color',leg_col(1,:));
lineh2 = findobj(objh(9:10),'type','line'); set(lineh2,'MarkerSize',0.1); set(lineh2,'Color',leg_col(2,:));
lineh3 = findobj(objh(11:12),'type','line'); set(lineh3,'MarkerSize',0.1); set(lineh3,'Color',leg_col(3,:));

%Change some line colors to grey
linei = findobj(objh(13:18),'type','line'); set(linei,'Color',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);

xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;
xlim([Mmin Mmax]); ylim([10^-5 5*10^-2]); 


%% Compare RSQSim and Stochastic Event Catalogs Interoccurrence times

figure(5);

tiledlayout(2,2,"TileSpacing","compact");

txt_opt=["(b)","(c)","(d)"]; binw=5; ymax=60000; xmax=100;

otago_rsqsimcatalog=readtable('catalogs_rsqsim/otago_rsqsim_catalog/eqs_otago_1e6yr_partial_mw.csv');

%only use events on Otago faults in interoccurrence times!
otago_indx=find(otago_rsqsimcatalog.num_otago_fault>0); 
rsqsim_inter_occ_times=zeros(length(otago_indx)-1,1);

for rr=1:length(rsqsim_inter_occ_times)
     %Interoccurrence time for Otago Fault events
     rsqsim_inter_occ_times(rr)=(otago_rsqsimcatalog.t0(otago_indx(rr+1))-otago_rsqsimcatalog.t0(otago_indx(rr)))/(3600*365*24);

end


nexttile

%Plot histogram for RSQSim occurrence times
histogram(rsqsim_inter_occ_times,'BinWidth',binw,'facealpha',0.5,'facecolor',[1 0 0]); hold on
xlabel('interoccurrence time (yrs)'); xlim([0 xmax]);
ylabel('count'); ylim([0 ymax]); legend(['rsqsim' newline 'n = ' num2str(length(rsqsim_inter_occ_times))]); axis square

t4=title('(a)','Fontsize',12,'fontweight','normal');
set(t4, 'horizontalAlignment', 'left');
set(t4, 'units', 'normalized');
h1 = get(t4, 'position');
set(t4, 'position', [-0.2 h1(2) h1(3)]);
    
for ss=1:3 %loop through catalogs by interevent time distribution
   
    nexttile
    
    inter_occ_times=cell(3,1);
    
    if ss==1; tmp_stochastic_catalog=StochasticEventCatalog_poisson; 
    elseif ss==2; tmp_stochastic_catalog=StochasticEventCatalog_bpt; 
    elseif ss==3; tmp_stochastic_catalog=StochasticEventCatalog_weibull; end

    for tt=1:length(tmp_stochastic_catalog) %loop through catalog by MFD
        inter_occ_times{tt}=zeros(length(tmp_stochastic_catalog{tt})-1,1);
        
        for kk=2:length(tmp_stochastic_catalog{tt})
            inter_occ_times{tt}(kk)=tmp_stochastic_catalog{tt}(kk,2)-tmp_stochastic_catalog{tt}(kk-1,2);
        end
   
    histogram(inter_occ_times{tt},'BinWidth',binw,'facealpha',0.5); hold on

    end
    
    xlabel('interoccurrence time (yrs)'); xlim([0 xmax]);
    ylabel('count'); ylim([0 ymax]); axis square
    legend(['segmented-char,' newline 'n = ' num2str(length(tmp_stochastic_catalog{1}))],...
        ['combined-char,' newline 'n = ' num2str(length(tmp_stochastic_catalog{2}))],...
        ['combined-gr,' newline 'n = ' num2str(length(tmp_stochastic_catalog{3}))],'fontsize',7);
    
    t4=title(txt_opt(ss),'Fontsize',12,'fontweight','normal');
    set(t4, 'horizontalAlignment', 'left');
    set(t4, 'units', 'normalized');
    h1 = get(t4, 'position');
    set(t4, 'position', [-0.2 h1(2) h1(3)]);
    
end

set(gcf,'Position',[680 301 652 676]);

clear otago_rsqsimcatalog

%% Compare (a) NSHM (all branches) and RSQSim, and (b) NSHM, RSQSim, and stochastic catalogs

%derive counts of magnitudes for different bins

tmp_mag_range=[6.9:0.1:mag_range_nshm2022(end)-0.2]; mag_counts=zeros(length(tmp_mag_range),2);
tmp_mag_counts=zeros(1,3); b_weight=[0.24 0.53 0.23];

fntsize=13;

for kk=1:length(tmp_mag_range)-1
    
    for mm=1:num_branches
 
        total_rup_rate=sum(geo_branch_rupture_rates{mm}.AnnualRate);
        mag_indx=find(geo_branch_rupture_rates{mm}.magnitude>tmp_mag_range(kk) & geo_branch_rupture_rates{mm}.magnitude<tmp_mag_range(kk+1));

        if isempty(mag_indx)==0
            tmp_mag_counts(1,mm)=length(mag_indx)*sum(geo_branch_rupture_rates{mm}.AnnualRate(mag_indx))/total_rup_rate;
        end
    end
    %weight magnitude count for geologic deformation logic tree branches
    mag_counts(kk,1)=tmp_mag_counts(1)*b_weight(1)+tmp_mag_counts(2)*b_weight(2)+tmp_mag_counts(3)*b_weight(3);
    
    %mag_count for rsqsim_catalog
    mag_counts(kk,2)=length(find(otago_rsqsimcatalog_mw69(:,1)>tmp_mag_range(kk) & otago_rsqsimcatalog_mw69(:,1)<tmp_mag_range(kk+1)));

end

%normalise results betweeen 0-1
mag_counts(:,1)=mag_counts(:,1)./max(mag_counts(:,1));
mag_counts(:,2)=mag_counts(:,2)./max(mag_counts(:,2));

figure(6);

% MFD plots for all logic tree branches

tiledlayout(1,2,'TileSpacing','compact');

nexttile

%Plot magnitude counts as bar plots


yyaxis right;  

h1=bar(tmp_mag_range,mag_counts); 
alphaValue = 0.3;
set(h1(1), 'FaceColor', [0.5 0.5 0.5]); set(h1(2), 'FaceColor', [1 0 0]);
set(h1, 'FaceAlpha', alphaValue,'BarWidth',1); 
ylim([0 1]); ylabel('mag count (normalizied)'); set(gca,'YColor',[0 0 0]);

yyaxis left;

%plot weighted mean IFM solution MFD and weighted geologic results
semilogy(mag_range_nshm2022,w_mag_rate,'k-','LineWidth',1.5,'LineStyle','-'); hold on
semilogy(mag_range_nshm2022,geo_b_mag_rate(:,4),'Color',[0.4 0.4 0.4],'LineWidth',1.5,'LineStyle','-'); hold on

%plot rsqsim catalog
semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'r-','LineWidth',1.5); hold on

legend({'IFM weighted geologic','RSQSim','IFM mag count','RSQsim mag count'},'Location','southoutside');

xlim([tmp_mag_range(1)-0.05 tmp_mag_range(end)]); ylim([5*10^-6 7*10^-3]);

xlabel('Magnitude'); set(gca,'YColor',[0 0 0]); 
ylabel('Annual rate of exceedance'); set(gca,'FontSize',fntsize); axis square; grid on

set(gcf, 'Position',[218 203 515 515])


t2=title('(a)','fontsize',fntsize+3,'fontweight','normal');
set(t2, 'horizontalAlignment', 'left');
set(t2, 'units', 'normalized');
h2 = get(t2, 'position');
set(t2, 'position', [-0.15 h2(2) h2(3)]);

nexttile

semilogy(mag_range_nshm2022,w_mag_rate,'k-','LineWidth',1.5); hold on
semilogy(mag_range_nshm2022,geo_b_mag_rate(:,4),'LineWidth',1.5,'Color',[0.4 0.4 0.4]); hold on
semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'r-','LineWidth',1.5); hold on

%plot stochastic catalog (Poisson only)
col_opt=vertcat([1,0.5,0],[0 0 1],[0.2 0.5 0.2]);

for kk=1:3
 semilogy(all_mag_range_GR,fault_syncatRate_poi(:,kk),'Color',col_opt(kk,:),'LineWidth',1.5); hold on; %MFD from poi catalog
end

axis([tmp_mag_range(1)-0.05 tmp_mag_range(end) 5*10^-6 7*10^-3]);

legend({'IFM weighted all','IFM weighted geologic','RSQsim','segmented-char','combined-char','combined-GR'},'Location','Southoutside');

xlabel('Magnitude');ylabel('Annual rate of exceedance'); set(gca,'FontSize',fntsize); axis square; grid on

t2=title('(b)','fontsize',fntsize+3,'fontweight','normal');
set(t2, 'horizontalAlignment', 'left');
set(t2, 'units', 'normalized');
h2 = get(t2, 'position');
set(t2, 'position', [-0.15 h2(2) h2(3)]);


set(gcf, 'Position',[218 203 1102 515])

%% Compare RSQSim and Instrumental Catalog

figure(7);

fntsize=12;

%option to plot ins_MFD (note this only shows events up to M~5)
%1 - plots ins MFD
%2 - ins MFD does't plot
ins_mfd_opt=2;

tiledlayout(1,2,"TileSpacing","compact");
nexttile

if runs2==0
    %fix samples with 0 moment rate to 0.01 so can plot in pdf
    sampleMoRate_rsqsim_70(find(sampleMoRate_rsqsim_70(:,1)==0),1)=1;
    empty_count_rsqsim=length(find(sampleMoRate_rsqsim_70(:,1)==1));
end

semilogy(rsqsim_mag_range,squeeze(sampleAnnualRate_rsqsim_70),'Color',[0.6 0.6 0.6 0.5],'LineWidth',0.5); hold on; %MFD from catalog samples
p6=semilogy(rsqsim_mag_range,sampleAnnualRate_rsqsim_70(1,:),'Color',[0.6 0.6 0.6 0.5],'LineWidth',0.5); hold on; %MFD from catalog samples
p7=semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'k-','LineWidth',1.5); hold on
p8=semilogy(ins_mag_range,ins_AnnualRate,'b-','LineWidth',1.5); hold on; %MFD from G_R_event catalog

set(gca,'fontsize',fntsize); 

if ins_mfd_opt==1
    legend([p7 p8],{[sprintf('RSQSim \n empty samples = '), num2str(empty_count_rsqsim)],'Instrumental'},'Location','southwest');hold on;
    axis([ins_Mmin Mmax 10^-5 10^1]);
else
    legend([p7 p6],{'RSQSim catalog',[sprintf('70 year catalog subsamples,\n empty samples = '), num2str(empty_count_rsqsim)]},'Location','southwest');hold on;
    axis([5 Mmax 10^-5 10^1]);
end
xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;

t1=title('(a)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile
  
%fit kernel distribution
pd_i_rsqsim = fitdist(log10(sampleMoRate_rsqsim_70),'kernel');
x_i_rsqsim = 0:0.05:max(log10(sampleMoRate_rsqsim_70));
y_i_rsqsim = pdf(pd_i_rsqsim,x_i_rsqsim);

plot(x_i_rsqsim,y_i_rsqsim,'k-','linewidth',1.3);hold on
plot([log10(ins_mo_rate) log10(ins_mo_rate)],[0 1],'Color',[1 0.6 0],'linestyle','--','linewidth',1.3);
plot([log10(sum(input_mo_rate)) log10(sum(input_mo_rate))],[0 1],'b--','linewidth',1.3);
 
legend('RSQSim','Instrumental','NZ CFM','Location','northwest');hold on;

set(gca,'fontsize',fntsize); axis square; xlabel('log10 (moment rate, Nm/yr)'); ylabel('f(x)'); 
ylim([0 0.6]); xlim([0 19]);

runs2=1;


t1=title('(b)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.15 h1(2) h1(3)]);

set(gcf,'Position',[440 459 708 338]);

%% Plot URZ forecasts with RSQSim and Stochastic Event Catalog

col_opt=vertcat([1,0.5,0],[0 0 1],[0.2 0.8 0.2]);
rate_count=zeros(3,1); runs2=0;

figure(8);

tiledlayout(2,2,'TileSpacing','tight'); xmax=[33 21]; xmin=[-2 14.8];

mag_check=5; text_opt=["(a)","(b)"]; label_opt=["M>=5 count","log10(Moment Rate)"];
loc_opt=["Southeast","Southeast"];

for kk=1:2

    nexttile

    if kk==1 %distributions based on annual rate of M>5 events
        mag_indx1=find(urz_mag_range==mag_check);%mag indx for URZ forecasts
        mag_indx2=find(all_mag_range_GR==mag_check);%mag indx for NZ CFM forecasts
        
        %calc ecdf and 5th and 95th % for nzcfm based catalogs
        [f1,x1]=ecdf(sampleAnnualRate_rsqsim_50(:,mag_indx2)*50);
        [f2,x2]=ecdf(sampleAnnualRate(:,mag_indx2,1)*50);
         x1_5_95=prctile(sampleAnnualRate_rsqsim_50(:,mag_indx2)*50,[5,95]);
         x2_5_95=prctile(sampleAnnualRate(:,mag_indx2,1)*50,[5,95]);

    else %distributions based on annual rate on moment rate in 50 year samples
        
        %force non-zero moment rates to avoid areas with log10 function
        if runs3==0
            sampleMoRate_rsqsim_50(find(sampleMoRate_rsqsim_50==0))=1;
            sampleMoRate(find(sampleMoRate(:,1)==0),1)=1;
        end

        [f1,x1]=ecdf(log10(sampleMoRate_rsqsim_50));
        [f2,x2]=ecdf(log10(sampleMoRate(:,1)));
        x1_5_95=prctile(log10(sampleMoRate_rsqsim_50(:,mag_indx2)),[5,95]);
        x2_5_95=prctile(log10(sampleMoRate(:,1)),[5,95]);

        urz_mo_rate_samples=zeros(length(urz_mag_rate_sample),2);
      
        for ff=1:2 %for each URZ forecast

            for ss=1:length(urz_mag_rate_sample)%for each sample
             
             mo_rate_mag_bin=zeros(length(urz_mag_range)-1,1);

                 for ii=1:length(urz_mag_range)-1 %for each mag increment
                     inc_rate=urz_mag_rate_sample(ss,ii,ff)-urz_mag_rate_sample(ss,ii+1,ff);%incremental rate for each mag bin
                     mo_rate_mag_bin(ii)=10^(1.5*mean([urz_mag_range(ii),urz_mag_range(ii+1)])+9.05)*inc_rate;
                 end
                    urz_mo_rate_samples(ss,ff)=sum(mo_rate_mag_bin);%moment rate of each sample
            end
        end
        
    end
    
    %Plot 5-95% percentiles and ecdf of nzcfm based catalog moment rate distributions
    plot([x1_5_95(1),x1_5_95(1)],[0 1],'LineWidth',1,'Color','k','LineStyle','--'); hold on
    plot([x1_5_95(2),x1_5_95(2)],[0 1],'LineWidth',1,'Color','k','LineStyle','--'); hold on
    plot([x2_5_95(1),x2_5_95(1)],[0 1],'LineWidth',1,'Color',[0.5 0.5 0.5 1],'LineStyle','--'); hold on
    plot([x2_5_95(2),x2_5_95(2)],[0 1],'LineWidth',1,'Color',[0.5 0.5 0.5 1],'LineStyle','--'); hold on
    p2=plot(x1,f1,'LineWidth',1.3,'Color','k'); hold on
    p3=plot(x2,f2,'LineWidth',1.3,'Color',[0.5 0.5 0.5 1]); hold on
    
    
    %Plot mean rate and moment rate in URZ forecasts
    if kk==1
        for jj=1:2
            p1(jj)=plot([urz_mag_rate(mag_indx1,jj)*50,urz_mag_rate(mag_indx1,jj)*50],[0 1],'--','LineWidth',1.3,'Color',[col_opt(jj,:),0.9]);hold on
        end
    else
        for jj=1:2
            p1(jj)=plot([log10(urz_mo_rate(jj)) log10(urz_mo_rate(jj))],[0 1],'--','LineWidth',1.3,'Color',[col_opt(jj,:),0.9]);hold on
        end
    end
    %}
    set(gca,'FontSize',fntsize+1); ylabel('F(x)'); xlabel(label_opt(kk)); xlim([xmin(kk) xmax(kk)]);
    
    
    if kk==1
        legend([p1(1),p1(2),p2,p3],{['neg-binom: ',num2str(urz_mag_rate(mag_indx1,1)*50,2)],['Poisson: ',num2str(urz_mag_rate(mag_indx1,2)*50,2)],...
            ['RSQSim: ',num2str(mean(sampleAnnualRate_rsqsim_50(:,mag_indx2))*50,2)],...
           ['Stochastic: ',num2str(mean(sampleAnnualRate(:,mag_indx2))*50,2)]},'FontSize',fntsize,'Location',loc_opt(kk));

    else
              legend([p1(1),p1(2),p2,p3],{['neg-binom:',newline,num2str(urz_mo_rate(1),2),' Nm/yr'],['Poisson:',newline,num2str(urz_mo_rate(2),2),' Nm/yr'],...
            ['RSQSim:',newline,num2str(mean(sampleMoRate_rsqsim_50),3),' Nm/yr'],...
            ['Stochastic:',newline,num2str(mean(sampleMoRate(:,1)),3),' Nm/yr']},'FontSize',fntsize-2.5,'Location',loc_opt(kk));

    end
 
    axis square; grid on; box on

    t1=title(text_opt(kk),'fontsize',fntsize+2,'fontweight','normal');
    set(t1, 'horizontalAlignment', 'left');
    set(t1, 'units', 'normalized');
    h1 = get(t1, 'position');
    set(t1, 'position', [-0.15 h1(2) h1(3)]);

end

nexttile %Plot URZ with RSQSim

%plot RSQSim
semilogy(rsqsim_mag_range,squeeze(sampleAnnualRate_rsqsim_50),'Color',[0.5 0.5 0.5 0.3],'LineWidth',0.5); hold on; %MFD from catalog samples
p2=semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'k-','LineWidth',1.2); hold on

%plot overall URZ MFD, derive % of RSQSim samples with moment rate higher than URZ forecast
for ff=1:2
    p1(ff)=semilogy(urz_mag_range,urz_mag_rate(:,ff),'-v','Color',col_opt(ff,:),'LineWidth',1.2,'MarkerIndices',[1+ff:5:length(urz_mag_range)]);hold on
    rate_count(ff)=(length(find(urz_mo_rate(ff)<sampleMoRate_rsqsim_50)))*100/1000; % %of RSQSim samples that had higher M>5 rate than DSM catalog
end

ylabel('Annual frequency of exceedance'); xlabel('Magnitude');xlim([5 8]); ylim([5*10^-5 2]);
set(gca,'fontsize',fntsize+1);
legend([p1(1), p1(2), p2],{['URZ-negbinom, rate test: ', num2str(rate_count(1),3),'%'],...
        ['URZ-poi, rate test: ',num2str(rate_count(2),2),'%'],...
        ['hybrid, rate test: ',num2str(rate_count(3),2),'%'],'RSQSim',...
},'location','northeast','fontsize',fntsize-1)
axis square; grid on

t1=title('(c)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.15 h1(2) h1(3)]);

%note samples taken from Weibull catalogs with aperioodicity of 2 and segmented-char on-fault MFD
if model_opt~=2
    disp('set model_opt to 2 for concistency with sampling')
end

nexttile 

semilogy(all_mag_range_GR,squeeze(sampleAnnualRate(:,:,1)),'Color',[0.5 0.5 0.5 0.3],'LineWidth',0.5); hold on; %MFD from catalog samples
p3=semilogy(all_mag_range_GR,fault_syncatRate_wbl(:,1),'Color',[0.5 0.5 0.5 1],'LineWidth',1.2); hold on; %MFD from poi catalog

%plot overall URZ MFD
for ff=1:2
       p1(ff)=semilogy(urz_mag_range,urz_mag_rate(:,ff),'-v','Color',col_opt(ff,:),'LineWidth',1.2,'MarkerIndices',[1+ff:5:length(urz_mag_range)]);hold on
       rate_count(ff)=(length(find(urz_mo_rate(ff)<sampleMoRate(:,1))))*100/1000; % %of stochastic event catalog samples that had higher M>5 rate than DSM catalog
end

ylabel('Annual frequency of exceedance'); xlabel('Magnitude');xlim([5 8]); ylim([5*10^-5 2]);
set(gca,'fontsize',fntsize+1);
    legend([p1(1), p1(2), p3],{['URZ-negbinom, rate test: ', num2str(rate_count(1),2),'%'],...
        ['URZ-poi, rate test: ',num2str(rate_count(2),2),'%'],...
        ['Stochastic event catalog']},'location','northeast','fontsize',fntsize-1)
axis square; grid on

t1=title('(d)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.15 h1(2) h1(3)]);

set(gcf,'Position', [440 126 653 572]); runs3=1;


%% Plot (a) NSHM with paleoseismic constraints and (b) RSQSim and stochastic event catlogs

%Note model_opt should be set to 2
load('catalogs_stochastic/model2/catalog_poi_10kyr_samples')

inc_mo_rate=zeros(length(mag_range_nshm2022_combined)-1,3);

for jj=1:length(mag_range_nshm2022_combined)-1
    mag_inc_rate=nshm2022_combined_rate(jj,:)-nshm2022_combined_rate(jj+1,:);
    inc_mo_rate(jj,:)=mag_inc_rate.*10.^(((mag_range_nshm2022_combined(jj)+0.025)*b)+c);
end

fntsize=12;

figure(9);

tiledlayout(1,2,"TileSpacing","compact");

nexttile

%plot NSHM inversion MFD (NOT combined with URZ)
semilogy(mag_range_nshm2022,w_mag_rate,'LineWidth',1,'Color',[0 0 0]);hold on %weighted rupture rates
semilogy(mag_range_nshm2022,geo_b_mag_rate(:,4),'LineWidth',1,'Color',[0.4 0.4 0.4 0.5]);hold on %rupture rates for weighed geologic rates

%plot NSHM URZ (negative binomial model)
semilogy(urz_mag_range,urz_mag_rate(:,1),'LineWidth',1.5,'Color',[0.6 0 0.6]);hold on;

%plot NSHM inversion MFD (combined rates with URZ)
semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,1),'LineWidth',1.5,'Color',[0 0 0]);hold on %combined weighted rupture rates
semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,2),'LineWidth',1.5,'Color',[0.4 0.4 0.4]);hold on %combined rupture rates for weighed geologic rates

%Errorbar estimate from paleoseismology constraint as derived through otago_eq_timings_simulations with 95% uncertainity
errorbar(paleoseismic_mag_constraint(:,1),(paleoseismic_rate_constraint(:,1)/time_window_length),(2*paleoseismic_rate_constraint(:,2)/time_window_length),(2*paleoseismic_rate_constraint(:,2)/10000),...
    paleoseismic_mag_constraint(:,2),paleoseismic_mag_constraint(:,2),"o","MarkerFaceColor",[1 1 1],'LineWidth',1.5,...
    'Color',[0.3 0.3 0.3],"MarkerEdgeColor",[0.3 0.3 0.3]);

legend({['IFM weighted', newline, 'mean: ', num2str(ifm_otago_fault_stats{1}(3),3), ' Nm/yr'],...
    ['IFM weighted', newline, 'geologic mean ', num2str(ifm_otago_fault_stats{2}(4,3),3), ' Nm/yr'],...
    ['URZ (negbinom): ', num2str(urz_mo_rate(1),3), ' Nm/yr'],...
    ['URZ-IFM wm', newline, 'combined: ', num2str(sum(inc_mo_rate(:,1)),3), ' Nm/yr'],...
    ['URZ-IFM geol', newline, 'combined: ', num2str(sum(inc_mo_rate(:,2)),3), ' Nm/yr'],...
    ['paleoseismic rate']},'Location','northeast','Fontsize',fntsize-3);hold on;

xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); xlim([6.5,8.1]); ylim([10^-5 2*10^-2]); 
grid on; axis square; set(gca,'FontSize',fntsize);

t1=title('(a)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile

%plot stochastic catalog (Poisson only)
col_opt=vertcat([1,0.5,0],[0 0 1],[0.2 0.5 0.2]);

paleo_mag_lowerb=paleoseismic_mag_constraint(:,1)-paleoseismic_mag_constraint(:,2);
paleo_mag_upperb=paleoseismic_mag_constraint(:,1)+paleoseismic_mag_constraint(:,2);

paleo_rate_lowerb=paleoseismic_rate_constraint(:,1)-2*paleoseismic_rate_constraint(:,2);
paleo_rate_upperb=paleoseismic_rate_constraint(:,1)+2*paleoseismic_rate_constraint(:,2);

nsamples=length(sampleAnnualRate_10kyr(:,1,:,1)); 
paleoseis_mag=randi([paleo_mag_lowerb*200 paleo_mag_upperb*200],nsamples,1)/200; paleosmag_rate=zeros(nsamples,1);

for kk=1:3
   semilogy(all_mag_range_GR,squeeze(sampleAnnualRate_10kyr(:,:,kk)),'Color',[col_opt(kk,:),0.3],'LineWidth',0.12); hold on; %10 kyr subsamples
    
   %for each catalog sample, randomly sample a mag between constraints
   %then see if catalog sample is within the rate range
    for nn=1:nsamples
        paleosmag_indx=find(round(paleoseis_mag(nn),3)==round(all_mag_range_GR,3));
        paleosmag_rate(nn)=sampleAnnualRate_10kyr(nn,paleosmag_indx,kk); %rate for randomly sampled mag constraint
    end
   
   paleoseis_check(kk)=length(find(paleosmag_rate>=paleo_rate_lowerb/time_window_length & paleosmag_rate<=paleo_rate_upperb/time_window_length));
end

for kk=1:3
   pp(kk)=semilogy(all_mag_range_GR,fault_syncatRate_poi(:,kk),'Color',col_opt(kk,:),'LineWidth',1.5); hold on; %MFD from poi catalog
end

%plot RSQSIM analysis
p4=semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'r-','LineWidth',1.5); hold on
semilogy(all_mag_range_GR,squeeze(sampleAnnualRate_rsqsim_10kyr),'Color',[1 0 0 0.3],'LineWidth',0.12); hold on; %10 kyr subsamples

for nn=1:nsamples
        paleosmag_indx=find(round(paleoseis_mag(nn),3)==round(all_mag_range_GR,3));
        paleosmag_rate(nn)=sampleAnnualRate_rsqsim_10kyr(nn,paleosmag_indx); %rate for randomly sampled mag constraint
end

paleoseis_check(kk+1)=length(find(paleosmag_rate>=paleo_rate_lowerb/time_window_length & paleosmag_rate<=paleo_rate_upperb/time_window_length));

%Errorbar estimate from paleoseismology constraint as derived through otago_eq_timings_simulations with 95% uncertainity
p5=errorbar(paleoseismic_mag_constraint(:,1),(paleoseismic_rate_constraint(:,1)/time_window_length),(2*paleoseismic_rate_constraint(:,2)/time_window_length),(2*paleoseismic_rate_constraint(:,2)/time_window_length),...
    paleoseismic_mag_constraint(:,2),paleoseismic_mag_constraint(:,2),"o","MarkerFaceColor",[1 1 1],'LineWidth',1.5,...
    'Color',[0.2 0.2 0.2],"MarkerEdgeColor",[0.2 0.2 0.2]);

legend([pp(1),pp(2),pp(3),p4,p5],{['segmented-char-poi: ',num2str(MoRate_StochasticEventCatalog_all_poi(1),3), ' Nm/yr: ',num2str(paleoseis_check(1)/10,'%.0f'), '%'],...
    ['combined-char-poi: ', num2str(MoRate_StochasticEventCatalog_all_poi(2),3), ' Nm/yr: ', num2str(paleoseis_check(2)/10,'%.0f'), '%'],...
    ['combined-GR-poi: ', num2str(MoRate_StochasticEventCatalog_all_poi(3),3), ' Nm/yr: ',num2str(paleoseis_check(3)/10,'%.0f'), '%'],...
    ['RSQSim: ', num2str(rsqsim_MoRate_catalog,3), ' Nm/yr: ',num2str(paleoseis_check(4)/10,'%.0f'), '%'],...
    ['paleoseismic rate'],...
    },'fontsize',fntsize-3,'Location','northeast');hold on;

xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); xlim([6.5,8.1]); ylim([10^-5 2*10^-2]); 
grid on; axis square; set(gca,'FontSize',fntsize);

t1=title('(b)','fontsize',fntsize+2,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.15 h1(2) h1(3)]);

set(gcf,'Position',[318 378 830 419]);

%% Compare all catalogs

%option to plot ins_MFD (note this only shows events up to M~5)
%1 - plots ins MFD
%2 - ins MFD does't plot
ins_mfd_opt=2;

figure(10);

mstyl=["v","^","*"]; marker_space{1}=[1:50:length(all_mag_range_GR)];
marker_space{2}=[25:50:length(all_mag_range_GR)];

%plot stochastic catalog (Poisson only)
col_opt=vertcat([1,0.5,0],[0 0 1],[0.2 0.5 0.2]);


for kk=1:3
 semilogy(all_mag_range_GR,fault_syncatRate_poi(:,kk),'Color',col_opt(kk,:),'LineWidth',1.5); hold on; %MFD from poi catalog
end

%plot RSQSIM analysis
semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'r-','LineWidth',1.5); hold on

%plot NSHM inversion MFD (NOT combined with URZ)
%semilogy(mag_range_nshm2022,w_mag_rate,'k-','LineWidth',1.5);hold on %weighted rupture rates
semilogy(mag_range_nshm2022,geo_b_mag_rate(:,4),'LineWidth',1.5,'Color',[0.4 0.4 0.4]);hold on %rupture rates for weighted geologic mean

%plot NSHM DSM (for negative binomial URZ)
semilogy(urz_mag_range,urz_mag_rate(:,1),'LineWidth',1.5);hold on;

%Uncomment, and change legend option if want to plot this instead
%plot NSHM inversion MFD (combined rates with URZ)
semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,1),'Color',[0.5 0 1],'LineWidth',1.5);hold on %weighted rupture rates
semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,2),'b-','LineWidth',1.5);hold on %rupture rates for central geologic mean


%Errorbar estimate from paleoseismology constraint as derived through otago_eq_timings_simulations with 95% uncertainity
errorbar(7,(paleoseismic_rate_constraint(:,1)/10000),(2*paleoseismic_rate_constraint(:,2)/10000),(2*paleoseismic_rate_constraint(:,2)/10000),0.2,0.2,...
    "o","MarkerFaceColor",[1 1 1],'LineWidth',1.5,...
    'Color',[0.3 0.3 0.3],"MarkerEdgeColor",[0.3 0.3 0.3]);

%plot instrumental MFD (if ins_mdf_opt=1)
if ins_mfd_opt==1
    semilogy(ins_mag_range,ins_AnnualRate,'b-','LineWidth',1.5); hold on; %MFD from G_R_event catalog
    xlim([ins_Mmin,max(mag_range_nshm2022)]); ylim([10^-5 10^1]);
else
    xlim([6,8.1]); ylim([10^-5 2*10^-2]);
end

set(gca,'fontsize',fntsize+1);  

legend({['segmented-char-poi: ',num2str(MoRate_StochasticEventCatalog_all_poi(1),3), ' Nm/yr'],...
    ['combined-char-poi: ', num2str(MoRate_StochasticEventCatalog_all_poi(2),3), ' Nm/yr'],...
    ['combined-GR-poi: ', num2str(MoRate_StochasticEventCatalog_all_poi(3),3), ' Nm/yr'],...
    ['RSQSim: ', num2str(rsqsim_MoRate_catalog,3), ' Nm/yr'],...
    ['2022 NSHM inversion', newline, 'weighted mean: ', num2str(ifm_otago_fault_stats{1}(3),3), ' Nm/yr'],...
    ['2022 NSHM inversion', newline, 'geologic dm: ', num2str(ifm_otago_fault_stats{2}(4,3),3), ' Nm/yr'],...
    ['2022 NSHM URZ: ', num2str(urz_mo_rate(1),3), ' Nm/yr'],...
    ['2022 NSHM combined: ', num2str(sum(inc_mo_rate(:,1)),3), ' Nm/yr'],...
    ['10 kyr paleoseismic rate'],...
    %['Instrumental: ', num2str(ins_mo_rate,3), ' Nm/yr']},...
   },'Location','southwest');hold on;


xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;
set(gcf,'Position',[680 496 611 481]);

%% Plot to illustrate moment rate reduction associated with area weighted moment

figure(11);

tiledlayout(2,2,'TileSpacing','compact');

nexttile

semilogy(mag_range_nshm2022,mag_rate,'r-','LineWidth',1.2);hold on
semilogy(mag_range_nshm2022,w_mag_rate,'k-','LineWidth',1.2);

legend('All IFM ruptures','Scaled by area within Otago')
axis([6.9 8.2  5*10^-6 5*10^-3]);

xlabel('Magnitude');ylabel('Annual rate of exceedance'); axis square; grid on

t1=title('(a)','fontsize',12,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile

%Comparing RSQSim after downweighting magnitudes
semilogy(rsqsim_mag_range,rsqsim_mfd_rate,'r-','LineWidth',1.5); hold on;

semilogy(rsqsim_mag_range,rsqsim_partial_mfd_rate,'k-','LineWidth',1.5); hold on;
legend('All RSQSim catalog','Scaled by moment within Otago');

axis([5 8 5*10^-5 5*10^-1]);

xlabel('Magnitude');ylabel('Annual rate of exceedance'); axis square; grid on

t1=title('(b)','fontsize',12,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.15 h1(2) h1(3)]);

% MFD plots for all logic tree branches
%(c) show for geologic deformation branches
%(d) show for geodetic deformation branches

col_opt=vertcat([0 0 1],[1 0 1],[1,0.5,0],[0.2 0.5 0.2]);

title_opt=["(c)","(d)"]; txt_opt=["geol weighted mean","geod weighted mean"];
tmp_width=[1 1 1 1.5];

tmp_mag_range=[6.9:0.1:mag_range_nshm2022(end)-0.4]; 

for hh=1:2
    

    if hh==1 %plot geologic rates
        tmp_rates=geo_b_mag_rate; tmp_branch_indx=geo_branch_indx; tmp_branch_rupture_rates=geo_branch_rupture_rates;
    else    %plot geodetic rates
        tmp_rates=ged_b_mag_rate; tmp_branch_indx=ged_branch_indx; tmp_branch_rupture_rates=ged_branch_rupture_rates;
    end
    
    nexttile
    
    yyaxis right
    
    mag_counts=zeros(length(tmp_mag_range),length(num_branches));

    for kk=1:length(tmp_mag_range)-1
        
        for ii=1:num_branches
        
            total_rup_rate=sum(tmp_branch_rupture_rates{ii}.AnnualRate);

            %find number of rupture magnitudes in each bin
            mag_indx=find(tmp_branch_rupture_rates{ii}.magnitude>tmp_mag_range(kk) & tmp_branch_rupture_rates{ii}.magnitude<tmp_mag_range(kk+1));

            %weight magnitude count by rupture rate
            if isempty(mag_indx)==0
                mag_counts(kk,ii)=length(mag_indx)*sum(tmp_branch_rupture_rates{ii}.AnnualRate(mag_indx))/total_rup_rate;
            end

        end %end ii loop for branches
    end %end kk loop for mag_count

    b1=bar(tmp_mag_range,mag_counts,'stacked');
    alphaValue = 0.3; set(b1, 'FaceAlpha', alphaValue);
    
    b1(1).FaceColor=col_opt(1,:); b1(2).FaceColor=col_opt(2,:); b1(3).FaceColor=col_opt(3,:);

    ylabel('Weighted Count'); set(gca,'YColor',[0 0 0]);

    yyaxis left

    semilogy(mag_range_nshm2022,w_mag_rate,'k-','LineWidth',1.5,'LineStyle','-'); hold on
    tmp_txt=cell(num_branches,1);

    legend_txt=["weighted mean rupture rates"];
    
    for ii=1:num_branches+1
    
        semilogy(mag_range_nshm2022,tmp_rates(:,ii),'Color',col_opt(ii,:),'LineWidth',tmp_width(ii),'LineStyle','-','Marker', 'none'); hold on
        %obtain deformation model, B-n pair, and weighting for each branch
        if ii<=num_branches
            tmp_txt{ii}=strjoin({cell2mat(logic_tree_branches.deformation_model(tmp_branch_indx(ii),:)),...
                cell2mat(logic_tree_branches.bN_pair(tmp_branch_indx(ii),:)),num2str(b_weight(ii))},' ');
        else
            tmp_txt{ii}=txt_opt(hh);
        end

    legend_txt=[legend_txt,tmp_txt{ii}];
    
    end

    legend(legend_txt,'Location','southoutside','fontsize',8);

    axis([6.8 8 5*10^-6 7*10^-3]); set(gca,'YColor',[0 0 0]);

    xlabel('Magnitude');ylabel('Annual rate of exceedance'); axis square; grid on

    t2=title(title_opt(hh),'fontsize',12,'fontweight','normal');
    set(t2, 'horizontalAlignment', 'left');
    set(t2, 'units', 'normalized');
    h2 = get(t2, 'position');
    set(t2, 'position', [-0.15 h2(2) h2(3)]);

end

set(gcf,'Position',[770 156 878 722])