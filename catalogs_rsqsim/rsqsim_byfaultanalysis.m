%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Analysis of individual faults in Otago-based RSQSim Catalog %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%should be consistent with mu in rupture_patches.py and rigidity
%in catalogs_stocastic/fault_recurrence_parameters
rigidity=3e10;
load('catalog_rsqsim_statistics','rsqsim_mag_range');

%Index events by fault in catalog based on the rupture patches output
addpath('otago_rsqsim_catalog'); addpath('otago_rsqsim_catalog/fault_catalog');
num_simu_rsqsim=1e6; yearsec=365.25*24*3600;

c=9.05; b=1.5; %scaling between moment and magnitude

%load NZ CFM info about faults
orb_faults=readtable('OtagoRangeBasinFaults.xlsx','Sheet','segmented'); 
nzcfm_fault_id=orb_faults.Fault_ID;

%load faults assessed in RSQSIM catalog
otago_faults=readtable('otago_fault_list_20240428.csv'); 

%% Loop through all fault catalogs to obtain moment rate and interevent time distribution

fault_mo=zeros(length(otago_faults.name),2); fault_cov=zeros(length(otago_faults.name),2);
fault_int=cell(length(otago_faults.name),1); fault_sr_int=cell(length(otago_faults.name),1);

rsqsim_fault_stats=zeros(length(otago_faults.name),3);

for ii=1:length(otago_faults.name)
    
    catalog=readtable([char(otago_faults.name(ii)) '.csv']);
    
    %moment rate of fault in catalog
    fault_mo(ii,1)=sum(10.^(catalog.partial_mw*b+c))/num_simu_rsqsim;
    
    %find index of Otago fault in NZ CFM info, and derive NZ CFM slip rate
    tmp_indx=find(nzcfm_fault_id==otago_faults.Fault_ID(ii));
    
    fault_mo(ii,2)=orb_faults.SR_pref(tmp_indx)*10^-3*orb_faults.Area_m2(tmp_indx)*rigidity;
    
    fault_int{ii}=zeros(length(catalog.t0)-1,1);
    
    %derive interevent time distributions for all fault events
    for jj=1:length(catalog.t0)-1
        
        fault_int{ii}(jj)=(catalog.t0(jj+1)-catalog.t0(jj))/yearsec;
    end
    
    %derive CoV for all fault events
    fault_cov(ii,1)=std(fault_int{ii})/mean(fault_int{ii});

    %derive interevent time distributions for surface_rup events
    sr_indx=find(catalog.surface_rupture==1);
    
    %only derive interevent time distributions if fault has hosted >=5 surface ruptures
    if length(sr_indx)>=5
        
        fault_sr_int{ii}=zeros(length(sr_indx-1),1);
        for jj=1:length(sr_indx)-1
            fault_sr_int{ii}(jj)=(catalog.t0(sr_indx(jj+1))-catalog.t0(sr_indx(jj)))/yearsec;
        end
        
        %derive CoV for all surface rupturing events
        fault_cov(ii,2)=std(fault_sr_int{ii})/mean(fault_sr_int{ii});
    end
    
    % multifaults events for M>6.9 ruptures for each fault
    rsqsim_fault_stats(ii,1)=length(find(catalog.num_fault>1 & catalog.mw>6.9))/length(find(catalog.mw>6.9));

    %average and std magnitude of event that fault partipcates in (and contributes seismic moment equivalent to 6.7 event)
    rsqsim_fault_stats(ii,2)=mean(catalog.mw(find(catalog.partial_mw>6.7)));
    rsqsim_fault_stats(ii,3)=std(catalog.mw(find(catalog.partial_mw>6.7)));
end

%correction for Titri fault which has different extents in the NZ CFM and RSQSim
tmp_indx1=42; %index for 'titricombined' in otago_faults
tmp_indx2=find(nzcfm_fault_id==804);  %index for Tiri Central in NZ CFM
%add additional moment from Titri central into NZ CFM moment rate
fault_mo(tmp_indx1,2)=fault_mo(tmp_indx1,2)+orb_faults.SR_pref(tmp_indx2)*10^-3*orb_faults.Area_m2(tmp_indx1)*rigidity;

%% Derive moment rate of faults in NZ CFM following IFM slip rate reduction

addpath([mydir(1:idcs(end)-1),'/nshm_inversion/mfd_analysis/by_fault']);

%add nshm fault stats from weighted mean logic tree branch (5) and weighted geolgoc logic tree branches(9)
%Files are created in nshm_inversion/mfd_analysis/nshm_otago_inversion_results_by_fault.m
load nshm_fault_stats5; nshm_fault_stats5=nshm_fault_stats; clear nshm_fault_stats
load nshm_fault_stats9; nshm_fault_stats9=nshm_fault_stats; clear nshm_fault_stats

otago_faults_nshm=readtable([mydir(1:idcs(end)-1),'/nshm_inversion/mfd_analysis/otago_fault_list_20240428.csv']);
nshm_fault_nzcfm_mo_rate=zeros(height(otago_faults_nshm),1);

for ii=1:height(otago_faults_nshm)
    
    tmp_indx=find(otago_faults_nshm.Fault_ID(ii)==orb_faults.Fault_ID);
    %only use 90% of slip rate due to moment released at M<Mmin (Gerstenbgerger et al 2024)
    %correction for area already made (see orb_fault_geometries/faultgeometries.m)
    nshm_fault_nzcfm_mo_rate(ii)=orb_faults.SR_pref(tmp_indx)*10^-3*0.9*orb_faults.Area_m2(tmp_indx)*rigidity;
end

%remove path so doesn't get confused by fault catalog to search for
rmpath([mydir(1:idcs(end)-1),'/nshm_inversion/mfd_analysis/by_fault']);


%% Plot distribution of fault magnitudes for RSQSim and the IFM
figure(1);
num_fault=43;

tiledlayout(2,1,'TileSpacing','compact')

title_opt=["(a)","(b)"];

for ii=1:2

    nexttile
    if ii==1
        nshm_fault_stats=nshm_fault_stats9; %load results from central geologic deformation model branch
    elseif ii==2
        nshm_fault_stats=nshm_fault_stats5; %load results from weighted mean results
    end

    %plot nshm mags
    errorbar([1:1:num_fault],nshm_fault_stats([1:1:41,44,45],3),nshm_fault_stats([1:1:41,44,45],3)-nshm_fault_stats([1:1:41,44,45],4),nshm_fault_stats([1:1:41,44,45],5)-nshm_fault_stats([1:1:41,44,45],3),...
        "o","MarkerFaceColor",[1 1 1],'LineWidth',1.2,'Color',[0 0 0],"MarkerEdgeColor",[0 0 0]);hold on
    %plot for RSQSim
    errorbar([1:1:num_fault],rsqsim_fault_stats([1:1:41,43,44],2),rsqsim_fault_stats([1:1:41,43,44],3),...
        "o","MarkerFaceColor",[1 1 1],'LineWidth',1.2,'Color',[1 0 0],"MarkerEdgeColor",[1 0 0]);hold on;

    xticks([1:1:num_fault]);
    xticklabels(otago_faults.name([1:1:41,43,44]));

    ax=gca; ax.XAxis.FontSize = 7.5; xtickangle(60); ylim([6.6 8.35]);

    legend('IFM','RSQSim','Location','southeast');


    t1=title(title_opt(ii),'fontsize',12,'fontweight','normal');
    set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
    h1 = get(t1, 'position'); set(t1, 'position', [-0.1 h1(2) h1(3)]);
end

set(gcf,'Position',[430 89 686 708]);

%% Plot all fault analysis figures 

axislimits=[14 16.5];

figure(2);

tiledlayout(2,2,'TileSpacing','compact')

%plot RSQSim moment rate comparison
nexttile

plot([axislimits(1) axislimits(2)],[axislimits(1) axislimits(2)],'k--'); hold on
plot(log10(fault_mo(:,2)),log10(fault_mo(:,1)),'kx','MarkerSize',12,'LineWidth',1.5);hold on
axis([axislimits(1) axislimits(2) axislimits(1)  axislimits(2)]); axis square; hold on; 

xlabel('NZ CFM Moment Rate (log Nm/yr)'),ylabel('RSQSim Catalog Moment Rate (log Nm/yr)');
grid on; set(gca,'fontsize',11)

t1=title('(a)','fontsize',12,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile

%plot NSHM moment rate comparison (for central logic tree branch)
plot([axislimits(1) axislimits(2)],[axislimits(1) axislimits(2)],'k--'); hold on
plot(log10(nshm_fault_nzcfm_mo_rate),log10(nshm_fault_stats9(:,1)),'kx','MarkerSize',12,'LineWidth',1.5);hold on
axis([axislimits(1) axislimits(2) axislimits(1)  axislimits(2)]); axis square; hold on; 

xlabel('NZ CFM Moment Rate (log Nm/yr)'),ylabel(['NZ NSHM 2022 IFM weighted geologic',char(10), 'mean Moment Rate (log Nm/yr)']);
grid on; set(gca,'fontsize',11)

t1=title('(b)','fontsize',12,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile

%bar plot for the number of ruptures in each inversion solution

load([mydir(1:idcs(end)-1),'/nshm_inversion/mfd_analysis/nshm_otago_inversion_results'],'geo_branch_rupture_rates',...
    'ged_branch_rupture_rates','logic_tree_branches','geo_branch_indx','ged_branch_indx','num_branches');

tmp_txt_geo=cell(1,3);tmp_txt_ged=cell(1,3);

for mm=1:num_branches

    geo_rup_count(mm)=height(geo_branch_rupture_rates{mm});
    ged_rup_count(mm)=height(ged_branch_rupture_rates{mm});
    
    tmp_txt_geo{mm}=strcat(cell2mat(logic_tree_branches.deformation_model(geo_branch_indx(mm),:)),...
        cell2mat(logic_tree_branches.bN_pair(geo_branch_indx(mm),:)));
    tmp_txt_ged{mm}=strcat(cell2mat(logic_tree_branches.deformation_model(ged_branch_indx(mm),:)),...
        cell2mat(logic_tree_branches.bN_pair(ged_branch_indx(mm),:)));
end


b1=bar(horzcat(geo_rup_count,ged_rup_count));
alphaValue = 0.3; set(b1, 'FaceAlpha', alphaValue);
xticklabels([tmp_txt_ged,tmp_txt_geo]); ax1=gca; 
ax1.XAxis.FontSize = 8.5; xtickangle(45);

xlabel('logic tree branch','fontsize',11); ylabel('number of ruptures'); axis square; 

t1=title('(c)','fontsize',12,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

nexttile

%CoV Histogram for all events on each fault

histogram(fault_cov(:,1),'BinWidth',0.25); hold on
tmp_indx=find(fault_cov(:,2)>0); histogram(fault_cov(tmp_indx,2),'BinWidth',0.25);
xlabel('CoV'),ylabel('Count');ylim([0 20]); xlim([0 8]);set(gca,'fontsize',11); axis square

t1=title('(d)','fontsize',12,'fontweight','normal'); legend('all events','surface rupturing events')
set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

set(gcf,'position',[562 87 846 710]);


%% Fault specific analysis

%fault_select=["nwcardronanorth","nwcardronasouth","titrisouth"];
%fault_select_nshm=["NW Cardrona North","NW Cardrona South","Titri South"];
fault_select=["akatore","dunstan","pisa"];
fault_select_nshm=["Akatore","Dunstan","Pisa"];

window_length=70000; %time window sampling   

p_mag_rate=zeros(length(rsqsim_mag_range),length(fault_select));
mag_rate=zeros(length(rsqsim_mag_range),length(fault_select));
%Legend1=cell(length(fault_select),1); Legend2=cell(length(fault_select),1);
col_opt=vertcat([1 0 0],[0 0 0],[0 0 1]);%may need to update depending on number of faults selected

label_opt=["(a)","(b)","(c)","(d)","(e)","(f)"]; count=0;

tmp_mag_range=[6.9:0.1:8.4];

figure(3);

tiledlayout(length(fault_select),2,'TileSpacing','compact')

for jj=1:length(fault_select)
    
    %select random time interval to sample
    tmp=rand(1)*(num_simu_rsqsim-window_length); 
    time_int=[tmp tmp+window_length];
    
    nexttile; count=count+1;
    
    tmp_indx(jj)=find(fault_select(jj)==otago_faults.name);
    
    %open fault specific catalog
    catalog=readtable([char(fault_select(jj)) '.csv']);
    %correct for NSHM moment-magnitude scaling
    catalog.mw=(log10(catalog.m0)-c)/b; 

    %derive annual frequency of exceedance for each mag in RSQSim catalog
    for kk=1:length(rsqsim_mag_range)  
        mag_rate(kk,jj)=length(find(catalog.mw >= rsqsim_mag_range(kk)))/(num_simu_rsqsim);
        p_mag_rate(kk,jj)=length(find(catalog.partial_mw >= rsqsim_mag_range(kk)))/(num_simu_rsqsim);  
    end
    
    addpath([mydir(1:idcs(end)-1),'/nshm_inversion/mfd_analysis/by_fault']);
    
    %derive annual frequency of exceedance for each mag in NSHM 2022 IFM
    tmp_table=readtable(strjoin(fault_select_nshm(jj),'.csv'));%read in fault data
    
    %remove path so doesn't get confused by fault catalog to search for
    rmpath([mydir(1:idcs(end)-1),'/nshm_inversion/mfd_analysis/by_fault']);
    nshm_mag_range=[6.5:0.05:8.35];
    
    for kk=1:length(nshm_mag_range)
      
        %MFD rate for full rupture magntitudes
        rup_indx1{kk}=find(tmp_table.rup_mw>nshm_mag_range(kk));
     
        mfd_rate1(kk,:)=[sum(tmp_table.weighted_mean_rate(rup_indx1{kk})),sum(tmp_table.weighted_geologic_mean(rup_indx1{kk})),...
            sum(tmp_table.weighted_geodetic_mean(rup_indx1{kk}))];
     
        %MFD rate for partial magnitudes
        rup_indx2{kk}=find(tmp_table.rup_weighted_mw>nshm_mag_range(kk));
        mfd_rate2(kk,:)=[sum(tmp_table.weighted_mean_rate(rup_indx2{kk})),sum(tmp_table.weighted_geologic_mean(rup_indx2{kk})),...
            sum(tmp_table.weighted_geodetic_mean(rup_indx2{kk}))];
     
    end
    
    %total rate for weighted geologic mean rate
    total_rup_rate=sum(tmp_table.weighted_geologic_mean); mag_counts=zeros(length(tmp_mag_range),2);
    
    %derive counts for IFM magntiude bins
    for kk=1:length(tmp_mag_range)-1
        
       mag_indx=find(tmp_table.rup_mw>tmp_mag_range(kk) & tmp_table.rup_mw<tmp_mag_range(kk+1));
        %Weight bins by rupture rate
       if isempty(mag_indx)==0
            mag_counts(kk,1)=length(mag_indx)*sum(tmp_table.weighted_geologic_mean(mag_indx))/total_rup_rate;
       end

    %mag_count for rsqsim_catalog
    mag_counts(kk,2)=length(find(catalog.mw>tmp_mag_range(kk) & catalog.mw<tmp_mag_range(kk+1)));

    end
    
    %normalise results betweeen 0-1
    mag_counts(:,1)=mag_counts(:,1)./max(mag_counts(:,1));
    mag_counts(:,2)=mag_counts(:,2)./max(mag_counts(:,2));

    yyaxis right

    h1=bar(tmp_mag_range,mag_counts); 
    alphaValue = 0.3;
    set(h1(1), 'FaceColor', [0.5 0.5 0.5]); set(h1(2), 'FaceColor', [1 0 0]);
    set(h1, 'FaceAlpha', alphaValue,'BarWidth',0.8); 
    ylim([0 1]); ylabel('mag count (normalizied)'); set(gca,'YColor',[0 0 0]);

    yyaxis left

    %plot RSQSim fault MFD
    semilogy(rsqsim_mag_range,mag_rate(:,jj),'Color',col_opt(1,:),'LineWidth',1.5,'LineStyle','-','Marker','none'); hold on %plot total MFD
    semilogy(rsqsim_mag_range,p_mag_rate(:,jj),'Color',[col_opt(1,:) 0.5],'LineWidth',1.5,'LineStyle','-','Marker','none'); hold on %plot MFD for partial moment release
    
    %plot NSHM, for weighted geologic results
    semilogy(nshm_mag_range,mfd_rate1(:,2),'Color',col_opt(2,:),'LineWidth',1.5,'LineStyle','-','Marker','none'); hold on %plot total MFD
    semilogy(nshm_mag_range,mfd_rate2(:,2),'Color',[col_opt(2,:) 0.5],'LineWidth',1.5,'LineStyle','-','Marker','none'); hold on %plot MFD for partial moment release
    
    
    xlabel('Magnitude');ylabel('Annual rate of exceedance'); set(gca,'YColor',[0 0 0]);
        
    lgd1=legend({'RSQSim: full mw','RSQSim: area-weighted mw','NSHM IFM: full mw',...
        'NSHM IFM: area weighted mw','nshm mag count','rsqsim mag count'},'FontSize',6,'Location','Westoutside');hold on
    lgd1.Title.String=fault_select(jj);
        
    xlim([6.5 8.35] ); ylim([5*10^-6 10^-3]); 
    axis square; grid on; set(gca,'Fontsize',8);
        
     t1=title(label_opt(count),'fontsize',12,'fontweight','normal');
     set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
     h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);
        

    nexttile; count=count+1;
    
    %find events in time interval
    t_indx=find(catalog.t0/yearsec>=time_int(1) & catalog.t0/yearsec<time_int(2));
    
    %find surface rupturing in time interval
    sr_t_indx=find(catalog.surface_rupture(t_indx)==1);
    
    %plot timeline of events
    
    %plot surface rupturing events only
    p1=scatter(catalog.t0(t_indx(sr_t_indx))/yearsec,catalog.partial_mw(t_indx(sr_t_indx)),'o','SizeData',35,'LineWidth',0.7,'MarkerEdgeColor',col_opt(3,:),'MarkerFaceColor',col_opt(3,:)); hold on
    p1.MarkerFaceAlpha = 0.7;
    %plot all events
    p2=scatter(catalog.t0(t_indx)/yearsec,catalog.partial_mw(t_indx),'Marker','o','SizeData',20,'LineWidth',1.3,'MarkerEdgeColor',col_opt(3,:),'MarkerEdgeAlpha',0.5); hold on
    
    label_tmp=vertcat(fault_select(jj),strcat('Mean RI: ',num2str(mean(fault_sr_int{tmp_indx(jj)}),2),' years') , strcat('CoV: ',num2str(fault_cov(tmp_indx(jj),2),2)));
    
    xlabel('Event time (years)');ylabel('Magnitude');%legend(fault_select(1:end));
        
    lgd2=legend([p2 p1],{'all events','sr event'},'FontSize',6);hold on
    lgd2.Title.String=(label_tmp);    
        
     axis([time_int(1) time_int(2) 5 rsqsim_mag_range(end)+0.5]); 
     axis square; grid on; set(gca,'Fontsize',8);
        
      t1=title(label_opt(count),'fontsize',10,'fontweight','normal');
      set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
        h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);
        box on
    
end

set(gcf,'position',[250 87 617 710]);

%% Plots results from different faults in same plot 

fault_select=["nwcardronasouth","nwcardronanorth","pisa"];
%fault_select=["akatore","titrisouth","titricombined"];
%select time interval interested in
time_int=[6.0*1e5 6.75*1e5];    

p_mag_rate=zeros(length(rsqsim_mag_range),length(fault_select));
mag_rate=zeros(length(rsqsim_mag_range),length(fault_select));
Legend1=cell(length(fault_select),1); Legend2=cell(length(fault_select),1);
col_opt=vertcat([1 0 0],[0 0 0],[0 0 1]);%length of this must equal number of faults selected

figure(3);

tiledlayout(1,2,'TileSpacing','compact')

nexttile

for jj=1:length(fault_select)
    
    tmp_indx(jj)=find(fault_select(jj)==otago_faults.name);
    
    %open fault specific catalog
    catalog=readtable([char(fault_select(jj)) '.csv']);
    
    %derive annual frequency of exceedance for each mag
    for kk=1:length(rsqsim_mag_range)  
        mag_rate(kk,jj)=length(find(catalog.mw >= rsqsim_mag_range(kk)))/(num_simu_rsqsim);
        p_mag_rate(kk,jj)=length(find(catalog.partial_mw >= rsqsim_mag_range(kk)))/(num_simu_rsqsim);  
    end
    
    %plot fault MFD
    semilogy(rsqsim_mag_range,mag_rate(:,jj),'Color',col_opt(jj,:),'LineWidth',1.2,'LineStyle','--'); hold on %plot total MFD
    p2=semilogy(rsqsim_mag_range,p_mag_rate(:,jj),'Color',col_opt(jj,:),'LineWidth',1.2,'LineStyle','-'); hold on %plot MFD for partial moment release
    
    Legend1{jj}=p2;%store Matlab line values for legend
    
    if jj==length(fault_select)
        xlabel('Magnitude');ylabel('Annual rate of exceedance');
        
        legend_tmp=Legend1{1};
        for mm=2:length(fault_select)
            legend_tmp=[legend_tmp,Legend1{mm}];
        end
        
        legend(legend_tmp,fault_select(1:end));hold on
        
        axis([rsqsim_mag_range(1) rsqsim_mag_range(end) 10^-5 10^-1]); 
        axis square; grid on; set(gca,'Fontsize',11);
        
        t1=title('(a)','fontsize',12,'fontweight','normal');
        set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
        h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);
        
    end

end

nexttile
    
for jj=1:length(fault_select)   

    catalog=readtable([char(fault_select(jj)) '.csv']);
    
    %find events in time interval
    t_indx=find(catalog.t0/yearsec>=time_int(1) & catalog.t0/yearsec<time_int(2));
    
    %find surface rupturing in time interval
    sr_t_indx=find(catalog.surface_rupture(t_indx)==1);
    
    %plot timeline of events
    
    %plot surface rupturing events only
    p1=scatter(catalog.t0(t_indx(sr_t_indx))/yearsec,catalog.partial_mw(t_indx(sr_t_indx)),'o','SizeData',35,'LineWidth',0.7,'MarkerEdgeColor',col_opt(jj,:),'MarkerFaceColor',col_opt(jj,:)); hold on
    p1.MarkerFaceAlpha = 0.7;
    %plot all events
    p2=scatter(catalog.t0(t_indx)/yearsec,catalog.partial_mw(t_indx),'Marker','o','SizeData',20,'LineWidth',1.3,'MarkerEdgeColor',col_opt(jj,:),'MarkerEdgeAlpha',0.5); hold on
    Legend2{jj}=p2;%store Matlab line values for legend
    
    label_tmp(jj)=strcat(fault_select(jj), ' Mean RI: ',num2str(mean(fault_sr_int{tmp_indx(jj)}),2),' years, CoV: ',num2str(fault_cov(tmp_indx(jj),2),2));
    
       if jj==length(fault_select)
        xlabel('Event time (years)');ylabel('Magnitude');%legend(fault_select(1:end));
        
        legend_tmp=Legend2{1};
        for mm=2:length(fault_select)
            legend_tmp=[legend_tmp,Legend2{mm}];
        end
        
        legend(legend_tmp,label_tmp(1:end));hold on
        
        
        axis([time_int(1) time_int(2) 5 rsqsim_mag_range(end)+0.5]); 
        axis square; grid on; set(gca,'Fontsize',11);
        
        t1=title('(b)','fontsize',12,'fontweight','normal');
        set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
        h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);
        box on
        end 
    
end

set(gcf,'position',[680 533 865 444]);