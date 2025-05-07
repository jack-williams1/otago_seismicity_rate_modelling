%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  MFD Plots from 2022 NSHM Inversion for %%%%%
%%%%%%%%  ruptures within Otago polygon %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%read in MFD data using inversion results from Solvis that have been filtered by orb_area_polygons.geojson
mydir  = pwd; idcs   = strfind(mydir,'/');

addpath([mydir(1:idcs(end)-1),'/solvis/OtagoFaults']);

%branches explored, these must match file names in logic_tree_branches
%All branches are time independent and non-stationarity =1

geo_branches=["SW52ZXJzaW9uU29sdXRpb246MTEzMDUz",... %dm: geologic, n-b pair: 0.823-2.7
        "SW52ZXJzaW9uU29sdXRpb246MTEzMDMy",... %dm: geologic, n-b pair: 0.959-3.4
        "SW52ZXJzaW9uU29sdXRpb246MTEzMDM5"];%dm: geologic, n-b pair: 1.089-4.6

ged_branches=["SW52ZXJzaW9uU29sdXRpb246MTEzMDY3",... %dm geodetic, n-b pair: 0.823-2.7
    "SW52ZXJzaW9uU29sdXRpb246MTEzMDgw",...%dm geodetic, n-b pair: 0.959-3.4
    "SW52ZXJzaW9uU29sdXRpb246MTEzMDYz"]; %dm geodetic, n-b pair: 1.089-4.6

num_branches=length(geo_branches); %number of branches assessing 

for ii=1:num_branches
    
    geo_branch_rupture_rates{ii}=readtable(strjoin([geo_branches(ii),'_Mmod_Rtrim.csv'],""));
    ged_branch_rupture_rates{ii}=readtable(strjoin([ged_branches(ii),'_Mmod_Rtrim.csv'],""));

end

%read in geologic and geodetic dm weighted mean rates
%rates created in solvis/scripts/combine_branches.py
geo_mean_rates=readtable('geo_rates_w.csv'); ged_mean_rates=readtable('ged_rates_w.csv'); 
addpath([mydir(1:idcs(end-1)-1),'/catalogs_stochastic/model2']); load('orb_fault_parameters','dm');

%This is NOT the same minimum modelled magnitude used in the inversion (M 6.9, Gerstenberger et al 2024)
%It is less to reflect that the seismic moment of some ruptures has been scaled down
Mmin_nshm2022=6.5; 
Mmax_nshm2022=round(table2array(max(geo_mean_rates(:,2))),2);
mag_range_nshm2022=Mmin_nshm2022:dm:Mmax_nshm2022; 
c=9.05; b=1.5; %scaling between moment and magnitude

%add in other info
addpath([mydir(1:idcs(end-1)-1),'/nshm_urz']); load nshm_urz_analysis

%% Calculate magnitude exceedance rate based from geologic and geodetic weighted mean rate of each rupture

%cells for results from all logic tree branches
geo_w_mag_rate=zeros(length(mag_range_nshm2022),2); ged_w_mag_rate=zeros(length(mag_range_nshm2022),2);
geo_lb_mag_rate=zeros(length(mag_range_nshm2022),num_branches); ged_lb_mag_rate=zeros(length(mag_range_nshm2022),num_branches);

%logic tree branch weightings for different n-b pairs (Gerstenberger et al 2024)
b_weight=[0.24 0.53 0.23];

for ii=1:length(mag_range_nshm2022)
    
    %M>m rate for entire ruptures
    tmp_indx1=find(geo_mean_rates.mag>mag_range_nshm2022(ii));
    geo_w_mag_rate(ii,1)=sum(geo_mean_rates.w_rate(tmp_indx1));
    
    tmp_indx1=find(ged_mean_rates.mag>mag_range_nshm2022(ii));
    ged_w_mag_rate(ii,1)=sum(ged_mean_rates.w_rate(tmp_indx1));
    
    %M>m rate for ruptures when scaled to whether they intersect with Otago
    tmp_indx2=find(geo_mean_rates.w_mag>mag_range_nshm2022(ii));
    geo_w_mag_rate(ii,2)=sum(geo_mean_rates.w_rate(tmp_indx2));

    tmp_indx2=find(ged_mean_rates.w_mag>mag_range_nshm2022(ii));
    ged_w_mag_rate(ii,2)=sum(ged_mean_rates.w_rate(tmp_indx2));
    
    %Repeat for individual logic tree branch solutions
    for jj=1:num_branches
       
       tmp_indx3=find(geo_branch_rupture_rates{jj}.weighted_magnitude>mag_range_nshm2022(ii));
       tmp_indx4=find(ged_branch_rupture_rates{jj}.weighted_magnitude>mag_range_nshm2022(ii));
       geo_lb_mag_rate(ii,jj)=sum(geo_branch_rupture_rates{jj}.AnnualRate(tmp_indx3));
       ged_lb_mag_rate(ii,jj)=sum(ged_branch_rupture_rates{jj}.AnnualRate(tmp_indx4));
    end
    

end

%% Make MFD Plots
% Plot weighted rupture rates with total rupture magntiude, and rupture magnitude in Otago

figure(100);

semilogy(mag_range_nshm2022,geo_w_mag_rate(:,1),'r-','LineWidth',1.2);hold on
semilogy(mag_range_nshm2022,geo_w_mag_rate(:,2),'k-','LineWidth',1.2);

legend('All ruptures','Scaled by area within Otago')
axis([Mmin_nshm2022 Mmax_nshm2022 5*10^-6 5*10^-3]);

xlabel('Magnitude');ylabel('Annual rate of exceedance'); axis square; grid on

% MFD plots for all logic tree branches
figure(101);

col_opt=vertcat([0 0 1],[1 0 1],[1,0.5,0],[0.2 0.5 0.2]);
tiledlayout(1,2,'TileSpacing','compact');

title_opt=["(a)","(b)"]; txt_opt=["geol weighted mean","geod weighted mean"];
tmp_width=[1 1 1 1.5];

for hh=1:2
    
    nexttile
    
    if hh==1 %plot geologic rates
        semilogy(mag_range_nshm2022,geo_w_mag_rate(:,2),'k-','LineWidth',1.5); hold on
        tmp_rates=geo_lb_mag_rate;  
    else    %plot geodetic rates
        semilogy(mag_range_nshm2022,ged_w_mag_rate(:,2),'k-','LineWidth',1.5); hold on
        tmp_rates=ged_lb_mag_rate;  
    end

    for ii=1:num_branches
    
        semilogy(mag_range_nshm2022,tmp_rates(:,ii),'Color',col_opt(ii,:),'LineWidth',tmp_width(ii)); hold on
    
    end

    legend([txt_opt(hh),{"n-b: 2.7-0.823","n-b: 3.4-0.959","n-b: 4.6-1.089"}],'Location','Northeast');

    axis([Mmin_nshm2022 Mmax_nshm2022 5*10^-6 7*10^-3]);

    xlabel('Magnitude');ylabel('Annual rate of exceedance'); axis square; grid on

    t2=title(title_opt(hh),'fontsize',12,'fontweight','normal');
    set(t2, 'horizontalAlignment', 'left');
    set(t2, 'units', 'normalized');
    h2 = get(t2, 'position');
    set(t2, 'position', [-0.15 h2(2) h2(3)]);

end

set(gcf, 'Position',[218 324 936 394])

%% Assess proportion of multifault ruptures and moment rate

%read in information about each rupture's section
%section info created in solvis/scripts/get_rupture_traces_unscaled.py

ifm_otago_fault_stats=cell(3,1);

for hh=1:3

    if hh==1 %weighted mean geologic rates
        fault_sec_info=readtable('fault_sections_geo_w.csv');
        rup_indx=geo_mean_rates.rup_id; %read in ruptures from weighted inversion results
        rup_rate=geo_mean_rates.w_rate;%rupture rates 
        rup_wmag=geo_mean_rates.w_mag;

    elseif hh==2 %weighted mean geodetic rates
        fault_sec_info=readtable('fault_sections_ged_w.csv');
        rup_indx=ged_mean_rates.rup_id; %read in ruptures from weighted inversion results
        rup_rate=ged_mean_rates.w_rate;%rupture rates 
        rup_wmag=ged_mean_rates.w_mag;
    end

    if hh~=3

        num_fault=zeros(length(rup_indx),1); num_fault_weighted=zeros(length(rup_indx),1);     
        %loop through each rupture

        for jj=1:length(rup_indx) 
            rup_sect_indx=find(rup_indx(jj)==fault_sec_info.RuptureIndex);%index rupture in sections table
            rup_parent_faults=fault_sec_info.ParentID(rup_sect_indx);%participant parent faults of each rupture
            num_fault(jj)=length(unique(rup_parent_faults));%number of faults in each rupture
        end
            
        %proportion of ruptures that are multi faults
        ifm_otago_fault_stats{hh}(1)=length(find(num_fault>1))/length(rup_indx);
        %proportion of ruptures that are multi faults once weighted by rate
        ifm_otago_fault_stats{hh}(2)=sum(rup_rate(find(num_fault>1)))/sum(rup_rate);
            
        w_rupture_mo=10.^((b.*rup_wmag)+c); %seismic moment of each rupture
        w_rupture_mo_r=w_rupture_mo.*rup_rate; %moment rate of each rupture
        ifm_otago_fault_stats{hh}(3)=sum(w_rupture_mo_r);%total moment rate of all ruptures

    end

    if hh==3 %geologic logic tree branches

        fault_sec_info=readtable('fault_sections_geo_w.csv');

        for ii=1:num_branches
            rup_indx=geo_branch_rupture_rates{ii}.RuptureIndex;%read in ruptures from geologic logic tree branches
            rup_rate=geo_branch_rupture_rates{ii}.AnnualRate;%rupture rates
            rup_wmag=geo_branch_rupture_rates{ii}.weighted_magnitude;
    
            num_fault=zeros(length(rup_indx),1); num_fault_weighted=zeros(length(rup_indx),1); 

            %for each rupture
            for jj=1:length(rup_indx) 
        
                rup_sect_indx=find(rup_indx(jj)==fault_sec_info.RuptureIndex);%index rupture in sections table
                rup_parent_faults=fault_sec_info.ParentID(rup_sect_indx);%participant parent faults of each rupture
                num_fault(jj)=length(unique(rup_parent_faults));%number of faults in each rupture
            end
    
            ifm_otago_fault_stats{hh}(ii,1)=length(find(num_fault>1))/length(rup_indx);
            ifm_otago_fault_stats{hh}(ii,2)=sum(rup_rate(find(num_fault>1)))/sum(rup_rate);
        
            %determine logic tree branch moment rate
            rupture_mo=10.^((b.*rup_wmag)+c); 
            rupture_mo_r=rupture_mo.*rup_rate;
            ifm_otago_fault_stats{hh}(ii,3)=sum(rupture_mo_r);

    end %end ii loop for num branches
    
    end %if hh=3 condition

end

%% Combine with URZ

%IFM not a complete description of seismicity

mag_range_nshm2022_combined=[4.95:0.05:Mmax_nshm2022];
nshm2022_combined_rate=zeros(length(mag_range_nshm2022_combined),2);

for ii=1:length(mag_range_nshm2022_combined)
    
    if mag_range_nshm2022_combined(ii)<6.7
       %rate entirely described by the URZ (negative binomial model)
       mag_indx=ii; %this only works when both mag range increments are 0.05 and start at M4.95!
       nshm2022_combined_rate(ii,:)=urz_mag_rate(mag_indx,1);
       
    elseif mag_range_nshm2022_combined(ii)>=6.7 && mag_range_nshm2022_combined(ii)<=8.0  
       %rate  combination of the URZ (negative binomial model) and IFM
       mag_indx1=ii; %this only works when both mag range increments are 0.05 and start at M4.95!
       mag_indx2=find(round(mag_range_nshm2022,3)==round(mag_range_nshm2022_combined(ii),3));
       
       nshm2022_combined_rate(ii,1)=0.8*geo_w_mag_rate(mag_indx2,2)+0.2*urz_mag_rate(mag_indx1,1); %combined with geologic dm mean
       nshm2022_combined_rate(ii,2)=0.8*ged_w_mag_rate(mag_indx2,2)+0.2*urz_mag_rate(mag_indx1,1); %combined with geodetic dm mean
       
    elseif mag_range_nshm2022_combined(ii)>8.0 
       mag_indx2=find(round(mag_range_nshm2022,3)==round(mag_range_nshm2022_combined(ii),3));
       nshm2022_combined_rate(ii,1)=geo_w_mag_rate(mag_indx2,2);
       nshm2022_combined_rate(ii,2)=ged_w_mag_rate(mag_indx2,2);
       
    end
    
end

%Plot MFD, note plots rates only for: (1) weighted logic tree branch, (2)
%weighted mean geologic dm and (3) geodetic dm

figure(102);

p1=semilogy(urz_mag_range,urz_mag_rate(:,1),'LineWidth',1.5,'Color',[1 0 0]);hold on
p2=semilogy(mag_range_nshm2022,geo_w_mag_rate(:,2),'LineWidth',0.7,'Color',[1 0 1 0.7]);hold on
p3=semilogy(mag_range_nshm2022,ged_w_mag_rate(:,2),'LineWidth',0.7,'Color',[0.2 0.5 0.2]);hold on
p4=semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,1),'LineWidth',1.5,'Color',[1 0 1]);hold on
p5=semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,2),'LineWidth',1.5,'Color',[0.2 0.5 0.2]);hold on


legend([p1,p2,p4,p3,p5],{"URZ","geol weighted mean","geol weighted mean-URZ combined",...
        "geod weighted mean","geod weighted mean-URZ combined"},...
        'Location','Northeast');

axis([Mmin_nshm2022 Mmax_nshm2022 10^-5 5*10^-3]);
xlabel('Magnitude');ylabel('Annual rate of exceedance'); axis square; grid on


%% Save analysis

save('nshm_otago_inversion_results','mag_range_nshm2022','geo_w_mag_rate','ged_w_mag_rate','nshm2022_combined_rate','ifm_otago_fault_stats',...
    'mag_range_nshm2022_combined','num_branches','geo_branch_rupture_rates','ged_branch_rupture_rates','b_weight');
