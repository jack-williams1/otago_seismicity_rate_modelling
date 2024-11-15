%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  MFD Plots from 2022 NSHM Inversion for %%%%%
%%%%%%%%  ruptures within Otago polygon %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%read in MFD data using inversion results from Solvis that have been filtered by area
mydir  = pwd; idcs   = strfind(mydir,'/');

%read in weighted mean rupture rates from all logic tree branches
addpath([mydir(1:idcs(end)-1),'/solvis/WORK']); orb_ruptures=readtable('otago_ruptures_weighted.csv');


%read in results for different logic tree branches
addpath([mydir(1:idcs(end)-1),'/solvis']);
logic_tree_branches=readtable('nshm_logictreebranch_lookuptable.txt','ReadVariableNames', false);
logic_tree_branches=mergevars(logic_tree_branches,["Var4","Var5"]);

logic_tree_branches.Properties.VariableNames={'branch','deformation_model',...
    'time_dependence','bN_pair', 'scaling','Nonstationarity scaling'};

%branches explored, these must match file names in logic_tree_branches
%All branches are time independent and non-stationarity =1

geo_branches=["U2NhbGVkSW52ZXJzaW9uU29sdXRpb246MTIwNzA3",... %dm: geologic, n-b pair: 0.823-2.7
        "U2NhbGVkSW52ZXJzaW9uU29sdXRpb246MTIwNjk4",... %dm: geologic, n-b pair: 0.959-3.4
        "U2NhbGVkSW52ZXJzaW9uU29sdXRpb246MTIwNjc2"];%dm: geologic, n-b pair: 1.089-4.6

ged_branches=["U2NhbGVkSW52ZXJzaW9uU29sdXRpb246MTIwNzE3",... %dm geodetic, n-b pair: 0.823-2.7
    "U2NhbGVkSW52ZXJzaW9uU29sdXRpb246MTIwNzUz",...%dm geodetic, n-b pair: 0.959-3.4
    "U2NhbGVkSW52ZXJzaW9uU29sdXRpb246MTIwNzI4"]; %dm geodetic, n-b pair: 1.089-4.6

num_branches=length(geo_branches); %number of branches assessing 
geo_branch_indx=zeros(num_branches,1); ged_branch_indx=zeros(num_branches,1);

for ii=1:num_branches
    
    geo_branch_rupture_rates{ii}=readtable(strjoin(['otago_ruptures_',geo_branches(ii),'.csv'],""));
    ged_branch_rupture_rates{ii}=readtable(strjoin(['otago_ruptures_',ged_branches(ii),'.csv'],""));
    geo_branch_indx(ii)=find(geo_branches(ii)==logic_tree_branches.branch); %look up index for branch in logic tree branch table
    ged_branch_indx(ii)=find(ged_branches(ii)==logic_tree_branches.branch); %look up index for branch in logic tree branch table
end

%This is NOT the same minimum modelled magnitude used in the inversion (M 6.9, Gerstenberger et al 2024)
%It is less to reflect that the seismic moment of some ruptures has been scaled down
Mmin_nshm2022=6.5; dm=0.005; %dm should be consistent with dm in catalgs_stochastic/fault_recurrence_paramters.m
Mmax_nshm2022=round(max(orb_ruptures.Magnitude),2);
mag_range_nshm2022=Mmin_nshm2022:dm:Mmax_nshm2022; 
c=9.05; b=1.5; %scaling between moment and magnitude

%add in other info
addpath([mydir(1:idcs(end-1)-1),'/nshm_urz']); load nshm_urz_analysis

%% Calculate magnitude exceedance rate based from weighted mean rate of each rupture

%cells for results from all logic tree branches
mag_rate=zeros(length(mag_range_nshm2022),1); w_mag_rate=zeros(length(mag_range_nshm2022),1); 
geo_b_mag_rate=zeros(length(mag_range_nshm2022),num_branches+1); ged_b_mag_rate=zeros(length(mag_range_nshm2022),num_branches+1);

%logic tree branch weightings for different n-b pairs (Gerstenberger et al 2024)
b_weight=[0.24 0.53 0.23];

for ii=1:length(mag_range_nshm2022)
    
    %M>m rate for entire ruptures
    tmp_indx1=find(orb_ruptures.Magnitude>mag_range_nshm2022(ii));
    mag_rate(ii)=sum(orb_ruptures.rate_weighted_mean(tmp_indx1));
    
    %M>m rate for ruptures when scaled to whether they intersect with Otago
    tmp_indx2=find(orb_ruptures.weighted_magnitude>mag_range_nshm2022(ii));
    w_mag_rate(ii)=sum(orb_ruptures.rate_weighted_mean(tmp_indx2));
    
    for jj=1:num_branches
       
       tmp_indx3=find(geo_branch_rupture_rates{jj}.weighted_magnitude>mag_range_nshm2022(ii));
       tmp_indx4=find(ged_branch_rupture_rates{jj}.weighted_magnitude>mag_range_nshm2022(ii));
       geo_b_mag_rate(ii,jj)=sum(geo_branch_rupture_rates{jj}.AnnualRate(tmp_indx3));
       ged_b_mag_rate(ii,jj)=sum(ged_branch_rupture_rates{jj}.AnnualRate(tmp_indx4));
    end
    
    %weighted average for geologic logic tree branches
    geo_b_mag_rate(ii,jj+1)=geo_b_mag_rate(ii,1)*b_weight(1)+geo_b_mag_rate(ii,2)*b_weight(2)+geo_b_mag_rate(ii,3)*b_weight(3);
    ged_b_mag_rate(ii,jj+1)=ged_b_mag_rate(ii,1)*b_weight(1)+ged_b_mag_rate(ii,2)*b_weight(2)+ged_b_mag_rate(ii,3)*b_weight(3);

end

%% Make MFD Plots
% Plot weighted rupture rates with total rupture magntiude, and rupture magnitude in Otago

figure(100);

semilogy(mag_range_nshm2022,mag_rate,'r-','LineWidth',1.2);hold on
semilogy(mag_range_nshm2022,w_mag_rate,'k-','LineWidth',1.2);

legend('All ruptures','Scaled by area within Otago')
axis([Mmin_nshm2022 Mmax_nshm2022 5*10^-6 5*10^-3]);

xlabel('Magnitude');ylabel('Annual rate of exceedance'); axis square; grid on

figure(101);

col_opt=vertcat([0 0 1],[1 0 1],[1,0.5,0],[0.2 0.5 0.2]);
tiledlayout(1,2,'TileSpacing','compact');

title_opt=["(a)","(b)"]; txt_opt=["geol weighted mean","geod weighted mean"];
tmp_width=[1 1 1 1.5];
for hh=1:2
    
    nexttile

    % MFD plots for all logic tree branches

    semilogy(mag_range_nshm2022,w_mag_rate,'k-','LineWidth',1.5); hold on
    tmp_txt=cell(num_branches,1);

    legend_txt=["weighted mean rupture rates"];
    
    if hh==1 %plot geologic rates
        tmp_rates=geo_b_mag_rate; tmp_branch_indx=geo_branch_indx; tmp_branch_rupture_rates=geo_branch_rupture_rates;
    else    %plot geodetic rates
        tmp_rates=ged_b_mag_rate; tmp_branch_indx=ged_branch_indx; tmp_branch_rupture_rates=ged_branch_rupture_rates;
    end

    for ii=1:num_branches+1
    
        semilogy(mag_range_nshm2022,tmp_rates(:,ii),'Color',col_opt(ii,:),'LineWidth',tmp_width(ii)); hold on
        %obtain deformation model, B-n pair, and weighting for each branch
        if ii<=num_branches
            tmp_txt{ii}=strjoin({cell2mat(logic_tree_branches.deformation_model(tmp_branch_indx(ii),:)),...
                cell2mat(logic_tree_branches.bN_pair(tmp_branch_indx(ii),:)),num2str(tmp_branch_rupture_rates{ii}.weight(1))},' ');
        else
            tmp_txt{ii}=txt_opt(hh);
        end

    legend_txt=[legend_txt,tmp_txt{ii}];
    
    end

    legend(legend_txt,'Location','Northeast');

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
fault_sec_info=readtable('otago_fault_sections.csv');

ifm_otago_fault_stats=cell(3,1);

for hh=1:3

    if hh==1
            rup_indx=orb_ruptures.RuptureIndex; %read in ruptures from weighted inversion results
            rup_rate=orb_ruptures.rate_weighted_mean;%rupture rates
            num_fault=zeros(length(rup_indx),1); num_fault_weighted=zeros(length(rup_indx),1); 

            %for each rupture
            for jj=1:length(rup_indx) 
        
                rup_sect_indx=find(rup_indx(jj)==fault_sec_info.RuptureIndex);%index rupture in sections table
                rup_parent_faults=fault_sec_info.ParentID(rup_sect_indx);%participant parent faults of each rupture
                num_fault(jj)=length(unique(rup_parent_faults));%number of faults in each rupture
            end
            
            %proportion of ruptures that are multi faults
            ifm_otago_fault_stats{hh}(1)=length(find(num_fault>1))/length(rup_indx);
            %proportion of ruptures that are multi faults once weighted by rate
            ifm_otago_fault_stats{hh}(2)=sum(rup_rate(find(num_fault>1)))/sum(rup_rate);
            
            w_rupture_mo=10.^((b.*orb_ruptures.weighted_magnitude)+c); %seismic moment of each rupture
            w_rupture_mo_r=w_rupture_mo.*orb_ruptures.rate_weighted_mean; %moment rate of each rupture
            ifm_otago_fault_stats{hh}(3)=sum(w_rupture_mo_r);%total moment rate of all ruptures
    else
    
    for ii=1:num_branches

        if hh==2
            rup_indx=geo_branch_rupture_rates{ii}.RuptureIndex;%read in ruptures from geologic logic tree branches
            rup_rate=geo_branch_rupture_rates{ii}.AnnualRate;%rupture rates
            rup_wmag=geo_branch_rupture_rates{ii}.weighted_magnitude;
        elseif hh==3
            rup_indx=ged_branch_rupture_rates{ii}.RuptureIndex;%read in ruptures from geodetic logic tree branches
            rup_rate=ged_branch_rupture_rates{ii}.AnnualRate;%rupture rates
            rup_wmag=ged_branch_rupture_rates{ii}.weighted_magnitude;
        end
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

    %weighted mean results for different logic tree branches
    ifm_otago_fault_stats{hh}(4,1)=ifm_otago_fault_stats{hh}(1,1)*b_weight(1)+ifm_otago_fault_stats{hh}(2,1)*b_weight(2)+ifm_otago_fault_stats{hh}(3,1)*b_weight(3);
    ifm_otago_fault_stats{hh}(4,2)=ifm_otago_fault_stats{hh}(1,2)*b_weight(1)+ifm_otago_fault_stats{hh}(2,2)*b_weight(2)+ifm_otago_fault_stats{hh}(3,2)*b_weight(3);
    ifm_otago_fault_stats{hh}(4,3)=ifm_otago_fault_stats{hh}(1,3)*b_weight(1)+ifm_otago_fault_stats{hh}(2,3)*b_weight(2)+ifm_otago_fault_stats{hh}(3,3)*b_weight(3);
    
    end
end

%% Combine with URZ

%IFM not a complete description of seismicity

mag_range_nshm2022_combined=[4.95:0.05:Mmax_nshm2022];
nshm2022_combined_rate=zeros(length(mag_range_nshm2022_combined),3);

for ii=1:length(mag_range_nshm2022_combined)
    
    if mag_range_nshm2022_combined(ii)<6.9
       %rate entirely described by the URZ (negative binomial model)
       mag_indx=ii; %this only works when both mag range increments are 0.05 and start at M4.95!
       nshm2022_combined_rate(ii,:)=urz_mag_rate(mag_indx,1);
       
    elseif mag_range_nshm2022_combined(ii)>=6.9 && mag_range_nshm2022_combined(ii)<=8.0  
       %rate  combination of the URZ (negative binomial model) and IFM
       mag_indx1=ii; %this only works when both mag range increments are 0.05 and start at M4.95!
       mag_indx2=find(round(mag_range_nshm2022,3)==round(mag_range_nshm2022_combined(ii),3));
       
       nshm2022_combined_rate(ii,1)=0.8*w_mag_rate(mag_indx2)+0.2*urz_mag_rate(mag_indx1,1); %combined with IFM weighted mean
       nshm2022_combined_rate(ii,2)=0.8*geo_b_mag_rate(mag_indx2,4)+0.2*urz_mag_rate(mag_indx1,1); %combined with geologic dm mean
       nshm2022_combined_rate(ii,3)=0.8*ged_b_mag_rate(mag_indx2,4)+0.2*urz_mag_rate(mag_indx1,1); %combined with geodetic dm mean
       
    elseif mag_range_nshm2022_combined(ii)>8.0 
       mag_indx2=find(round(mag_range_nshm2022,3)==round(mag_range_nshm2022_combined(ii),3));
       nshm2022_combined_rate(ii,1)=0.8*w_mag_rate(mag_indx2);
       nshm2022_combined_rate(ii,2)=0.8*geo_b_mag_rate(mag_indx2,4);
       nshm2022_combined_rate(ii,3)=0.8*ged_b_mag_rate(mag_indx2,4);
       

    end
    
end

%Plot MFD, note plots rates only for: (1) weighted logic tree branch, (2)
%weighted mean geologic dm and (3) geodetic dm

figure(102);

p1=semilogy(urz_mag_range,urz_mag_rate(:,1),'LineWidth',1.5,'Color',[1 0 0]);hold on
p2=semilogy(mag_range_nshm2022,w_mag_rate,'LineWidth',0.7,'Color',[0 0 0 0.7]);hold on
p3=semilogy(mag_range_nshm2022,geo_b_mag_rate(:,4),'LineWidth',0.7,'Color',[1 0 1 0.7]);hold on
p4=semilogy(mag_range_nshm2022,ged_b_mag_rate(:,4),'LineWidth',0.7,'Color',[0.2 0.5 0.2]);hold on
p5=semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,1),'k-','LineWidth',1.5);hold on
p6=semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,2),'LineWidth',1.5,'Color',[1 0 1]);hold on
p7=semilogy(mag_range_nshm2022_combined,nshm2022_combined_rate(:,3),'LineWidth',1.5,'Color',[0.2 0.5 0.2]);hold on


legend([p1,p2,p5,p3,p6,p4,p7],{"URZ","weighted mean","weighted mean-URZ combined",...
        "geol weighted mean","geol weighted mean-URZ combined",...
        "ged weighted mean","ged weighted mean-URZ combined"},...
        'Location','Northeast');

axis([Mmin_nshm2022 Mmax_nshm2022 10^-5 5*10^-3]);
xlabel('Magnitude');ylabel('Annual rate of exceedance'); axis square; grid on


%% Save analysis

save('nshm_otago_inversion_results','mag_range_nshm2022','mag_rate','geo_b_mag_rate','ged_b_mag_rate','w_mag_rate','nshm2022_combined_rate','ifm_otago_fault_stats',...
    'mag_range_nshm2022_combined','num_branches','logic_tree_branches','geo_branch_indx','ged_branch_indx' ,'geo_branch_rupture_rates','ged_branch_rupture_rates','b_weight');