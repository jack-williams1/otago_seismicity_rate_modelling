%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Fault-specfic MFD Plots from 2022 NSHM Inversion for %%%%%
%%%%%%%%  ruptures within Otago polygon %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%read in MFD data using inversion results from Solvis that have been filtered by area
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
geo_branch_indx=zeros(num_branches,1); ged_branch_indx=zeros(num_branches,1);
b_weight=[0.24 0.53 0.23]; c=9.05; b=1.5; %scaling between moment and magnitude

for ii=1:num_branches
    
    geo_branch_rupture_rates{ii}=readtable(strjoin([geo_branches(ii),'_Mmod_Rtrim.csv'],""));
    ged_branch_rupture_rates{ii}=readtable(strjoin([ged_branches(ii),'_Mmod_Rtrim.csv'],""));

end

%read in geologic and geodetic dm weighted mean rates and section info
%rates created in solvis/scripts/combine_branches.py
geo_mean_rates=readtable('geo_rates_w.csv'); ged_mean_rates=readtable('ged_rates_w.csv');

geo_sec_info=readtable('fault_sections_geo_w.csv'); ged_sec_info=readtable('fault_sections_ged_w.csv');

addpath([mydir(1:idcs(end-1)-1),'/catalogs_stochastic/model2']); load('orb_fault_parameters','dm');
%read in section area from  solvis/scripts/get_section_area.py
sec_area=readtable('otago_section_area.csv');

%check is consistent with catalog used for RSQSim
%(except for Titri Combined vs. Titri North and Central)
otago_faults=readtable('otago_fault_list_20240428.csv'); num_fault=height(otago_faults);
nshm_fault_stats=zeros(num_fault,5);


%% Find MFD information for each fault

table_headers=["rup_indx","rup_mw","rup_weighted_mw","num_fault","bN_0.823-2.7",...
    "bN_0.959-3.4","bN_1.089-4.6","weighted_mean"];

for ii=1:2 %loop between geologic and then geodetic inversion solutions

    %section info created in solvis/scripts/get_rupture_traces_unscaled.py
    if ii==1
        fault_sec_info=readtable('fault_sections_geo_w.csv');
        branch_rupture_rates= geo_branch_rupture_rates;
        mean_rates=geo_mean_rates;
        folder_dir='by_fault_geo';
    else
        fault_sec_info=readtable('fault_sections_ged_w.csv');
        branch_rupture_rates= ged_branch_rupture_rates;
        mean_rates=ged_mean_rates;
        folder_dir='by_fault_ged';
    end

    fault_list=unique(fault_sec_info.ParentName);
    
    for jj=1:length(fault_list) %for each fault

        b_count=0;

        if isempty(find(string(fault_list(jj))==string(otago_faults.name)))==0 %check is Otago Fault
        
        %index fault within otago_faults.name
        count=find(string(fault_list(jj))==string(otago_faults.name));

        fault_indx=find(string(fault_sec_info.ParentName)==string(fault_list(jj)));
        %find all ruptures that fault participates in
        rup_indx=unique(fault_sec_info.RuptureIndex(fault_indx));
        
        fault_rup_rate=zeros(length(rup_indx),8);
        fault_rup_rate(:,1)=rup_indx; 
        
        %find magnitude, and area-weighted magnitude, for each rupture 
        for kk=1:length(rup_indx)
            %find times rupture is indexed for each fault
            tmp_indx=find(rup_indx(kk)==fault_sec_info.RuptureIndex(fault_indx));
             %get section ID's of faults involved in the rupture
             flt_sections=fault_sec_info.section(fault_indx(tmp_indx));
             flt_sec_area=zeros(length(flt_sections),1);
            
             for mm=1:length(flt_sections)
                 flt_sec_area(mm)=sec_area.Var2(find(flt_sections(mm)==sec_area.Var1));
             end
            
             %get rupture's moment
             rup_mo=10^(fault_sec_info.Magnitude(fault_indx(tmp_indx(1)))*b+c);
             %area weighted magnitude of fault
             fault_rup_rate(kk,3)=(log10(rup_mo*sum(flt_sec_area)/fault_sec_info.Area_m_2_(fault_indx(tmp_indx(1))))-c)/b;
             %find ruptures total magnitude
             fault_rup_rate(kk,2)=fault_sec_info.Magnitude(fault_indx(tmp_indx(1))); 
            
             %find number of participating faults in each rupture
              parent_indx=find(rup_indx(kk)==fault_sec_info.RuptureIndex);
              fault_rup_rate(kk,4)=length(unique(fault_sec_info.ParentID(parent_indx)));
            end
        
        %find rupture rate from different logic tree branches
        for kk=1:num_branches+1
            
            if kk<=3
                 b_count=b_count+1; 
                 rup_rate=branch_rupture_rates{b_count}.AnnualRate;%rupture rates
                 rup_id=branch_rupture_rates{b_count}.RuptureIndex;%rupture indxes
            
            elseif kk ==4%for weighted solution
                rup_rate=mean_rates.w_rate;%rupture rates
                rup_id=mean_rates.rup_id;%rupture indxes
                 
            end
              
                for mm=1:length(rup_indx)
                    %find ruptures associated with fault in each inversion result
                    tmp_indx=find(rup_indx(mm)==rup_id);
                    
                    if isempty(tmp_indx)==0
                        fault_rup_rate(mm,kk+4)=rup_rate(tmp_indx);
                    end
                end 

        end %end kk loop for branches
        
    %derive fault mo rate in inversion for mean geologic logic tree branch
    nshm_fault_stats(count,1)=sum(10.^(fault_rup_rate(:,3).*b+c).*fault_rup_rate(:,8));
  
    %dervive proportion of multifault events for selected logic tree branch
    multi_fault_indx=find(fault_rup_rate(:,4)>1);
    nshm_fault_stats(count,2)=sum(fault_rup_rate(multi_fault_indx,8))/sum(fault_rup_rate(:,8));
    
    %weighted mean magnitude of ruptures that fault participates in
    nshm_fault_stats(count,3)=sum(fault_rup_rate(:,2).*fault_rup_rate(:,8))/sum(fault_rup_rate(:,8));
    
    %derive min and max of ruptures that fault participates in
    nshm_fault_stats(count,4)=min(fault_rup_rate(:,2));
    nshm_fault_stats(count,5)=max(fault_rup_rate(:,2));
    fault_table=array2table(fault_rup_rate,'VariableNames',table_headers);
    writetable(fault_table, [folder_dir,'/',char(fault_list(jj)),'.csv']);  
    
   end %end if statement to only select Otago faults
   end %end jj loop for each fault


save(strcat(folder_dir,'/nshm_fault_stats'),'nshm_fault_stats');

clear nshm_fault_stats

end %end ii loop for geol vs geod rates

%% Plot Fault MFD

%select fault
fault_select='Akatore'; 

mag_range=[6.5:0.05:8.2]; 
mdf_rate1=zeros(length(mag_range),2); mdf_rate2=zeros(length(mag_range),2);

col_opt1=vertcat([0 0 1],[1 0 1]); p={};
col_opt2=vertcat([0 0 1 0.5],[1 0 1 0.5]); 

figure(1);

for ii=1:2

    if ii==1
        folder_dir='by_fault_geo';
    else
        folder_dir='by_fault_ged';
    end

    tmp_table=readtable([folder_dir,'/',fault_select,'.csv']);

    %Derive rates for each magnitude
    for kk=1:length(mag_range)
      
     %MFD rate for full rupture magntitudes
     tmp_indx1=find(tmp_table.rup_mw>mag_range(kk));
     
     mdf_rate1(kk,ii)=sum(tmp_table.weighted_mean(tmp_indx1));
     
     %MFD rate for partial magnitudes
     tmp_indx2=find(tmp_table.rup_weighted_mw>mag_range(kk));
     mdf_rate2(kk,ii)=sum(tmp_table.weighted_mean(tmp_indx2));

    end

    p{ii}=semilogy(mag_range,mdf_rate1(:,ii),'Color',col_opt1(ii,:),'LineWidth',1.5);hold on  
        semilogy(mag_range,mdf_rate2(:,ii),'Color',col_opt2(ii,:),'LineWidth',1.5);hold on   

end

legend_txt=["weighed geologic mean","weighted geodetic mean"]; 

legend([p{1} p{2}],legend_txt,'Location','Northeast');
axis([min(mag_range) max(mag_range) 2*10^-6 5*10^-4]); axis square

xlabel('Magnitude'); ylabel('Annual frequency'); grid on;

%% Plot multiple fault MFD [for a given inversion deformation model]

fault_select=["Titri Central","Titri South","Titri North"]; %select fault to plot MFD
dm_select=1; %1 = geologic, 2 = geodetic

if dm_select==1
    folder_dir='by_fault_geo';
elseif dm_select==2
   folder_dir='by_fault_ged';
end

figure(2);

%needs to match length(fault_select)
col_opt1=vertcat([0 0 1],[1 0 1],[1 0.5 0]);
col_opt2=vertcat([0 0 1 0.5],[1 0 1 0.5],[1 0.5 0 0.5]); 

mag_range=[6.5:0.05:8.2];
mdf_rate1=zeros(length(mag_range),1); mdf_rate2=zeros(length(mag_range),1);

for ii=1:length(fault_select)
    
    tmp_table=readtable(strjoin([folder_dir,'/',fault_select(ii),'.csv'],''));

    for kk=1:length(mag_range)
      
       %MFD rate for full rupture magntitudes (geologic mean only)
       tmp_indx1=find(tmp_table.rup_mw>mag_range(kk));
       mdf_rate1(kk,:)=sum(tmp_table.weighted_mean(tmp_indx1));
     
       %MFD rate for partial magnitudes (geologic mean only)
       tmp_indx2=find(tmp_table.rup_weighted_mw>mag_range(kk));
       mdf_rate2(kk,:)=sum(tmp_table.weighted_mean(tmp_indx2));

    end
       p{ii}=semilogy(mag_range,mdf_rate1,'Color',col_opt1(ii,:),'LineWidth',1.5);hold on  
        semilogy(mag_range,mdf_rate2,'Color',col_opt2(ii,:),'LineWidth',1.5);hold on   
end

legend([p{1} p{2} p{3}],fault_select,'Location','Northeast'); grid on;
axis([min(mag_range)+0.3 max(mag_range)-0.2 1*10^-6 2*10^-4]); axis square

xlabel('Magnitude'); ylabel('Annual frequency'); 

%% Plot mean geologic RI and slip rate by fault

weighted_mean_ri=zeros(num_fault,1); nzcfm_sr=zeros(num_fault,1);

%load NZ CFM info about faults
orb_faults=readtable([mydir(1:idcs(end-1)-1),'/OtagoRangeBasinFaults.xlsx'],'Sheet','segmented'); 
nzcfm_fault_id=orb_faults.Fault_ID;

for ii=1:height(otago_faults)
    
    %read NZ CFM fault slip rate
    tmp_indx=find(otago_faults.Fault_ID(ii)==orb_faults.Fault_ID);
    nzcfm_sr(ii)=orb_faults.SR_pref(tmp_indx);

    %index rupture rates for each fault
    tmp_table=readtable(strcat('by_fault_geo/',char(otago_faults.name(ii)),'.csv'));
    %sum rates from all ruptures that fault partipates in
    weighted_mean_ri(ii,1)=sum(tmp_table.weighted_mean); %weighted mean ri for assessed geologic logic tree branches

    
end

figure(3);

%paleoseismic constraints for Otago used in the IFM from Tables 2&4 in Coffey et al (2024)
%1) Akatore (RVD method)
%2) Dunstan (Timings method)
%3) Titri (RVD method)

ri_constraints=vertcat([11500,7230,77200],[4520,2560,15500],[10530,5330,38900]);
ri_constraints_indx=[1,9,43];

yyaxis right

%plt
p1=bar([1:1:num_fault],nzcfm_sr);
p1.FaceAlpha=0.5; p1.FaceColor=[1 0 0];

ylabel('NZ CFM slip rate (mm/yr)'); ylim([0 1.5]); set(gca,'YColor','k');
yyaxis left

p2=errorbar(ri_constraints_indx,log10(1./ri_constraints(:,1)),log10(1./ri_constraints(:,1))-log10(1./ri_constraints(:,2)),log10(1./ri_constraints(:,3))-log10(1./ri_constraints(:,1)),...
     "o","MarkerFaceColor",[1 1 1],'LineWidth',1.3,'Color',[0 0 1],"MarkerEdgeColor",[0 0 1]);hold on
%
p3=plot([1:1:num_fault],log10(weighted_mean_ri(:,1)),'kx','LineWidth',1.2,"MarkerSize",8.5);hold on
xticks([1:1:num_fault]); xticklabels(otago_faults.name); xtickangle(60);

ylabel('log_{10} participation rate'); set(gca,'YColor','k');
legend([p1 p2 p3],{'NZ CFM slip rate','paleoseismic constraint','IFM weighted geologic'},'location','eastoutside');

set(gcf,'Position',[440 277 843 421]);
