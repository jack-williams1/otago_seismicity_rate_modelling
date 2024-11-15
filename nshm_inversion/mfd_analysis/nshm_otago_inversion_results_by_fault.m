%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Fault-specfic MFD Plots from 2022 NSHM Inversion for %%%%%
%%%%%%%%  ruptures within Otago polygon %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
b_weight=[0.24 0.53 0.23]; c=9.05; b=1.5; %scaling between moment and magnitude

for ii=1:num_branches
    
    geo_branch_rupture_rates{ii}=readtable(strjoin(['otago_ruptures_',geo_branches(ii),'.csv'],""));
    ged_branch_rupture_rates{ii}=readtable(strjoin(['otago_ruptures_',ged_branches(ii),'.csv'],""));
    geo_branch_indx(ii)=find(geo_branches(ii)==logic_tree_branches.branch); %look up index for branch in logic tree branch table
    ged_branch_indx(ii)=find(ged_branches(ii)==logic_tree_branches.branch); %look up index for branch in logic tree branch table
end

%read in section info
fault_sec_info=readtable('otago_fault_sections.csv'); sec_area=readtable('otago_section_area.csv');

%check is consistent with catalog used for RSQSim
%(except for Titri Combined vs. Titri North and Central)
otago_faults=readtable('otago_fault_list_20240428.csv'); num_fault=height(otago_faults);

fault_list=unique(fault_sec_info.ParentName);
nshm_fault_stats=zeros(num_fault,5);


%% Find MFD information for each fault

%select logic tree branch indx. 5 = weighted mean; 9 = weighted geol mean
brnch_indx=5; 

table_headers=["rup_indx","rup_mw","rup_weighted_mw","num_fault","weighted_mean_rate",...
    "geologic_dm_bN_0.823-2.7","geologic_dm_bN_0.959-3.4","geologic_dm_bN_1.089-4.6","weighted_geologic_mean",...
    "geodetic_dm_bN_0.823-2.7","geodetic_dm_bN_0.959-3.4","geodetic_dm_bN_1.089-4.6","weighted_geodetic_mean"];
  
for ii=1:length(fault_list) %for each fault

    b_count1=0; b_count2=0;

    if isempty(find(string(fault_list(ii))==string(otago_faults.name)))==0 %check is Otago Fault
        
        %index fault within otago_faults.name
        count=find(string(fault_list(ii))==string(otago_faults.name));
        
        fault_indx=find(string(fault_sec_info.ParentName)==string(fault_list(ii)));
        %find all ruptures that fault participates in
        rup_indx=unique(fault_sec_info.RuptureIndex(fault_indx));
        
        fault_rup_rate=zeros(length(rup_indx),13);
        fault_rup_rate(:,1)=rup_indx; 
        
        %find magnitude, and fault's area-weighted magnitude, for each rupture 
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
            fault_rup_rate(kk,3)=(log10(rup_mo*sum(flt_sec_area)/fault_sec_info.Area_m_2_(fault_indx(tmp_indx(1))))-9.05)/1.5;
            %find ruptures total magnitude
            fault_rup_rate(kk,2)=fault_sec_info.Magnitude(fault_indx(tmp_indx(1))); 
            
            %find number of particpating faults in each rupture
            parent_indx=find(rup_indx(kk)==fault_sec_info.RuptureIndex);
            fault_rup_rate(kk,4)=length(unique(fault_sec_info.ParentID(parent_indx)));
        end
        
        %find rupture rate from different logic tree branches
        for jj=1:9
            
            if jj==1 %for weighted solution
                rup_rate=orb_ruptures.rate_weighted_mean;%rupture rates
                rup_id=orb_ruptures.RuptureIndex;%rupture indxes

            elseif jj>1 && jj<=4
                 b_count1=b_count1+1; 
                 tmp_table=readtable(strjoin(['otago_ruptures_',geo_branches(b_count1),'.csv'],""));
                 rup_rate=tmp_table.AnnualRate;%rupture rates
                 rup_id=tmp_table.RuptureIndex;%rupture indxes
                 
            elseif jj>=6 && jj<9
                 b_count2=b_count2+1; 
                 tmp_table=readtable(strjoin(['otago_ruptures_',ged_branches(b_count2),'.csv'],""));
                 rup_rate=tmp_table.AnnualRate;%rupture rates
                 rup_id=tmp_table.RuptureIndex;%rupture indxes    
                 
            end
            
            if jj==5 || jj==9
                
                fault_rup_rate(:,jj+4)=fault_rup_rate(:,jj+1).*b_weight(1)+fault_rup_rate(:,jj+2).*b_weight(2)+fault_rup_rate(:,jj+3).*b_weight(3);
            else    
              for kk=1:length(rup_indx)
                    %find ruptures associated with fault in each inversion result
                    tmp_indx=find(rup_indx(kk)==rup_id);
                    
                    if isempty(tmp_indx)==0
                        fault_rup_rate(kk,jj+4)=rup_rate(tmp_indx);
                    end
              end 

            end
        end %end ii loop for branches
        
    %derive fault mo rate in inversion for selected logic tree branch

    nshm_fault_stats(count,1)=sum(10.^(fault_rup_rate(:,3).*b+c).*fault_rup_rate(:,brnch_indx));
  
    %dervive proportion of multifault events for selected logic tree branch
    multi_fault_indx=find(fault_rup_rate(:,4)>1);
    nshm_fault_stats(count,2)=sum(fault_rup_rate(multi_fault_indx,brnch_indx))/sum(fault_rup_rate(:,brnch_indx));
    
    %weighted mean magnitude of ruptures that fault participates in
    nshm_fault_stats(count,3)=sum(fault_rup_rate(:,2).*fault_rup_rate(:,brnch_indx))/sum(fault_rup_rate(:,brnch_indx));
    
    %derive min and max of ruptures that fault participates in
    tmp_indx=find(fault_rup_rate(:,brnch_indx)>0);
    nshm_fault_stats(count,4)=min(fault_rup_rate(tmp_indx,2));
    nshm_fault_stats(count,5)=max(fault_rup_rate(tmp_indx,2));
    fault_table=array2table(fault_rup_rate,'VariableNames',table_headers);
    writetable(fault_table, ['by_fault/',char(fault_list(ii)),'.csv']);  
    
    end %end if statement to only select Otago faults
end

save(strcat('by_fault/nshm_fault_stats',num2str(brnch_indx)),'nshm_fault_stats');

clear nshm_fault_stats

%% Plot Fault MFD for multiple branches

fault_select='Akatore'; %select fault to plot MFD

tmp_table=readtable(['by_fault/',fault_select,'.csv']);

mag_range=[6.5:0.05:8.2]; tmp_indx1={}; tmp_indx2={}; 
mdf_rate1=zeros(length(mag_range),3); mdf_rate2=zeros(length(mag_range),3);

%Derive rates for each magnitude
for kk=1:length(mag_range)
      
     %MFD rate for full rupture magntitudes
     tmp_indx1{kk}=find(tmp_table.rup_mw>mag_range(kk));
     
     mdf_rate1(kk,:)=[sum(tmp_table.weighted_mean_rate(tmp_indx1{kk})),sum(tmp_table.weighted_geologic_mean(tmp_indx1{kk})),...
            sum(tmp_table.weighted_geodetic_mean(tmp_indx1{kk}))];
     
     %MFD rate for partial magnitudes
     tmp_indx2{kk}=find(tmp_table.rup_weighted_mw>mag_range(kk));
     mdf_rate2(kk,:)=[sum(tmp_table.weighted_mean_rate(tmp_indx2{kk})),sum(tmp_table.weighted_geologic_mean(tmp_indx2{kk})),...
         sum(tmp_table.weighted_geodetic_mean(tmp_indx2{kk}))];     
end

%make MFD plot
figure(1);

legend_txt=["weighted mean","weighed geologic mean","weighted geodetic mean"]; p={};

col_opt1=vertcat([0 0 1],[1 0 1],[1 0.5 0]);
col_opt2=vertcat([0 0 1 0.5],[1 0 1 0.5],[1 0.5 0 0.5]); 

for mm=1:3
    
        p{mm}=semilogy(mag_range,mdf_rate1(:,mm),'Color',col_opt1(mm,:),'LineWidth',1.5);hold on  
        semilogy(mag_range,mdf_rate2(:,mm),'Color',col_opt2(mm,:),'LineWidth',1.5);hold on   
 
end

legend([p{1} p{2} p{3}],legend_txt,'Location','Northeast'); grid on;
axis([min(mag_range) max(mag_range) 2*10^-6 5*10^-4]); axis square

xlabel('Magnitude'); ylabel('Annual frequency'); 

%% Plot multiple fault MFDs for single branch (geolgic_dm currently selected)

fault_select=["Titri Central","Titri South","Titri North"]; %select fault to plot MFD

figure(2);

%col_opt length must match fault_select length
col_opt1=vertcat([0 0 1],[1 0 1],[1 0.5 0]);
col_opt2=vertcat([0 0 1 0.5],[1 0 1 0.5],[1 0.5 0 0.5]); 

mag_range=[6.5:0.05:8.2]; tmp_indx1={}; tmp_indx2={}; p={};
mdf_rate1=zeros(length(mag_range),1);  mdf_rate2=zeros(length(mag_range),1);

for ii=1:length(fault_select)
    
    tmp_table=readtable(strjoin(['by_fault/',fault_select(ii),'.csv'],''));

    for kk=1:length(mag_range)
      
     %MFD rate for full rupture magntitudes (geologic mean only)
       tmp_indx1{kk}=find(tmp_table.rup_mw>mag_range(kk));
        mdf_rate1(kk)=sum(tmp_table.weighted_geologic_mean(tmp_indx1{kk}));

     
     %MFD rate for partial magnitudes (geologic mean only)
     tmp_indx2{kk}=find(tmp_table.rup_weighted_mw>mag_range(kk));
      mdf_rate2(kk)=sum(tmp_table.weighted_geologic_mean(tmp_indx2{kk}));
    end

       p{ii}=semilogy(mag_range,mdf_rate1,'Color',col_opt1(ii,:),'LineWidth',1.2);hold on  
        semilogy(mag_range,mdf_rate2,'Color',col_opt2(ii,:),'LineWidth',1.5);hold on   
end


legend([p{1} p{2} p{3}],fault_select,'Location','Northeast'); grid on;
axis([min(mag_range)+0.3 max(mag_range)-0.2 1*10^-6 2*10^-4]); axis square

xlabel('Magnitude'); ylabel('Annual frequency'); 

%% Plot mean RI and slip rate by fault

weighted_mean_ri=zeros(num_fault,2); nzcfm_sr=zeros(num_fault,1);

%load NZ CFM info about faults
orb_faults=readtable([mydir(1:idcs(end-1)-1),'/OtagoRangeBasinFaults.xlsx'],'Sheet','segmented'); 
nzcfm_fault_id=orb_faults.Fault_ID;

for ii=1:height(otago_faults)
    
    %read NZ CFM fault slip rate
    tmp_indx=find(otago_faults.Fault_ID(ii)==orb_faults.Fault_ID);
    nzcfm_sr(ii)=orb_faults.SR_pref(tmp_indx);

    %index rupture rates for each fault
    tmp_table=readtable(strcat('by_fault/',char(otago_faults.name(ii)),'.csv'));
    %sum rates from all ruptures that fault partipates in
    weighted_mean_ri(ii,1)=sum(tmp_table.weighted_mean_rate); %weighted mean ri for all logic tree branches
    weighted_mean_ri(ii,2)=sum(tmp_table.weighted_geologic_mean); %weighted mean ri for assessed geologic logic tree branches

    
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
p4=plot([1:1:num_fault],log10(weighted_mean_ri(:,2)),'x',"MarkerEdgeColor",[0.4 0.4 0.4],'LineWidth',1.2,"MarkerSize",8.5);hold on
xticks([1:1:num_fault]); xticklabels(otago_faults.name); xtickangle(60);

ylabel('log_{10} participation rate'); set(gca,'YColor','k');
legend([p1 p2 p3 p4],{'NZ CFM slip rate','paleoseismic constraint','IFM weighted all','IFM weighted geologic'},'location','eastoutside');

set(gcf,'Position',[440 277 843 421]);