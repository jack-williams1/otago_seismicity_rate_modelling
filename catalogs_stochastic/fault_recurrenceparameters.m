%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Otago Range and Basin Fault Modelling Recurrence Parameters %%
%%%%%%% Considers: (1) NZCFM sources char (segmented), (2) combined NZCFM %%%
%%%%%%%   sources char, and (3) combined NZCFM sources G-R    %%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%Extract data from the NZCFM for individaul faults. This is the 'segmented' case
orb_faults=readtable('OtagoRangeBasinFaults.xlsx','Sheet','segmented');
num_fault=length(orb_faults.Name1);

num_simu=1*10^6;%length of catalog (years)
t_limit=1;%interval of catalog assessed

%Select recurrence model for Otago Faults:
%1) NZCFM slip rates - Aperiodicity = 0.5
%2) NZCFM slip rates - Aperiodicity = 2
%3) NZCFM slip rates - Aperiodicity = 4
%4) Geodetic slip rates - Aperiodicity = 2

recurrence_model=4;

if recurrence_model==1

    orb_sliprates=[orb_faults.SR_pref];
    aperiodicity=0.5;

elseif recurrence_model==2

    orb_sliprates=[orb_faults.SR_pref];
    aperiodicity=2;

elseif recurrence_model==3

    orb_sliprates=[orb_faults.SR_pref];
    aperiodicity=3;  

elseif recurrence_model==4
    %import geodetic slip rates. Need to run
    %nshm_geodetic_rates/extract_otago_sections.py first
    mydir  = pwd; idcs   = strfind(mydir,'/');
    addpath([mydir(1:idcs(end)-1),'/nshm_inversion']);
    tmp_table=readgeotable('otago_sections_geodetic_rates.geojson');
    unique_flts_id=unique(tmp_table.ParentID); orb_sliprates=zeros(length(unique_flts_id),1) ;
    
    %select slip rates from unique faults only
    for kk=1:length(unique_flts_id)
        tmp_indx=find(tmp_table.ParentID==unique_flts_id(kk),1,'first');

        tmp_sliprate=tmp_table.SlipRate(tmp_indx);
        orb_indx=find(tmp_table.ParentName(tmp_indx)==orb_faults.Name2);
        orb_sliprates(orb_indx,:)=tmp_sliprate;
    end

    aperiodicity=2;   
   
end

%Alternative option to mode if aperiodicity opt = 2:
%1: Use Eq A8 of Zoller (2008) where parameter deltaM2 of Z08 which = mean char mag (m_m) - max G-R event (m_max)  
%2: Is the same as for G-R mdf
    
aperiodicity_opt=1; if aperiodicity_opt==2; char_alpha=2; end

%% Slip rate range and magnitude for each fault from NZCFM

% segmented model

orb_area=orb_faults.Area_m2;%fault area
%Max mag for crustal fault scaling from Stirling et al (2024)
Mx_segmented = [log10(orb_area./10^6)+4.1,log10(orb_area./10^6)+4.2,log10(orb_area./10^6)+4.3]; 

% combined sources

orb_faults_comb=readtable('OtagoRangeBasinFaults.xlsx','Sheet','combined'); 
c_area=zeros(length(orb_faults_comb.NAME),1); w_slip_rate=zeros(length(orb_faults_comb.NAME),1);
num_comb=length(orb_faults_comb.NAME);

%index combined fault components to the NZCFM

for kk=1:length(orb_faults_comb.NAME)
    
    cflt_indx{kk}(1)=find(orb_faults_comb.OBJECTID_1(kk)==orb_faults.OBJECTID);
    
    %not all faults have >1 component, so only searches if subsequent components exist
    if isnan(orb_faults_comb.OBJECTID_2(kk))==0
        cflt_indx{kk}(2)=find(orb_faults_comb.OBJECTID_2(kk)==orb_faults.OBJECTID);
    end
    
    if isnan(orb_faults_comb.OBJECTID_3(kk))==0
        cflt_indx{kk}(3)=find(orb_faults_comb.OBJECTID_3(kk)==orb_faults.OBJECTID);
    end
    
    if isnan(orb_faults_comb.OBJECTID_4(kk))==0
        cflt_indx{kk}(4)=find(orb_faults_comb.OBJECTID_4(kk)==orb_faults.OBJECTID);
    end
    
    %Extract total area of combined source
    for mm=1:width(cflt_indx{kk})
        tmp_indx=cflt_indx{kk}(mm);
        %extract combined area of source
        c_area(kk)=orb_faults.Area_m2(tmp_indx)+c_area(kk);
        %for each componend, extract slip rates*area

        if recurrence_model==4
           tmp_slip_rate(mm,:)=orb_faults.Area_m2(tmp_indx)*orb_sliprates(tmp_indx);

        else
            tmp_slip_rate(mm,:)=orb_faults.Area_m2(tmp_indx)*orb_faults.SR_pref(tmp_indx);
        end                 
    end
    
    %extract area weighted average slip rate for min, pref, and max
    w_slip_rate(kk)= sum(tmp_slip_rate)/c_area(kk);
    clear tmp_slip_rate
    
end

Mx_combined = [log10(c_area./10^6)+4.1,log10(c_area./10^6)+4.2,log10(c_area./10^6)+4.3];

%% Logic tree options

rigidity = 3*10^10; % N/m^2, set so consistent with 2022 NSHM (Gerstenberger et al 2024)
Mmin = 5; Mmax=7.8;

%assymetric weighting for slip rate and Mmax as higher estimates are disproportionately high
Mmax_shift_weight = [0.2 0.7 0.1];

dm = 0.005;%discretized magnitude interval
mag_range_GR={};

%Parameter from Y&C 85 that governs magnitude increment over which fault
%hosts characteristic behaviour
DeltaM1_segmented = zeros(length(Mx_segmented),1);
DeltaM1_combined = zeros(length(Mx_combined),1);

%Parameter from Y&C 85 that governs magnitude increment between G-R and
%characteristic behaviour
DeltaM2_segmented = zeros(length(Mx_segmented),1); 
DeltaM2_combined = zeros(length(Mx_combined),1);

%Progressively scale Delta M1 and M2 values depending on size of Mmax>Mmin
%NOTE, with Mmin=5, DeltaM1 always equals 1 and DeltaM2 equals 0.5
for i=1:num_fault
    if Mx_segmented(i)>(Mmin+1.5)
        DeltaM1_segmented(i) = 1.0; DeltaM2_segmented(i) = 0.5; 
    elseif Mx_segmented(i)>(Mmin+1)
        DeltaM1_segmented(i) = 0.7; DeltaM2_segmented(i) = 0.3; 
    elseif Mx_segmented(i)>(Mmin+0.5)
        DeltaM1_segmented(i) = 0.3; DeltaM2_segmented(i) = 0.2;     
    else
        DeltaM1_segmented(i) = 0.2; DeltaM2_segmented(i) = 0.1;
    end
end

for i=1:num_comb
    if Mx_combined(i)>6.5
        DeltaM1_combined(i) = 1.0; DeltaM2_combined(i) = 0.5; 
    else
        DeltaM1_combined(i) = 0.5; DeltaM2_combined(i) = 0.2;
    end
end

%NZ crustal MFD (Gerstenberger et al 2024)
b = 0.96;

%random variation of b value with uncertainites from Gerstenbenger et al 2024
Bvalue_option = [b+0.14 b b-0.14];
Bvalue_weight = [0.17 0.66 0.17];


%Convert fault slip rate to m/yr
%For combined case, also correct for rupture weighting
orb_sliprates_m = orb_sliprates.*10^-3;

w_slip_rate_m = zeros(length(Mx_combined),1);
for i=1:num_comb
    w_slip_rate_m(i)=w_slip_rate(i).*10^-3.*orb_faults_comb.Weighting(i);
end


%% Logic tree set up for different Mmax, and b value options

count_case = 0;
RecurrenceVar_segmented = {}; RecurrenceVar_combined = {};

for ff = 1:num_fault
        for jj = 1:length(Bvalue_option)
            for kk = 1:length(Mmax_shift_weight)
                count_case = count_case + 1;                
                    RecurrenceVar_segmented{ff}(count_case,1:4) = [count_case Bvalue_option(jj) Mx_segmented(ff,kk) Bvalue_weight(jj)*Mmax_shift_weight(kk)];
            end
    end
   count_case = 0;
end

for ff = 1:num_comb
        for jj = 1:length(Bvalue_option)
            for kk = 1:length(Mmax_shift_weight)
                count_case = count_case + 1;                
                    RecurrenceVar_combined{ff}(count_case,1:4) = [count_case Bvalue_option(jj) Mx_combined(ff,kk) Bvalue_weight(jj)*Mmax_shift_weight(kk)];
            end
        end
   count_case = 0;
end

%% Run fault recurrence model for segment (char) and combined (char and G-R)

variable_name1=strcat('model',num2str(recurrence_model),'/orb_fault_segmented');
variable_name2=strcat('model',num2str(recurrence_model),'/orb_fault_combined');

%Run first for segmented case (ee=1), and then for combined case (ee=2)
%Note G-R parameters are calculated for segmented case, but these are not then used in stochastic event catalog

for ee=1:2
    
    %Save variables in cell arrays for each fault
    alpha_var1={}; fm_char_var1={}; fm_exp_var1={};
    alpha_bpt1={};  T={}; lambda={}; mag_range_GR={};%variables for BPT
    YC85SeismicMoment_exp={}; YC85SeismicMoment_char={}; %variables to store moment rates in
    
    if ee==1
        num_sources=num_fault; RecurrenceVar=RecurrenceVar_segmented; source_area=orb_area;
        DeltaM1=DeltaM1_segmented; DeltaM2=DeltaM2_segmented; slip_rates=orb_sliprates_m;
  
    else
        num_sources=num_comb;  RecurrenceVar=RecurrenceVar_combined; source_area=c_area;
        DeltaM1=DeltaM1_combined; DeltaM2=DeltaM2_combined;  slip_rates=w_slip_rate_m;
        
    end
    
    for ff= 1:num_sources %loop through each fault/multi_fault
        
        Mx=max(RecurrenceVar{ff}(:,3));
        mag_range_GR{ff} = ((Mmin-dm/2):dm:(Mx+0.15+dm/2))';
        
        SeismicMoment_char = zeros(height(RecurrenceVar{ff}),3);
        SeismicMoment_exp  = zeros(height(RecurrenceVar{ff}),3);

        weighted_GR_char = zeros(length(mag_range_GR{ff}),1);
        weighted_GR_exp  = zeros(length(mag_range_GR{ff}),1);   
  
        for jj = 1:height(RecurrenceVar{ff})
        
            % Input variables of para_GR:
            % 1) Rigidity (Nm)
            % 2) Fault plane area Af (m^2) for fault ff and case jj
            % 3) Slip rate S (m/year) for fault ff and case jj
            % 4) Gutenberg-Richter slope parameter b
            % 5) Maximum magnitude Mmax for fault ff
            % 6) Minimum magnitude Mmin
            % 7) Magnitude range for exponential distribution part DeltaM1
            % 8) Magnitude range for characteristic distribution part DeltaM2 
        
            para_GR = [rigidity source_area(ff) slip_rates(ff) RecurrenceVar{ff}(jj,2) RecurrenceVar{ff}(jj,3) Mmin DeltaM1(ff) DeltaM2(ff)];

            %Run YC1985 function where inputs are: 1) para_GR that includes SR and Mmax info,
            %2) mag range for GR relation, and 3) tick labels for a G-R plot
        
            %Outputs are activity rate and pdf and cdf for GR and characteristic behaviours
            [alpha_var(jj,1:3),fm_char_var(1:length(mag_range_GR{ff}),3*(jj-1)+1:3*jj),fm_exp_var(1:length(mag_range_GR{ff}),3*(jj-1)+1:3*jj)] = ...
                characteristic_magnitude_YC1985(para_GR,mag_range_GR{ff},[],[6 8 10^-5 10^-3]);
        
            %weighted mean of frequency for each magnitude for case where width =1 and width =2
        
            %Note weighted average cuts out at mmax-mshift (i,e lower bound of
            %Mmax). Not true MFD can't be used for moment rate calc
            weighted_GR_char(:,1) = weighted_GR_char(:,1) + RecurrenceVar{ff}(1,4)*fm_char_var(:,3*jj);
            weighted_GR_exp(:,1)  = weighted_GR_exp(:,1)  + RecurrenceVar{ff}(1,4)*fm_exp_var(:,3*jj);
          
            for kk = 1:length(mag_range_GR{ff})
                if fm_char_var(kk,3*jj) > 0
                    %moment rate for each case from assessing rate from each
                    %magnitude interval
                    SeismicMoment_char(jj) = SeismicMoment_char(jj) + (fm_char_var(kk,3*jj)-fm_char_var(kk+1,3*jj))*(10^(1.5*mag_range_GR{ff}(kk)+9.05));
                end
                
                if fm_exp_var(kk,3*jj) > 0
                    SeismicMoment_exp(jj) = SeismicMoment_exp(jj) + (fm_exp_var(kk,3*jj)-fm_exp_var(kk+1,3*jj))*(10^(1.5*mag_range_GR{ff}(kk)+9.05));
                end
            
            end
        
            YC85SeismicMoment_exp{ff}(jj,1)=SeismicMoment_exp(jj,1);
            YC85SeismicMoment_char{ff}(jj,1)=SeismicMoment_char(jj,1);
        
            if aperiodicity_opt==1
                alpha_bpt{ff}(jj,:)=[ones(1,2).*aperiodicity, RecurrenceVar{ff}(jj,4)];
            else
                %Derive aperiodicity from the b-value of a G-R relationship using Eq 21 of Zoller (2008);
                alpha_bpt{ff}(jj,1)=(RecurrenceVar{ff}(jj,2)/(3-RecurrenceVar{ff}(jj,2)))^0.5;
            
                if char_alpha==1
                    %Aperiodicity specific for characteristic mag using Eq A8 of Zoller (2008)
                    alpha_bpt{ff}(jj,2:3)=[10^((-3/2)*(DeltaM1+DeltaM2*0.5))*(RecurrenceVar{ff}(jj,2)/(3-RecurrenceVar{ff}(jj,2)))^0.5 RecurrenceVar{ff}(jj,4)];
                else
                    %Aperiodicity from b-value of a G-R relationship using Eq 21 of Zoller (2008);
                    alpha_bpt{ff}(jj,2:3)=[alpha_bpt{ff}(jj,1) RecurrenceVar{ff}(jj,4)];
                end
        
            end %end aperiodicity if statement
            
            %Mean recurrence time for BPT as inverse of activity rate for fault ff for (1) char and (2) GR
            T{ff}(jj,:)= [(1./(alpha_var(jj,1)+alpha_var(jj,2))) 1./(alpha_var(jj,3))];
        
         end  %end jj loop for different Mmax, b-value and slip rate variations
    
    
        %Derive shape parameter for use with inverse Gaussian functions
        %See eq.4 in Chapter 5 of Griffin (2021)
        lambda{ff}=[T{ff}(:,1)./alpha_bpt{ff}(:,1).^2,  T{ff}(:,2)./alpha_bpt{ff}(:,2).^2];

        % Assign fault ID's to data
        if ee==1
            alpha_var = horzcat(alpha_var,ones(height(alpha_var),1)*orb_faults.OBJECTID(ff));
        else
            alpha_var = horzcat(alpha_var,ones(height(alpha_var),1)*orb_faults_comb.Combined_ID(ff));
        end
   
        alpha_var1{ff} = alpha_var; fm_char_var1{ff} = fm_char_var; fm_exp_var1{ff} = fm_exp_var;
        
        clear alpha_var fm_char_var fm_exp_var 
        clear SeismicMoment_exp  SeismicMoment_char
    
    end % end ff loop for source

    if ee==1
    
    save(variable_name1,'RecurrenceVar','alpha_var1','fm_char_var1','fm_exp_var1','lambda','orb_faults','alpha_bpt','T',...
    'YC85SeismicMoment_exp','YC85SeismicMoment_char','num_fault','mag_range_GR');
    
    elseif ee==2
     
    save(variable_name2,'RecurrenceVar','alpha_var1','fm_char_var1','fm_exp_var1','lambda','orb_faults_comb','alpha_bpt','T',...
    'YC85SeismicMoment_exp','YC85SeismicMoment_char','num_comb','mag_range_GR','c_area','w_slip_rate_m');

    end
    
    
end %end ee loop

%% Save variables

variable_name3=strcat('model',num2str(recurrence_model),'/orb_fault_parameters');

save(variable_name3,'num_simu','t_limit','dm','rigidity','Mmin','Mmax','aperiodicity_opt');
