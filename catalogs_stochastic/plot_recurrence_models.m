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

%Select fault by (1) NZCFM ID or (2) Combined ID

sf_id=40; %NZCFM ID for Pisa
cf_id=[31 32]; %Combined ID for Pisa-Cluden and Pisa-Grandview
select_opt=[1 2]; 

col_opt=vertcat([0.87 0.32 0.08],[0.33 0.85 0.25]); 

%Plot median magnitude-frequency distributions for fault
%conisders both characteristic and G-R if combined source

figure(1);

for ss=1:length(select_opt)

select_opt=ss;

    if select_opt==1

        f_indx=find(orb_faults.OBJECTID==sf_id); syncat_opt=1; indx_opt=0;
        load('orb_fault_segmented','mag_range_GR','fm_char_var1','fm_exp_var1',...
        'RecurrenceVar','T','lambda');

        median_ri=(length(RecurrenceVar{f_indx})/2)+0.5;

        p1=semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,median_ri*3),'Color',[0.58,0.49,0.86],'LineWidth',0.5); hold on;

    elseif select_opt==2

        load('orb_fault_combined','mag_range_GR','fm_char_var1','fm_exp_var1',...
            'RecurrenceVar','T','lambda');

        for ii=1:length(cf_id)

            f_indx=find(orb_faults_comb.Combined_ID==cf_id(ii)); syncat_opt=[1,2]; indx_opt=1;

            p2=semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,median_ri*3),'Color',col_opt(ii,:),'LineStyle','-','LineWidth',0.5); hold on;
            p3=semilogy(mag_range_GR{f_indx},fm_exp_var1{f_indx}(:,median_ri*3),'Color',col_opt(ii,:),'LineStyle','-','LineWidth',0.5); hold on;
        end
    end
end

legend([p2 p3],{'comb-char','G-R'},...
  'FontSize',13,'Location','southwest'); hold on;

 set(gca,'fontsize',11); axis([Mmin Mmax-0.2 5*10^-6 5*10^-3]); hold on;
    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;


%% Plot fault's mfd for all b-value and Mmax cases

cf_id=30;%plot for the Dunstan Fault

char_col_opt=vertcat([1 0 0],[0.6 0.8 0.9]); gr_col_opt=vertcat([0 0 1],[0.2 0.5 0.2]);

figure(2);

for ii=1:length(cf_id)

   f_indx=find(orb_faults_comb.Combined_ID==cf_id(ii)); syncat_opt=[1,2]; indx_opt=1;

    for jj=1:height(RecurrenceVar{f_indx})
        tmp_width=(RecurrenceVar{f_indx}(jj,4))*3;
        p2=semilogy(mag_range_GR{f_indx},fm_char_var1{f_indx}(:,jj*3),'Color',char_col_opt(ii,:),'LineStyle','-','LineWidth',tmp_width); hold on;
        p3=semilogy(mag_range_GR{f_indx},fm_exp_var1{f_indx}(:,jj*3),'Color',gr_col_opt(ii,:),'LineStyle','-','LineWidth',tmp_width); hold on;
    end

end

legend([p2 p3],{'char','G-R'},'FontSize',13,'Location','southwest'); hold on;

 set(gca,'fontsize',13); axis([Mmin Mmax-0.2 5*10^-6 5*10^-2]); hold on;
    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;