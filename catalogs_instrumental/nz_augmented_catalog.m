%% Select Otago range and basin events from Augmented NZ Catalog and plot MFD %%

close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%Read in NZ augmented earthquake catalog (Rollins et al., 2022: https://geodata.nz/geonetwork/srv/api/records/00876699-4791-4528-a826-1d7ac3c77516) 
aug_catalog=readtable('NZNSHM2022_augmentedEQcatalogue_annotated_2Aug2022');

%Extract following parameters from augmented catalog:
%1) Public ID
%2) Year
%3) Month
%4) Day
%5) Longitude
%6) Latitude
%7) Depth (km)
%8) Magnitude
aug_catalog_mat=table2array(aug_catalog(:,[1 2 3 4 8 9 11 14]));

%filter out events with depths>dZ
dZ=30; %max depth to cont crust
tmp_indx=find(aug_catalog_mat(:,7)<dZ);
aug_catalog_mat=aug_catalog_mat(tmp_indx,:);

%filter for only events in Otago range and basin
addpath([mydir(1:idcs(end)-1) '/gis_files']);
spatialfilter=shaperead('orb_area_polygon.shp');

filter_s=zeros(length(aug_catalog_mat),1);

for ii=1:length(aug_catalog_mat)
    filter_s(ii)=inpolygon(aug_catalog_mat(ii,5),aug_catalog_mat(ii,6),spatialfilter(1).X(1:end-1),spatialfilter(1).Y(1:end-1));
   
end

%Catalog-1: for all Otago events
otago_aug_catalog_mat1=aug_catalog_mat(find(filter_s(:,1)==1),:);

%Catalog-2: for Otago events 1951-2021 (minimum amount of time in Otago with no M>5 events)

yearStart=1951; yearEnd=2021; catalog_duration=yearEnd-yearStart; ins_Mmin=2.0;ins_Mmax=5;

otago_aug_catalog_mat2 = otago_aug_catalog_mat1(find(otago_aug_catalog_mat1(:,2) >= yearStart &...
   otago_aug_catalog_mat1(:,2) <= yearEnd ),:);

%Catalog-3: for Otago events 1971-2021 (50-year period, allows comparison to DSM forecasts)

yearStart=1971; yearEnd=2021; catalog_duration=yearEnd-yearStart; ins_Mmin=2.0;ins_Mmax=5;

otago_aug_catalog_mat3 = otago_aug_catalog_mat1(find(otago_aug_catalog_mat1(:,2) >= yearStart &...
   otago_aug_catalog_mat1(:,2) <= yearEnd ),:);


mo_rate=sum(10.^(otago_aug_catalog_mat2(:,8).*1.5+9.05))/catalog_duration;

%% Plot catalog

mag_range = [ins_Mmin:0.05:ins_Mmax]; AnnualRate=zeros(length(mag_range),1); CumRate=zeros(length(mag_range),1);

for mm=1:length(mag_range)
    CumRate(mm) = length(find(otago_aug_catalog_mat1(:,8) >= mag_range(mm))); %cumulative events
    AnnualRate(mm) = length(find(otago_aug_catalog_mat2(:,8) >= mag_range(mm)))/(catalog_duration); %annual event rate 
end

%Fig option =1 plot mfd for cumulative event rate
%Fig option =2 plot mfd for annual rate 1951-2021
fig_option=2;

if fig_option ==1
     events=CumRate; ylabel_opt='Cumulative number of events'; ylim_opt=[0 100];
elseif fig_option ==2
     events=AnnualRate; ylabel_opt='Annual frequency of exceedance'; ylim_opt=[10^-2 10^1];
end

figure(23);

semilogy(mag_range,events,'b-','LineWidth',1.5); hold on; %MFD 

set(gca,'fontsize',13);  axis([ins_Mmin 6 ylim_opt(1) ylim_opt(2)]); hold on;
xlabel('Magnitude'); ylabel(ylabel_opt); grid on; axis square;

%% Model catalog G-R relationship using Weichart method

%NOTE THIS REQUIRES AN "OPIMISTIC" ESTIMATE FOR MMIN (M_MIN = 2.5) IN THIS CATALOG!
%INCLUDED HERE FOR DEMONSTRATION PURPOSES ONLY. DO NOT USE IN FORMAL ANALYSIS!
completeYear = [2.5 3.0 yearStart; 3.0 3.5 yearStart; 3.5 4.0 yearStart];

for kk = 1:height(completeYear)
    tmp2 = find(otago_aug_catalog_mat2(:,8) >= completeYear(kk,1) & otago_aug_catalog_mat2(:,8) < completeYear(kk,2) & completeYear(kk,3) <= otago_aug_catalog_mat2(:,2)); 
    countW(kk) = length(tmp2);
    timeW(kk)  = yearEnd - completeYear(kk,3) + 1;
    %clear tmp2

end

GRFitW = GRrelation_MLEWeichert_EQMATca(countW,2.75:0.5:3.75,0.5,timeW);

num_sample = 10000;
N0 = exp((log(GRFitW(3)) - 0.5*log(1+(GRFitW(4)/GRFitW(3))^2)) + sqrt(log(1+(GRFitW(4)/GRFitW(3))^2))*normrnd(0,1,num_sample,1));
b  = GRFitW(1) + GRFitW(2)*normrnd(0,1,num_sample,1);

log10N = log10(N0)*ones(1,length(ins_Mmin:0.1:ins_Mmax)) - b*(ins_Mmin:0.1:ins_Mmax);

for kk = 1:width(log10N)
      Stat_log10N(kk,1:7) = prctile(log10N(:,kk),[5 10 16 50 84 90 95]);
end

% Plot mfd of Otago Crustal Events with uncertainity modelling

if fig_option==2 %only plot modelled G-R if plotting annual rates
    semilogy(ins_Mmin:0.1:ins_Mmax,10.^(Stat_log10N(:,4)),'r-'); hold on %Plot extrapolated G-R curve
    semilogy(ins_Mmin:0.1:ins_Mmax,10.^(Stat_log10N(:,[ins_Mmin ins_Mmax])),'r--'); hold on %Plot extrapolated G-R curve uncertainity


    semilogy(mag_range,events,'b-','LineWidth',1.5); hold on; %MFD 

    set(gca,'fontsize',13);  axis([ins_Mmin ins_Mmax ylim_opt(1) ylim_opt(2)]); hold on;
    xlabel('Magnitude'); ylabel(ylabel_opt); grid on; axis square;
end

%% Save Parameters

save('orb_catalog','otago_aug_catalog_mat1','otago_aug_catalog_mat2','otago_aug_catalog_mat3','catalog_duration','ins_Mmin','ins_Mmax','Stat_log10N');

%Write catalog to csv file

tmp=array2table(otago_aug_catalog_mat2,'VariableNames',{'Public ID','Year','Month',...
    'Day','Longitude','Latitude','Depth (km)','Magnitude'});
writetable(tmp,'orb_catalog.csv');

