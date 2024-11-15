%% Script to load URZ files from Githib and create new files based on central coordinates of each cell

addpath('input_files');

%select forecast
%m = multiplacative hybid forecast
%fe = Poisson floor ensemble
%npfe = Negative binomial floor ensemble

forecast_select='m';
%load forecast
forecast=readtable(strcat(forecast_select,'.csv'));

%derive central coordinates 
centre_lat=(forecast.lat_min+forecast.lat_max)/2;
centre_lon=(forecast.lon_min+forecast.lon_max)/2;

rate_data=[centre_lon,centre_lat,forecast.rate]; rate_data=sortrows(rate_data,2,'ascend');
lograte_data=[centre_lon,centre_lat,log10(forecast.rate)]; lograte_data=sortrows(lograte_data,2,'ascend');

%Write tables to txt files that can be loaded as rasters onto QGIS
writetable(table(rate_data),strcat('input_files/',forecast_select,'.txt'),'Delimiter',' ','WriteVariableNames',0);
writetable(table(lograte_data),strcat('input_files/',forecast_select,'_log10.txt'),'Delimiter',' ','WriteVariableNames',0);

