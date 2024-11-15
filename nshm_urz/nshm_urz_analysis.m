%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analyse NSHM URZ (Iturrieta et al 2024) given spatial area 
%%%%%%%%%%%%%%%% of Otago within the URZ %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) Load URZ forecast rates from Github, select cell points within Otago,
%and combine cells to get total rate for M>5 events

%2) Obtain entire URZ rate, and use to model rate URZ rate varability using
%a negative binomial process

%3) Plots MFD of URZ, and random samples with events in NZ Integrated Catalog

%4) Save results

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1)); addpath([mydir(1:idcs(end)-1),'/orb_fault_geometries']);

%load Otago events from NZ integrated catalog 1971-2021 (Rollins et al 2022)
addpath([mydir(1:idcs(end)-1),'/catalogs_instrumental']); load('orb_catalog','otago_aug_catalog_mat3','ins_Mmin','ins_Mmax');

catalog_duration=50; ins_mag_range = [ins_Mmin:0.05:ins_Mmax]; insAnnualRate=zeros(length(ins_mag_range),1);

for mm=1:length(ins_mag_range)
    insAnnualRate(mm) = length(find(otago_aug_catalog_mat3(:,8) >= ins_mag_range(mm)))/(catalog_duration); %annual event rate 
end

%Annual rate of M>5 earthquakes for hybrid-URZ volume and b-value of 0.959
%in NZ to derive absolute rates (Rollins et al 2024)
Nval=4.83;

%load Otago Spatial Polygon
otago_polygon=shaperead([mydir(1:idcs(end)-1),'/gis_files/orb_area_polygon.shp']);

%Load forecast as per Poisson and Non-Poisson floor Multiplicative model
%forecast taken from the floor ensemble method @ https://github.com/pabloitu/nz_nshm2022_nonpoisson/tree/main/forecasts/negbinom_floors/src
addpath('input_files');

wgs84 = wgs84Ellipsoid("km");
forecast=["npfe_m","fe_m","m"];
num_forecast=length(forecast);

otago_forecast=cell(num_forecast,1); otago_dispersion=cell(num_forecast,1);
otago_rate=zeros(num_forecast,2);

%extract Otago data for each DSM forecast:
% (1) negative binomial-multiplacative floor (npfe_m)
% (2) poisson-multiplacative floor (fe_m)
% (3) multiplacative (m)

for ff=1:num_forecast

    forecast_m=readtable(strcat(forecast(ff),".csv")); %grid cells based on cell vertices
    forecast_central_point=readtable(strcat(forecast(ff),".txt")); %grid cells based on central point point of cell (created in 'find_centre_point.m')

    %rearrange npfe_m so cells are sorted by lat in ascending order like txt file
    forecast_m=sortrows(forecast_m,'lat_min','ascend');

    %get area of each latitude-longitude cell
    grid_area = areaquad(forecast_m.lat_min,forecast_m.lon_min,forecast_m.lat_max,forecast_m.lon_max,wgs84);

    %find grid cells in Otago polygon
    filter_s=zeros(height(forecast_central_point),1);

    for ii=1:length(filter_s)
        filter_s(ii)=inpolygon(forecast_central_point.Var1(ii),forecast_central_point.Var2(ii),otago_polygon.X(1:(end-1)),otago_polygon.Y(1:(end-1)));
    end
    
    %store forecast for each dsm for Otago in a cell array
    otago_forecast{ff}=forecast_central_point(find(filter_s==1),:);
    otago_dispersion{ff}=forecast_m.dispersion(find(filter_s==1),:);

    %Each cell contains rate*cell area, so sum to get total area rate
    otago_rate(ff,1)=sum(otago_forecast{ff}.Var3);

   
    %Rate is currently a normalised rate, so scale up rate by N value from
    %Rollins et al (2024)
    otago_rate(ff,2)=otago_rate(ff,1)*Nval;
    
end

%index check, if grid indexing is correct, grid area ~= polygon area
grid_area_otago=sum(grid_area(find(filter_s==1),:));
[otago_utm_x,otago_utm_y]=deg2utm(otago_polygon.Y(1:(end-1)),otago_polygon.X(1:(end-1)));
polygon_area_otago=polyarea(otago_utm_x,otago_utm_y);

%% Save and option to plot map of Otago points
finger_opt=0;

nz_coastline=shaperead([mydir(1:idcs(end)-1),'/gis_files/nz-coastline-polygon1.shp']);

if finger_opt==0

    figure(1);
    tiledlayout(num_forecast,1)

    for ff=1:num_forecast

        nexttile
        plot(forecast_central_point.Var1,forecast_central_point.Var2,'k*');hold on
        plot(otago_forecast{ff}.Var1,otago_forecast{ff}.Var2,'r*'); hold on
        plot(otago_polygon.X(1:(end-1)),otago_polygon.Y(1:(end-1)),'r-','LineWidth',1.1);hold on
        plot(nz_coastline(1).X(1:(end-1)),nz_coastline(1).Y(1:(end-1)),'k-','LineWidth',1);hold on
   
        ylim([-46.8,-44.2]); xlim([168.52 171.12]); axis square

    end
end

for ff=1:num_forecast
    writetable(otago_forecast{ff},strcat("input_files/",forecast(ff),"_otago.txt"),'Delimiter',' ','WriteVariableNames',0);
end

%% Find all cells in low strain rate J2 bins to get entire URZ rate

%select binning method: (3) 3 bin model, (4) 4 bin model
urz_polygon=3;

if urz_polygon==3
    urz_polygon=shaperead('hw_final_j2_2.shp');
    polygon_select=[1 2 3 4 5 6 7 9 10 13]; %polygons in shapefile NOT in lowest URZ (0) in 3 bin model:
elseif urz_polygon==4
    urz_polygon=shaperead('hw_final_j2_3.shp');
    polygon_select=[2 4 8 11];%polygons in shapefile NOT in lowest URZ (0) in 4 bin model
end

%Obtain dispersion for negative binomial URZ (i.e., forecast 1)
%Find points notin URZ low strain rate bin. Then use these points to identify points that ARE within the URZ low strain rate bin

tmp_table1=readtable(strcat(forecast(ff),".txt")); 
filter_z=zeros(height(tmp_table1),1);

for kk=1:length(polygon_select)
    
    for ii=1:length(filter_z)
    
        if filter_z(ii)==0   
            filter_z(ii)=inpolygon(tmp_table1.Var1(ii),tmp_table1.Var2(ii),...
                urz_polygon(polygon_select(kk)).X(1:(end-1)),urz_polygon(polygon_select(kk)).Y(1:(end-1)));
        end
    end
end

%Derive negative binomial forecast rate for low strain rate URZ
urz_points=tmp_table1(find(filter_z==0),:);
tmp_table2=readtable(strcat(forecast(1),".csv"));
tmp_table2=sortrows(tmp_table2,'lat_min','ascend');
urz_dispersion=tmp_table2.dispersion(find(filter_z==0),:); 
urz_rate=sum(urz_points.Var3)*Nval;

%gridded area of urz
grid_area_urz=sum(grid_area(find(filter_z==0),:));

figure_opt=0; %option to plot figure to check points sampled for low strain rate URZ

if figure_opt==1
    
    figure(2);

    plot(urz_points.Var1,urz_points.Var2,'rx');hold on
    
    for kk=1:length(urz_polygon)
        plot(urz_polygon(kk).X(1:(end-1)),urz_polygon(kk).Y(1:(end-1)),'b-');hold on
    end
    
    plot(nz_coastline(1).X(1:(end-1)),nz_coastline(1).Y(1:(end-1)),'k-','LineWidth',1);hold on
    plot(nz_coastline(4).X(1:(end-1)),nz_coastline(4).Y(1:(end-1)),'k-','LineWidth',1);hold on
    plot(nz_coastline(5).X(1:(end-1)),nz_coastline(5).Y(1:(end-1)),'k-','LineWidth',1);hold on
    
    axis equal
end


%% Randomly simulate 1000 other possible rate values for each forecast following their probbility distribution

numSamples = 1000;
rate_samples =zeros(numSamples,num_forecast);

for ff=1:num_forecast

    if ff==1
        %Using negative binomial distribution in terms of p (scale parameter or probability of success)
        %and r (shape parameter or number of successes or shape parameter)). Forecast is for all of URZ, 
        %so needs to be later scaled by area
        r=1/urz_dispersion(1);
        p=1/(1+urz_rate*urz_dispersion(1));

        rate_samples(:,ff) = nbinrnd(r, p, numSamples,1);
    
    elseif ff==2 %forecast modelled using a Poisson distribution
        
        rate_samples(:,ff)  = poissrnd( otago_rate(ff,2)*100,[numSamples 1])./100;

    elseif ff==3 %no uncertainty in this forecast

        rate_samples(:,ff)=otago_rate(ff,2);
    end
end

%% Build MFD for URZ and DSM forecasts and import historical seismicity info for Otago

urz_mmin=4.95; %magitude that URZ is forecasting for

urz_mag_range=[urz_mmin:0.05:8]; urz_mag_rate=zeros(length(urz_mag_range),num_forecast);
urz_mo_rate=zeros(num_forecast,1); urz_mag_rate_sample=zeros(numSamples,length(urz_mag_range),num_forecast);

b=0.959; %(after Rollins et al 2024)

for ff=1:num_forecast

    %MFD for URZ
    a=log10(otago_rate(ff,2))+b*urz_mmin;%derive the a-value to build G-R MFD for each forecast

    for kk=1:length(urz_mag_range)
        urz_mag_rate(kk,ff)=10^(a-b*urz_mag_range(kk));%cummulative rate
    end

    mo_rate_mag_bin=zeros(length(urz_mag_range),1);

    %Moment rate for URZ
    for kk=1:length(urz_mag_range)-1
        inc_rate=urz_mag_rate(kk,ff)-urz_mag_rate(kk+1,ff);%incremental rate for each mag bin
        mo_rate_mag_bin(kk)=10^(1.5*mean([urz_mag_range(kk),urz_mag_range(kk+1)])+9.05)*inc_rate;
    end

    urz_mo_rate(ff)=sum(mo_rate_mag_bin);

    for ll=1:numSamples
        
        if ff==1
            %sampled rates also need to be scaled to proportion of otago in URZ zone
            a=log10(rate_samples(ll,ff)*otago_rate(ff,2)/urz_rate)+b*urz_mmin;
        else
            a=log10(rate_samples(ll,ff))+b*urz_mmin;
        end

        for kk=1:length(urz_mag_range)
            urz_mag_rate_sample(ll,kk,ff)=10^(a-b*urz_mag_range(kk));
        end
    end
end

%derive MFD for Otago events in NZ Integrated catalog 1971-2021
catalog_duration=50;
ins_mag_range = [ins_Mmin:0.05:ins_Mmax]; AnnualRate=zeros(length(ins_mag_range),1);

for mm=1:length(ins_mag_range)
    AnnualRate(mm) = length(find(otago_aug_catalog_mat3(:,8) >= ins_mag_range(mm)))/(catalog_duration); %annual event rate 
end

%% Plot NSHM URZ Forecasts in MFD space and distribution of M>4.95 events

figure(3);

tiledlayout(1,2,"TileSpacing","compact")

fntsize=13;
nexttile

col_opt=vertcat([1,0.5,0],[0 0 1],[0.2 0.8 0.2]);

%plot samples first

for ff=1:num_forecast
    p2=semilogy(urz_mag_range,squeeze(urz_mag_rate_sample(:,:,ff)),'Color',[col_opt(ff,:),0.4],'LineWidth',0.7);hold on
end

%plot overall MFD
for ff=1:num_forecast
    p1(ff)=semilogy(urz_mag_range,urz_mag_rate(:,ff),'-v','Color',col_opt(ff,:),'LineWidth',1.5,'MarkerIndices',[1+ff:5:length(urz_mag_range)]);hold on
end

%plot instrumental record
p3=semilogy(ins_mag_range,insAnnualRate,'-*','Color',[0.4 0.4 0.4],'LineWidth',1.5,'MarkerIndices',[1:5:length(insAnnualRate)]);

set(gca,'FontSize',fntsize);
ylabel('Annual frequency of exceedance'); xlabel('Magnitude');xlim([ins_Mmin 8]); ylim([5*10^-5 2]);

legend([p1(1), p1(2), p1(3), p3],{['neg-binom (empty count: ',num2str(length(find(rate_samples(:,1)==0))),')'],...
['Poisson (empty count: ',num2str(length(find(rate_samples(:,2)==0))),')'],'Hybrid','NZ Integrated Catalog'},...
'Location','southwest','FontSize',fntsize-2)

axis square; grid on

%Plot distribution of M>4.95 events
nexttile

%first fit Kernel distribution to them
nb_dist=fitdist(urz_mag_rate_sample(:,1,1),'Kernel'); hold on
poi_dist=fitdist(urz_mag_rate_sample(:,1,2),'Kernel');

xx=[0:0.01:0.5];
pdf_nb=pdf(nb_dist,xx); pdf_poi=pdf(poi_dist,xx);

%Plot distributions
plot(xx,pdf_nb,'LineWidth',1.5,'Color',col_opt(1,:)); hold on
plot(xx,pdf_poi,'LineWidth',1.5,'Color',col_opt(2,:)); hold on
ymax=gca().YLim(2);
plot([urz_mag_rate_sample(1,1,3),urz_mag_rate_sample(1,1,3)],[0 ymax],'LineWidth',1.5,'Color',col_opt(3,:));

set(gca,'FontSize',fntsize); xlim([0 0.5]); ylabel('Probability'); xlabel('M>=4.95 annual rate')

legend('neg-binom','Poisson','Hybrid','FontSize',fntsize-2); axis square; grid on; box on

set(gcf,'Position',[440 309 796 390])

%% Save output

save('nshm_urz_analysis','num_forecast','urz_mag_range','urz_mag_rate','urz_mag_rate_sample','urz_mo_rate');
