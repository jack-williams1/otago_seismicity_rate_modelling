%% Otago earthquake chronology simulations

nsim_paleoseis=10000;

addpath('assymmetric_dist')

%the order of these earthquake ages MUST match their order in OtagoEarthquakeTimings.xlsx
otago_asymmetric_eq_indx=["NWCardrona_eq1","Dunstan_eqF_5event","Dunstan_eqE_6event","Dunstan_eqF_6event","NWCardrona_eq2"];

otago_eq_timings=readtable('OtagoEarthquakeTimings.xlsx');

simu_eq_timings=zeros(nsim_paleoseis,height(otago_eq_timings)); count=0;


for ii=1:height(otago_eq_timings)
    
    if otago_eq_timings.uncertainty(ii)==0
        
        count=count+1;
        
        %load earthquake ages with assymmetric distribution
        eq_data=table2array(readtable(strcat(otago_asymmetric_eq_indx(count),'.txt')));
        
        %derive cumulative probability that a certain earthquake age will be sampled
        eq_cum_prob=[zeros(length(eq_data(:,2)),1),eq_data(:,1)];

        for jj=1:length(eq_data(:,2))-1
            eq_cum_prob(jj+1,1)=eq_data(jj,2)+eq_cum_prob(jj,1);
        end
        
        ran_samp=rand(nsim_paleoseis,1).*sum(eq_data(:,2));
        
        %randonly sample eq age distribution by finding the min number 
        %between random sample and an ages cumulative probability
        for kk=1:nsim_paleoseis
            [~,indx]=min(abs(ran_samp(kk)-eq_cum_prob(:,1)));
            if count==1 || count==5 %correction for NW Cardrona ages
                simu_eq_timings(kk,ii)=2019+(-1*eq_cum_prob(indx,2));
            else
                simu_eq_timings(kk,ii)=2006+(-1*eq_cum_prob(indx,2));
            end
        end
        
    elseif ii==16 %tread Titri timing as a box car function
        simu_eq_timings(:,ii)=randi([otago_eq_timings.Lower_estimate_ka(ii),otago_eq_timings.Upper_estimate_ka(ii)],nsim_paleoseis,1);


    else %ages sampled from a normal distribution
        simu_eq_timings(:,ii)=normrnd(otago_eq_timings.Mean_ka(ii),otago_eq_timings.uncertainty(ii)*0.5,[nsim_paleoseis,1]);
    end
    
end

%% Simulate 10,000 random time windows between 0-20 ka

time_window=zeros(nsim_paleoseis,2); time_window_length=10000;
time_window(:,1)=randi([0,time_window_length],nsim_paleoseis,1);
time_window(:,2)=time_window(:,1)+time_window_length;

ran_LT1 = rand(nsim_paleoseis,1); ran_LT2 = rand(nsim_paleoseis,1);

num_earthquake=zeros(nsim_paleoseis,1);

for kk=1:nsim_paleoseis
    
    tmp_indx=find(time_window(kk,1)<simu_eq_timings(kk,:)' & time_window(kk,2)>simu_eq_timings(kk,:)' )';
    
    num_earthquake(kk)=width(tmp_indx);
    
    %If time window sampled Dunstan Fault event D (tmp_indx=10), 50% of these
    %counts are removed as event is only 'possible'
    if ran_LT1(kk)<0.5 && isempty(find(tmp_indx==10))==0
        num_earthquake(kk)=num_earthquake(kk)-1;
    end

   % 50% of time windows include Dunstan Fault Event E and F in 6 event scenario
   % (so remove Event F in 5 event scenario)
    if ran_LT2(kk)<0.5 && isempty(find(tmp_indx==11))==0
        num_earthquake(kk)=num_earthquake(kk)-1;
   % 50% of time windows include Dunstan Fault F in 5 event scenario
   % (so remove two events if both Events E and F in 6 event scenario sampled)
    elseif ran_LT2(kk)>=0.5 && isempty(find(tmp_indx==12))==0 && isempty(find(tmp_indx==13))==0 
        num_earthquake(kk)=num_earthquake(kk)-2;
     %(or remove 1 events if only Events E or F in 6 event scenario sampled)
    elseif ran_LT2(kk)>=0.5 && (isempty(find(tmp_indx==12))==0 || isempty(find(tmp_indx==13))==0)
        num_earthquake(kk)=num_earthquake(kk)-1;
    end

end

%% Save output

paleoseismic_rate_constraint=[mean(num_earthquake),std(num_earthquake)];
paleoseismic_mag_constraint=[7.1 0.2];
save('otago_paleoseismic_constraint','paleoseismic_rate_constraint','nsim_paleoseis','paleoseismic_mag_constraint','time_window_length')

%% Visualise simulated eq timings as check (optional)

tmp_indx=16;%select earthquake to check

pd=fitdist(simu_eq_timings(:,tmp_indx),'Kernel');

x=0:1:50000;
y=pdf(pd,x);
figure(1);
plot(x,y,'LineWidth',2);

