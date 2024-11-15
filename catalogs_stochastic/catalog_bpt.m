%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   SYNTHETIC CATALOG FOR BROWNIAN PASSAGE TIME INTEREVENT TIMES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%Select recurrence model for Otago Faults to run:
%1) NZCFM slip rates - Aperiodicity = 0.5
%2) NZCFM slip rates - Aperiodicity = 2
%3) NZCFM slip rates - Aperiodicity = 4
%4) Geodetic slip rates - Aperiodicity = 2

model_opt=2;

model_path=strcat('model',num2str(model_opt));
addpath(model_path)

%Save load state variable info too? (note this is a v.large parameter!)
save_opt=2; %don't save if opt=1

load orb_fault_parameters

disp(['Total duration of simulation: ',num2str(num_simu*t_limit),' years']);

%creation of 3 catalogs:
%1) segmented-char
%2) combined-char
%3) combined-GR

for ww=1:3 
     StochasticEventCatalog{ww}=zeros(4*num_simu,8);
end

%% Run simulations

clock_old = cputime; count_case=0; catalog_count=0;

cat_load_state=cell(3,1); cat_cum_disp_bpt=cell(3,1);

for ee=1:2

    if ee==1
        load('orb_fault_segmented','T','alpha_bpt','fm_char_var1','fm_exp_var1','RecurrenceVar','orb_faults','mag_range_GR','num_fault');
        syncat_opt=1; source_id=orb_faults.OBJECTID; num_sources=num_fault;
    elseif ee==2
        load('orb_fault_combined','T','alpha_bpt','fm_char_var1','fm_exp_var1','RecurrenceVar','orb_faults_comb','mag_range_GR','num_comb');
        %syncat opt: both char (1) and GR (2) explored
        syncat_opt=[1,2]; source_id=orb_faults_comb.Combined_ID; num_sources=num_comb;
    end
        
    % Weightings for slip rate-b-value-Mmax logic tree
    recurrence_para_weight = RecurrenceVar{1}(:,4)';
    % Increments between fault values in alpha parameter
    alpha_inc=length(recurrence_para_weight);

    % Combine cell arrays from each fault for use in syncat
    % Last column gives fault id
    alpha_comb = cell2mat(alpha_bpt');
    t_comb = cell2mat(T');
    %However, in this case just using it to make synthetic event catalog

    for ss=1:length(syncat_opt)

        global_count = 0; num_event_tmp=zeros(num_simu,1); catalog_count=catalog_count+1;
        tmp_count = 1; cum_disp_bpt=zeros(num_simu,num_sources); count_bpt=1;
    
        recurrence_type=syncat_opt(ss); %select recurrence_type
    
        %select on-fault MFD for recurrence model
        if ss==1
            fm_var=fm_char_var1;
        else
            fm_var=fm_exp_var1;
        end
    
        load_state=zeros(num_simu+1,num_sources); load_state(1,:)=rand(1,num_sources); %start load state variable somewhere random between 0-1

        for ii = 2:num_simu+1%Run event catalog for duration set by N 
    
            if ii == tmp_count*50000
                tmp_count = tmp_count + 1;
                clock_new = cputime;
                disp(['Current simulation cycle : ',num2str(ii*t_limit),' out of ',num2str(num_simu*t_limit),' & Required time is: ',...
                num2str(clock_new-clock_old), '. Catalog is ', num2str(catalog_count), ' out of 3']);
                clock_old = clock_new;
            end
  
            % Logic tree sampling
            ran_LT = rand(1,1);%1 random numbers between 0 and 1  
   
            % Recurrence parameters (slip rate, Bvalue, and Mmax shift)
            for jj = 1:length(recurrence_para_weight)
                if ran_LT <= sum(recurrence_para_weight(1:jj))
                    recurrence_para = jj; 
                break
                end
            end     
       
            %Sample mean recurrence time for each fault based on reccurrence para
            t_sample = t_comb((recurrence_para:9:length(t_comb)),ss);
            %sigma model calculated from its relationship with aperiodicity and mean activity rate
            %see appendix in powerpoint slides
            sigma =  (alpha_comb((recurrence_para:9:length(t_comb)),ss).^2.*((1./t_sample))).^0.5;
            
            %BPT model for earthquake occurrence 
            if max(load_state(ii-1,:))>=t_limit %record event
       
               f_indx = find(load_state(ii-1,:)>=t_limit); %index faults for which event occurred         

                for ff=1:length(f_indx) %for faults that ruptured
         
                    %sample magnitude from cdf created using Y&C 85 method
                    recurrence_info1 = fm_var{f_indx(ff)}(:,3*(recurrence_para-1)+2);
                    mw = mag_range_GR{f_indx(ff)}(find(rand(1,1)>recurrence_info1,1,'last')+1) + rand(1,1)*dm - dm/2;
                    r_area = 10^6*10^(mw-4.125); %magnitude-area scaling from Stirling et al 2023
                    sed = 10^(1.5*mw+9.05)/(r_area*rigidity);
         
                    % StochasticEventCatalog
                    % 1    ) Event number
                    % 2    ) Simulation cycle
                    % 3    ) Load state
                    % 4    ) Earthquake magnitude
                    % 5    ) NCZFM ID
                    % 6    ) Recurrence type (1 = characteristic versus 2 = exponential)
                    % 7    ) Reccurrece parameter
                    % 8    ) SED
             
                    StochasticEventCatalog{catalog_count}(count_bpt,1:8) = [count_bpt ii-1*t_limit load_state(ii-1,f_indx(ff)) mw source_id(f_indx(ff),1) recurrence_type recurrence_para sed];
  
                    cum_disp_bpt(ii,f_indx(ff))=cum_disp_bpt(ii-1,f_indx(ff))+sed; %Fault cumulative displacement
                    count_bpt=count_bpt+1;
                    load_state(ii,f_indx(ff))=0; %reset fault loading where e/q occurred
                end %end ff loop
        
                %index faults for which no event occurred and update load state
                if isreal(find(load_state(ii-1,:)<1))==1
                    n_f_indx = find(load_state(ii-1,:)<1);
                    load_state(ii,n_f_indx)=load_state(ii-1,n_f_indx)+(1./t_sample(n_f_indx))'+(sigma(n_f_indx).*normrnd(0,1,[length(n_f_indx),1]))';
                    cum_disp_bpt(ii,n_f_indx)=cum_disp_bpt(ii-1,n_f_indx);
                end
        
                num_event_tmp(ii,1) = length(f_indx); %Num_event for each simulation in catalog     
            
            else  %no event occurred 
                  %update load state of each fault following Eq.7 from Matthews et al (2002)
                  load_state(ii,:)=load_state(ii-1,:)+(1./t_sample)'+(sigma.*normrnd(0,1,[num_sources,1]))';
                  cum_disp_bpt(ii)=cum_disp_bpt(ii-1);
                  num_event_tmp(ii,1) = 0;
            
            end %end if statement for load state   
            
        end %end simulation cycle
    
        %Remove unnecessary zeros;
        StochasticEventCatalog{catalog_count}=StochasticEventCatalog{catalog_count}(1:sum(num_event_tmp),:);
        num_event{catalog_count}=num_event_tmp;
    
        cat_load_state{catalog_count}=load_state;%Store load state for each fault
        cat_cum_disp_bpt{catalog_count}=cum_disp_bpt;
    end %end ss loop
    
    %clear variables for segmented case
    clear alpha_var1 fm_char_var1 fm_exp_var1 RecurrenceVar mag_range_GR
    
end  %end ee loop

%% Save results

StochasticEventCatalog_bpt=StochasticEventCatalog;

variable_name=strcat('model',num2str(model_opt),'/catalog_bpt');

if save_opt==1
    save(variable_name,'StochasticEventCatalog_bpt','num_simu');
else %save load state info
    save(variable_name,'StochasticEventCatalog_bpt','cat_load_state','cum_disp_bpt','num_simu');
end