%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   SYNTHETIC CATALOG FOR POISSON INTEREVENT TIMES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%Select recurrence model for Otago Faults to run:
%1) NZCFM slip rates - Aperiodicity = 0.7
%2) NZCFM slip rates - Aperiodicity = 2
%3) NZCFM slip rates - Aperiodicity = 4
%4) Geodetic slip rates - Aperiodicity = 2

model_opt=2;

model_path=strcat('model',num2str(model_opt));
addpath(model_path)

load orb_fault_parameters

disp(['Total duration of simulation: ',num2str(num_simu*t_limit),' years']);

%creation of 3 catalogs:
%1) segmented-char
%2) combined-char
%3) combined-GR

for ww=1:3 
     StochasticEventCatalog{ww}=zeros(4*num_simu,7);
end

%% Run simulations

clock_old = cputime; count_case=0; catalog_count=0;

for ee=1:2
    
    if ee==1
        load('orb_fault_segmented','alpha_var1','fm_char_var1','fm_exp_var1','RecurrenceVar','orb_faults','mag_range_GR');
        syncat_opt=1; source_id=orb_faults.OBJECTID;
    elseif ee==2
        load('orb_fault_combined','alpha_var1','fm_char_var1','fm_exp_var1','RecurrenceVar','orb_faults_comb','mag_range_GR');
        %syncat opt: both char (1) and GR (2) explored
         syncat_opt=[1,2];  source_id=orb_faults_comb.Combined_ID;
    end
    
    % Weightings for slip rate-b-value-Mmax logic tree
    recurrence_para_weight = RecurrenceVar{1}(:,4)';
    % Increments between fault values in alpha parameter
    alpha_inc=length(recurrence_para_weight);
    
    % Combine cell arrays from each fault for use in syncat
    % Last column gives fault id
    alpha_var_comb = cell2mat(alpha_var1');
    
    %However, in this case just using it to make synthetic event catalog
    %make catalog
    
    for ss=1:length(syncat_opt)

        global_count = 0; num_event_tmp=zeros(num_simu,1);
        tmp_count = 1;  catalog_count=catalog_count+1;
    
        recurrence_type=syncat_opt(ss); %select recurrence_type: %1=char, %2=segmented
    
        %select on-fault MFD for recurrence model
        if recurrence_type==1
            fm_var=fm_char_var1;
        else
            fm_var=fm_exp_var1;
        end
    
        for ii = 1:num_simu %Run event catalog for duration set by N 
    
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
       
            % Characteristic events
            if recurrence_type == 1
                %activity rate for characteristic events, includes both 
                %exponential and non-exponential part of MFD
                %index from alpha_var variable at intervals of number of faults
                alpha_info_tmp = alpha_var_comb((recurrence_para:alpha_inc:length(alpha_var_comb)),1:2);
                alpha_info = sum(alpha_info_tmp,2);
            
            clear alpha_info_tmp
            
            %G-R events    
            elseif recurrence_type == 2
                %Chose activity rate and reccurence rate from exp MDF
                %Only one activty rate to consider
                alpha_info  = alpha_var_comb((recurrence_para:alpha_inc:length(alpha_var_comb)),3); 
            end
     
            %Occurrence times estimated from 1/alpha rate where alpha rate is
            %the fault's activity rate. In other words, the frequency of events
            %over Mmin, so 1/alpha is the recurrence interval of events
            %>Mmin
        
            for ff=1:length(alpha_info)
                %need inverse of activity rate (ie recurrence interval) to
                %set as mean in expinv
                alpha_info(ff,2)=1/alpha_info(ff,1);
            end
       
            %Poisson model for earthquake occurrence 
            occurrence_time = zeros(1,length(alpha_info));
            %Uses expinv, where event is taken if occurrence time is below 1
             %Needs to be done for each fault
            occurrence_time = expinv(rand(1,length(occurrence_time)),alpha_info(:,2)'); 
        
            if min(occurrence_time) <= t_limit
            
                f_indx = find(occurrence_time <=t_limit);%index faults for which event occurred
                global_count = global_count + length(f_indx);
            
                for rr=1:length(f_indx) %for number of faults with events in this simulation cycle
                    recurrence_info1 = fm_var{f_indx(rr)}(:,3*(recurrence_para-1)+2);
                    mw = mag_range_GR{f_indx(rr)}(find(rand(1,1)>recurrence_info1,1,'last')+1) + rand(1,1)*dm - dm/2;
                    
                    % StochasticEventCatalog
                    % 1    ) Event number
                    % 2    ) Simulation cycle
                    % 3    ) Occurrence time
                    % 4    ) Earthquake magnitude
                    % 5    ) NCZFM ID (if segmented) Combined ID
                    % 6    ) Recurrence type (1 = characteristic versus 2 = GR)
                    % 7    ) Reccurrece parameter
     
                    StochasticEventCatalog{catalog_count}(global_count-length(f_indx)+rr,1:7) = [global_count-length(f_indx)+rr ii*t_limit occurrence_time(f_indx(rr)) mw source_id(f_indx(rr),1) recurrence_type recurrence_para];
                end %end rr loop
            
                num_event_tmp(ii,1) = length(f_indx); %Num_event for each simulation in catalog 
          
            else  %no event occurred in simulation cycle
                  %Record if no events during that simulation
                num_event_tmp(ii,1) = 0;
   
            end
            
        end %end simulation cycle
    
        %Remove unnecessary zeros;
        StochasticEventCatalog{catalog_count}=StochasticEventCatalog{catalog_count}(1:sum(num_event_tmp),:);
        num_event{catalog_count}=num_event_tmp;
    
    end %end ss loop

    %end ss loop
    
    %clear variables for segmented case
    clear alpha_var1 fm_char_var1 fm_exp_var1 RecurrenceVar mag_range_GR
    
end %end ee loop

%% Save results

StochasticEventCatalog_poisson=StochasticEventCatalog;

variable_name=strcat('model',num2str(model_opt),'/catalog_poisson');
save(variable_name,'StochasticEventCatalog_poisson','num_simu');