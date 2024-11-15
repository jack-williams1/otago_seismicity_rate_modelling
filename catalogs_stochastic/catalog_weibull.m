%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   SYNTHETIC CATALOG FOR POISSON INTEREVENT TIMES %%%
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

load orb_fault_parameters

disp(['Total duration of simulation: ',num2str(num_simu*t_limit),' years']);

%creation of 3 catalogs:
%1) segmented-char
%2) combined-char
%3) combined-GR

clock_old = cputime; count_case=0;

% Run simulations

catalog_count=1;

for ee=1:2
    
    if ee==1
        load('orb_fault_segmented','T','alpha_bpt','fm_char_var1','fm_exp_var1','RecurrenceVar','orb_faults','mag_range_GR');
        syncat_opt=1; source_id=orb_faults.OBJECTID;
    elseif ee==2
        load('orb_fault_combined','T','alpha_bpt','fm_char_var1','fm_exp_var1','RecurrenceVar','orb_faults_comb','mag_range_GR');
        %syncat opt: both char (1) and GR (2) explored
         syncat_opt=[1,2];  source_id=orb_faults_comb.Combined_ID;
    end
    
    % Not exploring logic tree for Weibull model so only pick median values
    l_pick = round(length(RecurrenceVar{1})/2);
    % Increments between fault values in alpha parameter
    alpha_inc=length(RecurrenceVar{1});
    
    % Combine cell arrays from ea ch fault for use in syncat
    % Last column gives fault id
    alpha_tmp = cell2mat(alpha_bpt');%aperidoicity
    t_tmp = cell2mat(T');%activity rate
    
    for ss=1:length(syncat_opt)
    
        %select on-fault MFD for recurrence model
        if ss==1
            fm_var=fm_char_var1;
        else
            fm_var=fm_exp_var1;
        end
        
        %select median fault activity rates and aperiodicity
        alpha_comb=alpha_tmp((l_pick:alpha_inc:length(t_tmp)),ss);
        t_comb=t_tmp((l_pick:alpha_inc:length(t_tmp)),ss);
        
        num_fault=length(RecurrenceVar);
        tmp_beta_wbf = zeros(num_fault,1); tmp_tau_wbf = zeros(num_fault,1);
        eq_times=cell(num_fault,1);%cell to store eq_times
        tmp_count = 1; 
        
        %initial guesses needed for model opts 2,3, & 4
        if model_opt==3
            init_guess=2.55; %alpha=4
        else
            init_guess=3.15; %alpha=2
        end
    
        for ii=1:num_fault
            
            % Given mean and CoV of each fault derive Weibull shape parameters
            mu = round(t_comb(ii)); % Your mean value;
            sigma = round(t_comb(ii)*alpha_comb(ii)); % Your standard deviation value;
            
            if model_opt==1 %if cov<1.2, can use methods of moment to derive Weibull parameters

                % Define the equations for the moments in terms of shape and scale parameters
                equation1 = @(k, lambda) mu - lambda * gamma(1 + 1/k);
                equation2 = @(k, lambda) sigma - lambda * sqrt(gamma(1 + 2/k) - (gamma(1 + 1/k))^2);

                % Initial guess for shape (k) and scale (lambda) parameters
                
                if ss==3 && ii==29 || ss==2 && ii==21
                    initial_guess_k = 1;
                    initial_guess_lambda = 1;
                else
                    initial_guess_k = 2;
                    initial_guess_lambda = 2;
                end

                % Solve the equations using fsolve
                options = optimoptions('fsolve', 'Display', 'off');
                params = fsolve(@(x) [equation1(x(1), x(2)); equation2(x(1), x(2))], [initial_guess_k, initial_guess_lambda], options);

                % Estimated shape (beta) and scale (tau) parameters
                tmp_beta_wbf(ii) = params(1);
                tmp_tau_wbf(ii)= params(2);
                
            else %if cov>1.2, can use transcendatal equations to derive Weibull parameters
                
                % Define the transcendental equation for solving for the shape parameter (k)
                objective = @(k) abs(sigma - (gamma(1 + 2/k) - (gamma(1 + 1/k))^2)^0.5 / (gamma(1 + 1/k) - gamma(1 + 2/k))^(1.5) * mu);

                % Initial guess for k (you may need to adjust this based on your data)
                initialGuess = init_guess;

                % Optimize the objective function using fminunc
                options = optimoptions('fminunc', 'Display', 'off'); %set 'Display' to 'iter' to show results
                tmp_beta_wbf(ii) = fminunc(objective, initialGuess, options);

                % Estimate scale parameter lambda using the method of moments
                tmp_tau_wbf(ii) = mu / gamma(1 + 1/tmp_beta_wbf(ii));
            end
            
            eq_times{ii}=zeros(num_simu,5); jj=0; count=0;
            
            %create random earthquake times for as long as total interevent times is <num_simu
            % catalog created for 10% longer than num_simu to reduce edge effects (next step)
            
            while sum(eq_times{ii}(:,1))<num_simu+(num_simu*0.1) 
                  jj=jj+1; count=count+1;
                  tmp=wblrnd(tmp_tau_wbf(ii),tmp_beta_wbf(ii));%derive random earthquake interevent time from weibull
                  eq_times{ii}(jj,1)=tmp;
                  eq_times{ii}(jj+1,2)=eq_times{ii}(jj,1)+eq_times{ii}(jj,2);%add interevent times
            end
            
            %randomly sample num_simu window so not all earthquake records start at 0
            eq_times{ii}= eq_times{ii}(1:count,:);
            tmp_rand=randi([0,num_simu*0.1],1);
            tmp_indx=find(eq_times{ii}(:,2)>=tmp_rand & eq_times{ii}(:,2)<(tmp_rand+num_simu)); 
            eq_times{ii}= eq_times{ii}(tmp_indx,:);
            count=length(tmp_indx);%reset count to downsampled catalog
             
            %assign magnitude to each event by randomly indexing earthquake magnitude from fm_var cdf
            recurrence_info1 = fm_var{ii}(:,3*(l_pick-1)+2);%select median reccurence model cdf
            mag_rnd=rand(1,count);
            [~,mag_indx]=min((mag_rnd-recurrence_info1).^2);
            eq_times{ii}(:,3) = mag_range_GR{ii}(mag_indx)+ rand(1,1)*dm - dm/2;
            
            %add source id and recurrence parameter to each event
            eq_times{ii}(:,4:5)=[source_id(ii).*ones(count,1) ss.*ones(count,1)];
            
            if ii == tmp_count*5
                tmp_count = tmp_count + 1;
                clock_new = cputime;
                disp(['Current fault is ',num2str(ii),' out of ',num2str(num_fault),' & Required time is: ',...
                num2str(clock_new-clock_old), '. Catalog is ', num2str(catalog_count), ' out of 3']);
                clock_old = clock_new;
            end
 
        end %end ii loop for each fault
                    
         % StochasticEventCatalog
         % 1    ) Event number
         % 2    ) Simulation cycle
         % 3    ) Interevent time
         % 4    ) Earthquake magnitude
         % 5    ) NCZFM ID (if segmented) Combined ID
         % 6    ) Recurrence type (1 = characteristic versus 2 = GR)
        
        tmp_catalog=sortrows((cell2mat(eq_times)),2);
        event_id=[1:1:length(tmp_catalog)]'; 
        StochasticEventCatalog{catalog_count}= [event_id tmp_catalog(:,2) tmp_catalog(:,1) tmp_catalog(:,3:5)];
        
        tau_wbl{catalog_count}=tmp_tau_wbf; beta_wbl{catalog_count}=tmp_beta_wbf; 
        catalog_count=catalog_count+1;
        
           
    end %end ss loop for each catalog
          
    
    %clear variables for segmented case
    clear alpha_var1 fm_char_var1 fm_exp_var1 RecurrenceVar mag_range_GR
    
    
end %end ee loop

%% Save results, include alpha and beta parameters of wbl distributions

StochasticEventCatalog_weibull=StochasticEventCatalog;
variable_name=strcat('model',num2str(model_opt),'/catalog_weibull');

save(variable_name,'StochasticEventCatalog_weibull','num_simu','tau_wbl','beta_wbl');