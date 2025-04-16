
%%%% PLOT COMPARISON OF INTEREVENT TIMES FOR DIFFERENT APERIODICITIES %%%
%%%%    FOR A SINGLE FAULT. USE FOR PLOTTING FIGURE S9 IN MANUSCRIPT %%%%

tmp_limit=600; bylim=0.03; %set based on indivdual fault

%FIGURE ONLY WORKS WITH SELECT_OPT = 1 AND SO FAULT MUST BE SELECTED BY NZCFM ID
select_opt=1;

%Select fault by NZCFM ID
ff=184; %NZCFM ID for Akatore

model_opt=[1,2];

label_opt=vertcat(["(a)","(b)"],["(c)","(d)"]); location_opt=["Northwest","Northeast"];

col_opt=["k","r","m"];%set colours for plotting poisson, bpt, or weibull resultl

figure(203);
tiledlayout(length(model_opt),2,"TileSpacing","compact");

for ii=1:length(model_opt)
    
    nexttile

    %reset path to different model
    model_path=strcat('model',num2str(model_opt(ii)));
    addpath(model_path)
    
    load('orb_fault_segmented','mag_range_GR','fm_char_var1','fm_exp_var1',...
        'RecurrenceVar','T','lambda','orb_faults');
    
    load orb_fault_parameters; f_indx=find(orb_faults.OBJECTID==ff);

    load catalog_poisson; load catalog_bpt; load catalog_weibull

    %index fault events in segmented-char catalog
    flt_ruptures_poi=find(StochasticEventCatalog_poisson{1}(:,5)==ff);
    flt_ruptures_bpt=find(StochasticEventCatalog_bpt{1}(:,5)==ff);
    flt_ruptures_weibull=find(StochasticEventCatalog_weibull{1}(:,5)==ff);

    intertimes_bpt=zeros(length(flt_ruptures_bpt),1); intertimes_poi=zeros(length(flt_ruptures_poi),1);
    intertimes_wbl=zeros(length(flt_ruptures_weibull),1);
    
    %obtain empiricial interevent times from stochastic catalogs
    
    for kk=1:length(flt_ruptures_bpt)-1   
        intertimes_bpt(kk)=StochasticEventCatalog_bpt{1}(flt_ruptures_bpt(kk+1),2)-StochasticEventCatalog_bpt{1}(flt_ruptures_bpt(kk),2);   
    end

    for kk=1:length(flt_ruptures_poi)-1   
        intertimes_poi(kk)=StochasticEventCatalog_poisson{1}(flt_ruptures_poi(kk+1),2)-StochasticEventCatalog_poisson{1}(flt_ruptures_poi(kk),2);   
    end    

    for kk=1:length(flt_ruptures_weibull)-1   
        intertimes_wbl(kk)=StochasticEventCatalog_weibull{1}(flt_ruptures_weibull(kk+1),2)-StochasticEventCatalog_weibull{1}(flt_ruptures_weibull(kk),2);   
    end 
    
    intertimes={}; 
    intertimes={intertimes_poi,intertimes_bpt,intertimes_wbl};
    
    %obtain interevent time distributions from input probability distributions using median recurrence model
    
    exp_cdf=cdf('Exponential',1:tmp_limit,median(T{f_indx}(:,1)));
    exp_pdf=pdf('Exponential',1:tmp_limit,median(T{f_indx}(:,1)));
    
    %cdf of interevent times follows Inverse Gaussian Distribution for BPT, see Eq. 12 of Matthews et al (2002)
    bpt_cdf=cdf('InverseGaussian',1:tmp_limit,median(T{f_indx}(:,1)),median(lambda{f_indx}(:,1)));
    bpt_pdf=pdf('InverseGaussian',1:tmp_limit,median(T{f_indx}(:,1)),median(lambda{f_indx}(:,1)));    
    
    wbl_cdf=cdf('Weibull',1:tmp_limit,tau_wbl{1}(f_indx),beta_wbl{1}(f_indx));
    wbl_pdf=pdf('Weibull',1:tmp_limit,tau_wbl{1}(f_indx),beta_wbl{1}(f_indx));
    
    plot(1:tmp_limit,exp_cdf,col_opt(1),'Linewidth',1.2); hold on  
    plot(1:tmp_limit,bpt_cdf,col_opt(2),'Linewidth',1.2); hold on
    plot(1:tmp_limit,wbl_cdf,col_opt(3),'Linewidth',1.2); hold on 
    
    for dd=1:3
        [f,x]=ecdf(intertimes{dd});plot(x,f,'Color',col_opt(dd),'Linewidth',1.2,'LineStyle','--'); hold on
    end
    
    xlim([0 tmp_limit]); axis square
    xlabel ('interevent time (yrs)');
    
    legend('target poi','target bpt','target wbl',...
            'catalog poi','catalog bpt','catalog wbl','Location','Southeast');  
    
    t1=title(label_opt(ii,1),'fontsize',11,'fontweight','normal');
    set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
    h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

    nexttile %second column, plot empirical and theoritical hazard functions
    
    h_exp=exp_pdf./(1-exp_cdf); %theoritical hazard function following generic eq. 3 of Yakovlev et al (2006) 
    h_bpt=bpt_pdf./(1-bpt_cdf); 
    h_wbl=wbl_pdf./(1-wbl_cdf); 
    
    %plot theoritical hazard functions
    plot(1:tmp_limit,h_exp,col_opt(1),'Linewidth',1.4);hold on
    plot(1:tmp_limit,h_bpt,col_opt(2),'Linewidth',1.4);hold on
    plot(1:tmp_limit,h_wbl,col_opt(3),'Linewidth',1.4);hold on
    
    %derive cumulative hazard function (i.e. integral of hazard function)

    bin=30; %time period over which averaging hazard rate    
    edges=[1:bin:tmp_limit+bin];

    
    %calculate time average hazard rate from cumulative hazard function over
    %time period of bin following: https://www.itl.nist.gov/div898/handbook/apr/section1/apr123.htm
     
    for dd=1:3
            [f,x]=ecdf(intertimes{dd},'function','cumhazard'); 
            afr=zeros(length(edges)-1,1);
            
            for jj=1:(length(edges)-1)
                %find closest index to empirical cumulaltive hazard(t) and t
                [~,tmpindx1]=min(abs(edges(jj)-x));
                [~,tmpindx2]=min(abs(edges(jj+1)-x));
                afr(jj)=(f(tmpindx2)-f(tmpindx1))/bin; %time average hazard rate
            end  
            
            plot(edges(1:end-1),afr,'Color',col_opt(dd),'Linewidth',1.4,'LineStyle','--'); hold on
    end

    xlim([0 tmp_limit]); axis square; ylim([0 bylim])
    xlabel ('time since last event (yrs)'); ylabel('hazard rate');

    t1=title(label_opt(ii,2),'fontsize',11,'fontweight','normal');
    set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
    h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

    legend('target poi','target bpt','target wbl',...
            'catalog poi','catalog bpt','catalog wbl','Location',location_opt(ii));
    
    clear StochasticEventCatalog_poi; clear StochasticEventCatalog_bpt; clear StochasticEventCatalog_wbl
end

set(gcf,'Position',[579 176 652 621]);
