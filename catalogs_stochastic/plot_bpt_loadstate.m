%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot BPT Load State for two different aperiodicities %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp_duration=500; %length of catalog sampled
color_opt=vertcat([1 0 0],[0 0 1]);

aperiodicity_opt=[0.5,2]; label_opt=["(a)","(b)"]; 
model_opt=[1,2]; fntsize=10;

%Select fault by Combined ID
ff=30; %Combined ID for Dunstan

figure(204)
tiledlayout(2,1);

for mm=1:length(model_opt)

nexttile

    model_path=strcat('model',num2str(model_opt(mm)));
    addpath(model_path)
    load catalog_bpt; load orb_fault_combined
    f_indx=find(orb_faults_comb.Combined_ID==ff); 
    
    for ss=1:2
        plot(1:1:tmp_duration,cat_load_state{ss+1}(1:tmp_duration,f_indx),'Color',color_opt(ss,:),'Linewidth',1.1);hold on 
    end

    plot([1,tmp_duration],[1,1],'k--','Linewidth',1.2); plot([1,tmp_duration],[0,0],'k--','Linewidth',1.2);

    legend('char','G-R','Fontsize',11,'Location','southwest');

    xlabel('time (yrs)'); ylabel('load state'); set(gca,'fontsize',12); ylim([-2.5,1.5]);

    t1=title(strcat(label_opt(mm),' \alpha = ',num2str(aperiodicity_opt(mm))),...
    'fontsize',fntsize+2,'fontweight','normal');
    set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
    h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);

    rmpath(model_path)
end

%set(gcf,'Position',[680 584 936 393]); 
