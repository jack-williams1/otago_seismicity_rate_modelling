%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the number of empty counts in 70-year catalog subsamples %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


model_opt=[1,2,3];
intertime_txt=["poi","bpt","wbl"];
tiledlayout(1,3);
aperiodicity_opt=[0.5,2,4]; label_opt=["(a)","(b)","(c)"]; 

fntsize=7;

for ii=1:length(model_opt)
   
    nexttile

    model_path=strcat('model',num2str(model_opt(ii)));
    addpath(model_path)
    
    empty_count=zeros(length(intertime_txt),3);

    for jj=1:3  %loop through renewal processes
        load(strcat('catalog_',intertime_txt(jj),'_statistics'),strcat('sampleMoRate_',intertime_txt(jj)));
        samples=eval(strcat('sampleMoRate_',intertime_txt(jj)));
        
        for kk=1:3 %loop through interevent times
            empty_sample=find(samples(:,kk)==0); 
            empty_count(jj,kk)=100*length(empty_sample)/1000; 
            
            
        end
    end
    
    bar_cat=["seg-char poi","seg-char bpt","seg-gr wbl",...
                          "comb-char poi","comb-char bpt","comb-char wbl",...
                          "comb-gr poi","comb-gr bpt","comb-gr wbl"]';
        
    b=bar([empty_count(:,1)',empty_count(:,2)',empty_count(:,3)'],...
        'FaceColor','flat','FaceAlpha', 0.3);hold on
    b.CData([1,4,7],:) = vertcat([0 0 0],[0 0 0],[0 0 0]);
    b.CData([2,5,8],:) = vertcat([1 0 0],[1 0 0],[1 0 0]);
    b.CData([3,6,9],:) = vertcat([1 0 1],[1 0 1],[1 0 1]);

    tmp_xlim=get(gca,'Xlim'); ylim([0 50]); set(gca,'Fontsize',fntsize); axis square
    plot([tmp_xlim(1),tmp_xlim(2)],[5,5],'b--','LineWidth',1.5);
    xticklabels(gca,bar_cat); xtickangle(45); ylabel('empty counts (%)');


    t1=title(strcat(label_opt(ii),' \alpha = ',num2str(aperiodicity_opt(ii))),...
    'fontsize',fntsize+2,'fontweight','normal');
    set(t1, 'horizontalAlignment', 'left'); set(t1, 'units', 'normalized');
    h1 = get(t1, 'position'); set(t1, 'position', [-0.15 h1(2) h1(3)]);
    
    rmpath(model_path)
end

set(gcf,'Position',[440 495 581 203]);