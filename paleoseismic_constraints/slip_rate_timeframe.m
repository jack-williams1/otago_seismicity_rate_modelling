%% Plot bar graph showing slip rate timeframe for NZCFM Otago Fault slip rate estimates

close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

orb_faults=readtable('OtagoRangeBasinFaults.xlsx','Sheet','segmented');

srt_otago=zeros(height(orb_faults),1);

for kk=1:height(orb_faults)

    if string(orb_faults.SRT_pref(kk))=="1000s yrs"
        srt_otago(kk)=1000;
    elseif string(orb_faults.SRT_pref(kk))=="18 000"
        srt_otago(kk)=18000;
    elseif string(orb_faults.SRT_gen(kk))=="10,000s yrs"
         srt_otago(kk)=10000;
    elseif string(orb_faults.SRT_gen(kk))=="100,000s yrs"
         srt_otago(kk)=1e5;     
    elseif string(orb_faults.SRT_gen(kk))=="≥ 1 Ma" 
        srt_otago(kk)=1e6;
    end
end

xx=["1000","10000","18000","100,000","≥ 1 Ma"];

yy(1)=length(find(srt_otago==1000));
yy(2)=length(find(srt_otago==10000));
yy(3)=length(find(srt_otago==18000));
yy(4)=length(find(srt_otago==1e5));
yy(5)=length(find(srt_otago==1e6));

figure(1);
bar(xx,yy);
ylabel('count'); xlabel('Starting age (yrs)'); set(gca,'fontsize',13); grid on

set(gcf,'Position',[440 354 626 344]);