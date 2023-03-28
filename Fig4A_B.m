%% Relational code histograms for figure 4a and b

clear all
close all
clc;
% set the relevant paths:
path_to_toolboxes = 'E:\Rotem\MATLAB_ToolBoxes\';
addpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
rmpath(genpath(fullfile(path_to_toolboxes,'fieldtrip-20161028')));
rand('seed',sum(100*clock));

ElecGroups = {'ContentSelectiveElecs'};
DataPath = ['E:\Rotem\Adaptation Project\Scripts\ScriptsAndData_BrodayDvir_etal2002\PreprocessedData\']; %edit accordingly

StimLabels = {'Uma Thurman', 'Bruce Willis', 'Britney Spears', 'Aladdin', 'Barak Obama', 'Bart Simpson', ...
    'Leonardo DiCaprio', 'Seinfield', 'Oprah Winfery', 'Dr Evil', 'Donald Duck', 'Phoebe', 'Arnold Schawrtznager'...
    'Bill Clinton', 'Pisa', 'Walk of Fame', 'Golden gate bridge', 'World trade center', 'Giza Pyramid', ...
    'Eiffel tower', 'Hoover dam', 'Times square', 'The capital', 'Chicago bean', 'Taj Mahal', 'White House', ...
    'Lion King Cliff', 'Big Ben' };



TmpData= struct2cell(load([DataPath, 'ContentSelectiveElecs', '.mat' ]));
ROI_Struct = TmpData{1};       
nStimuli = numel(ROI_Struct.Stimuli);

MeanAllStimMat = nanmean(ROI_Struct.Data,4);



%% Bruce willis relational code example for bruce willis
MainStim = 2; %Bruce
Time = ROI_Struct.Time;

Ind3 = find(Time==0.3); % example time bin 0.3s after onset (~peak response)

MainStimResp = MeanAllStimMat(:,Ind3,MainStim);
Dist = nan(1, nStimuli);
for iStim = 1:(nStimuli)
    CurrStimResp = MeanAllStimMat(:,Ind3,iStim);
   clear r
    [r,p] = corr(MainStimResp ,  CurrStimResp);
    Dist(iStim) = 1-r;
end
Dist(MainStim)=[];
[DistsSorted IndDists] = sort(Dist, 'ascend');
StimLabels(MainStim)=[];

GClr = linspace(0,0.7,27);
BClr = linspace(0,1,27);
RClr = linspace(0, 0.3, 27);

for iStim =1:13
Stimclr(iStim,:)  = [RClr(10), GClr(10), BClr(10)];

end

for iStim =14:(nStimuli-1)
Stimclr(iStim,:)  = [RClr(end), GClr(end), BClr(end)];

end

figure;
set(gcf,'Position',[200 200 400 450]);

b = bar(DistsSorted, 'facecolor', 'flat','EdgeColor','none','BarWidth',0.9);
b.CData = Stimclr;
xticks(1:27);
xticklabels(StimLabels(IndDists));
xtickangle(90)
set(gca, 'FontSize', 10);
text(1,0.6, 'Faces', 'Color', [RClr(10), GClr(10), BClr(10)], 'FontSize', 12);
text(1,0.54, 'Places', 'Color', [RClr(end), GClr(end), BClr(end)], 'FontSize', 12);
ylabel('Distance 1-r');
ylim([0 0.8]);
title('Distances from Bruce image')

 

%% fig 4 b- average across sorted relational codes of all stimuli

IndTimeStart = find(Time==0 );
IndTimeEnd = find(Time==1.5 );
StimCategory = {'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F'...
    'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P'};
SortedDists=nan(28,27);
CategoryMat = nan(28,27);
for iStim = 1:nStimuli
    CurMainStimResp = nanmean(MeanAllStimMat(:,[IndTimeStart:IndTimeEnd],iStim),2);
    TmpDist=[];
    for jStim = 1:nStimuli
        CurrOtherStimResp = nanmean(MeanAllStimMat(:,[IndTimeStart:IndTimeEnd],jStim),2);
        [r,p] = corr(CurMainStimResp,CurrOtherStimResp);
        TmpDist(jStim) = 1-r;
    end
[SortedTmpDist, SortedInd ]= sort(TmpDist);
SortedTmpDist(1)=[];
SortedInd(1)=[];

SortedDists(iStim,:) = SortedTmpDist; 


 CategoryMat(iStim,:) =  strcmpi(StimCategory{iStim}, StimCategory(SortedInd)); 
 
end



MeanDists = mean(SortedDists);
StdDists = std(SortedDists);
StdDists = StdDists/(sqrt(28));
NumWithin = sum(CategoryMat);
ProportWithin = NumWithin/28;
WithinBars = ProportWithin.*MeanDists;
BetweenBars = MeanDists-WithinBars;

GClr = linspace(0.3,1,13);
BClr = linspace(0.2,0.6,13);
RClr = linspace(0.2, 0.3, 13);
clear Stimclr
for iStim =1:13
Stimclr(iStim,:)  = [RClr(iStim), GClr(iStim), BClr(iStim)];

end

Stimclr(14,:) = [0.2, 0.2, 0.2];


RClr = linspace(85/255,231/255,13);
GClr = linspace(37/255,208/255,13);
BClr = linspace(134/255, 245/255, 13);
iColor=1;
for iStim =15:27
Stimclr(iStim,:)  = [RClr(iColor), GClr(iColor), BClr(iColor)];
iColor=iColor+1;
end
y = [WithinBars', BetweenBars'];

figure;
set(gcf,'Position',[200 200 400 400]);
b = bar(y,'stacked', 'facecolor', 'flat','EdgeColor','none','BarWidth',0.9);
% b = bar(MeanDists, 'facecolor', 'flat','EdgeColor','none','BarWidth',0.9);
b(1).CData = [0.8500 0.3250 0.0980];
b(2).CData = [0.9290 0.6940 0.1250];
text(1,0.75, 'Within Category Distances', 'Color', [0.8500 0.3250 0.0980], 'FontSize', 12);
text(1,0.7, 'Between Category Distances', 'Color', [0.9290 0.6940 0.1250], 'FontSize', 12);
hold on;
errorbar(MeanDists, StdDists,'k','linestyle','none');
title('Average sorted relational codes')
ylabel('Distances (1-r)');
xlabel('Stimuli (sorted)')
set(gca, 'FontSize', 10);

%% plot mean within vs betweewn category dists bar
WithinDists = SortedDists(CategoryMat==1);
MeanWithingDists = mean(WithinDists);
WithinDistsSD = std(WithinDists);

BetweenDists = SortedDists(CategoryMat==0);
MeanBetweenDists = mean(BetweenDists);
BetweenDistsSD = std(BetweenDists);

figure;
set(gcf,'Position',[200 200 300 400])
b = bar([MeanWithingDists,MeanBetweenDists] , 'facecolor', 'flat','EdgeColor','none','BarWidth',0.9);
hold on;
errorbar([MeanWithingDists,MeanBetweenDists], [WithinDistsSD/sqrt(28),BetweenDistsSD/sqrt(28)],'k','linestyle','none');

b.CData = [0.8500 0.3250 0.0980;0.9290 0.6940 0.1250 ];
set(gca,'XTickLabel',{'Within ', 'Between ' });

ylabel('Distances (1-r)');

set(gca, 'FontSize', 11);