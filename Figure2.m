%% Figure 2 - example activity patterns of 4 stimuli (in refernce to Bruce Willis)
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
    'Leonardo DiCaprio', 'Seinfield', 'Oprah Winfery', 'Dr Evil', 'Donald Duck', 'Phoebe', 'Schawrtznager'...
    'Bill Clinton', 'Pisa', 'Walk of Fame', 'Golden gate bridge', 'World trade center', 'Giza Pyramid', ...
    'Eiffel tower', 'Hoover dam', 'Times square', 'The capital', 'Chicago bean', 'Taj Mahal', 'White House', ...
    'Lion King Cliff', 'Big Ben' };


%% Load datasets

clear CurROI ROI_Data TmpData ROI_Struct
CurROI = ElecGroups{1};
TmpData= struct2cell(load([DataPath, CurROI, '.mat' ]));
ROI_Struct = TmpData{1};

nRepeats=4;
nElec = numel(ROI_Struct.Elecs);
nTimeBins = numel(ROI_Struct.Time);
nStimuli = numel(ROI_Struct.Stimuli);
ROI_Data= ROI_Struct.Data;
MeanAcrossRepsData = nanmean(ROI_Data,4);
Time = ROI_Struct.Time;

%% Get Example stimuli indices
MainStim = 2; %Bruce
TimeWindow = [0.3]; %take example time 0.3s after stim onset
%%
GClr = linspace(0,1,nElec);
BClr = linspace(0,1,nElec);

for iElec=1:nElec
clr(iElec,:)  = [0, GClr(iElec), BClr(iElec)];

end
IndTime= find(Time==TimeWindow);

[B I] = sort(MeanAcrossRepsData(:,IndTime,MainStim ), 'descend');
%%
figure;
set(gcf,'Position',[200 200 800 500])

counter=1;
StimNames = {'Bruce',  'Uma', 'Donald', 'Times Square'};
MainStimResp = MeanAcrossRepsData(I,IndTime,MainStim);
for iStim =[ 2,1,11, 22] % (indices of chosen example stim)
   
     tmpStimResp = MeanAcrossRepsData(I,IndTime,iStim);
     %get elec bar colors sorted according to amp
     [tmp, StimSort] = sort(tmpStimResp);

     for Ind = 1:nElec
    StimClr(StimSort(Ind),:) = clr(Ind,:);
   
     end
 [r,p] = corr(MainStimResp, tmpStimResp);
    NormResp =(tmpStimResp-min(tmpStimResp))/(max(tmpStimResp)-min(tmpStimResp));

subplot(2,2,counter)

b = bar(NormResp, 'facecolor', 'flat','EdgeColor','none','BarWidth',0.9);
b.CData = StimClr;
ylim([0 1.1])
% xlabel('Elec #');
% ylabel('Normalized HFB amplitude')
set(gca, 'FontSize', 12)
title([StimNames{counter} ' r=' num2str(r)])
counter=counter+1;
box off
end


