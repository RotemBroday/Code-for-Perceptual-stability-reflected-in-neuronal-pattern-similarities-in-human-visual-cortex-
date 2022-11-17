%% Mean visual responses for visually-responsive contact groups
clear all
close all
clc;
% set the relevant paths:
path_to_toolboxes = 'E:\Rotem\MATLAB_ToolBoxes\';
addpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
rmpath(genpath(fullfile(path_to_toolboxes,'fieldtrip-20161028')));
rand('seed',sum(100*clock));

ElecGroups = {'ContentSelectiveElec', 'EVElec', 'FaceSelectiveElec', 'FrontoParVisElec'};
DataPath = ['E:\Rotem\Adaptation Project\Scripts\ScriptsAndData_BrodayDvir_etal2002\Data_PreSmoothing\']; %edit accordingly
Colors = {[255/255,195/255,10/255], [65/255,95/255,155/255], [178/255,34/255,34/255], [34/255,139/255,34/255]} ;
for iROI = 1:numel(ElecGroups)
clear CurROI ROI_Data MeanPerElectrode

%% Load datasets

CurROI = ElecGroups{iROI};
ROI_Data = struct2cell(load([DataPath, CurROI, 'Data.mat' ]));
ROI_Data = cell2mat(ROI_Data); % dimensions are n Elecs/ n Time points (1500) / n Stimuli (28*4 reps=112)
MeanPerElectrode = nanmean(ROI_Data,3); %Mean across all stimuli
Time = [-750:2:2248];
Time = Time/1000;
MeanPerElectrode= 10*log10(MeanPerElectrode); %log normalization



%% Plotting
figure;
 set(gcf,'Position',[100 100 500 400])
 xlim([-0.5 2.3]);


hold on; plot([0,0], [-0.29,2.8], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 4)
h1=shadedErrorBar(Time,nanmean(MeanPerElectrode,1),nanstd(MeanPerElectrode,0,1)./sqrt(size(MeanPerElectrode,1)),{'color',Colors{iROI} , 'LineWidth',5,'LineSmoothing','on'},1);

xlabel('Time [sec]')
ylabel ('HFB normalized amplitude [dB]')
hold on; plot([0,1.5], [-0.14, -0.14], 'Color', [0.5 0.5 0.5], 'LineWidth', 10)
title(CurROI)
hold on; text(0.2,-0.06, 'Stimulus on', 'Color', [0.3 0.3 0.3], 'FontSize', 16)
set(gca, 'FontSize', 14)

 if iROI==1
 ylim([-0.15 1.2]);
 elseif iROI==2
ylim([-0.15,2.7]);
 elseif iROI==3
ylim([-0.15,1.4]);
 elseif iROI==4
ylim([-0.15,0.4]);
 end
end
