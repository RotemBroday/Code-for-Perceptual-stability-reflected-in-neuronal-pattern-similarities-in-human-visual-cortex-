clear all
close all
clc;
% set the relevant paths:
path_to_toolboxes = 'E:\Rotem\MATLAB_ToolBoxes\';
addpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
rmpath(genpath(fullfile(path_to_toolboxes,'fieldtrip-20161028')));
rand('seed',sum(100*clock));


DataPath = ['E:\Rotem\Adaptation Project\Scripts\ScriptsAndData_BrodayDvir_etal2002\PreprocessedData\']; 
TmpData= struct2cell(load([DataPath, 'ContentSelectiveElecs.mat' ]));
ROI_Struct = TmpData{1};

%%
%comparing the amplitude, PV, RDM and relational code differences between
%time 0.2 and 1 sec
nRepeats=4;
nElec = numel(ROI_Struct.Elecs);
nTimeBins = numel(ROI_Struct.Time);
nStimuli = numel(ROI_Struct.Stimuli);
ROI_Data= ROI_Struct.Data;


RefTimeEarly = find(ROI_Struct.Time==0.2);
RefTimeLate = find(ROI_Struct.Time==1.2);

%% Get mean amplitude ratio 
MeanAcrossElecs = squeeze(nanmean(ROI_Data,1));
MeanAcrossReps = nanmean(MeanAcrossElecs,3);
MeanAcrossStim = squeeze(nanmean(MeanAcrossReps ,2));



EarlyVsLateRatio_Amplitude =MeanAcrossReps(RefTimeEarly,:)./MeanAcrossReps(RefTimeLate,:);
[EarlyVsLateRatio_Amplitude, I1] = rmoutliers(EarlyVsLateRatio_Amplitude);

[h, pAmp, ~, statsAmp] = ttest2(MeanAcrossReps(RefTimeEarly,:), MeanAcrossReps(RefTimeLate,:));



%% Get PV mean ratio 
DistPerStimMat = nan(nTimeBins, nTimeBins, nStimuli, nRepeats);
r_values = nan(nTimeBins, nTimeBins, nStimuli, nRepeats);
p_values = nan(nTimeBins, nTimeBins, nStimuli, nRepeats);

%loop over stimuli and individual repetitions
for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
for iStim = 1:nStimuli
    %loop over 1st dimension time bin
    for iTime1 = 1:nTimeBins
        clear tVec1
        tVec1 = ROI_Data(:, iTime1, iStim, iRepeat);
        %loop over 2nd dimension time bin
        for iTime2 = 1:nTimeBins
            clear tVec2
            tVec2 = nanmean(ROI_Data(:, iTime2, iStim,AvgOver),4);
        %get 1-pearson correlation
        [r,p] = corr(tVec1, tVec2, 'Type', 'Pearson', 'rows', 'complete');
        Dist = 1-r;
        %save in matrix
        DistPerStimMat(iTime1, iTime2, iStim, iRepeat) = Dist;
        r_values(iTime1, iTime2, iStim, iRepeat) = r;
        p_values(iTime1, iTime2, iStim, iRepeat) = p;
        end
    end
end
end

PV_Dist_MeanAcrossReps =squeeze(nanmean(DistPerStimMat,4));
PV_rVals_MeanAcrossReps = squeeze(nanmean(r_values,4));
RefTimeStimEnd = find(ROI_Struct.Time==1.5);


PV_rVals_Early =squeeze(PV_rVals_MeanAcrossReps(RefTimeStimEnd,RefTimeEarly,:));
PV_rVals_Late =squeeze(PV_rVals_MeanAcrossReps(RefTimeStimEnd,RefTimeLate,:));

clear I1 I2
EarlyVsLateRatio_PV = PV_rVals_Early ./PV_rVals_Late;
[EarlyVsLateRatio_PV, I1] = rmoutliers(EarlyVsLateRatio_PV );
[h, pPV, ~, statsPV] = ttest2(PV_rVals_Early, PV_rVals_Late);




%% relational code ratio

PairwiseDistMat = nan(nStimuli, nStimuli, nTimeBins, nRepeats);
r_values = nan(nStimuli, nStimuli, nTimeBins, nRepeats);
p_values = nan(nStimuli, nStimuli, nTimeBins, nRepeats);

for iTime = 1:nTimeBins
    for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
    for iStim1=1:nStimuli
        clear StimVec1
        StimVec1 = ROI_Data(:, iTime, iStim1, iRepeat);
        for iStim2 = 1:nStimuli
            clear StimVec2
            StimVec2 =nanmean(ROI_Data(:, iTime, iStim2, AvgOver),4);
            %get 1-pearson correlation
            [r,p] = corr(StimVec1, StimVec2, 'Type', 'Pearson', 'rows', 'complete');
            Dist = 1-r;
            %save in matrix
            PairwiseDistMat (iStim1,iStim2,iTime,iRepeat) = Dist;
            r_values(iStim1,iStim2,iTime,iRepeat) = r;
            p_values(iStim1,iStim2,iTime,iRepeat) = p;
            
        end
    end 
    end
end



%% Correlation relations code (distances vector) per stimuli across time bins
CorrPerItem = nan(nTimeBins, nTimeBins, nStimuli, iRepeat);
rvals_CorrPerItem = nan(nTimeBins, nTimeBins, nStimuli, iRepeat);
for iRepeat = 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
    for iStim = 1:nStimuli
        for iTime1 = 1:nTimeBins
            for iTime2 = 1:nTimeBins
                clear Vec1
                Vec1 = squeeze(PairwiseDistMat(iStim,:,iTime1, iRepeat));
                clear Vec2
                Vec2 = squeeze(nanmean(PairwiseDistMat(iStim,:,iTime2, AvgOver),4));
                [r,p] = corr(Vec1', Vec2');
                Dist = 1-r;
                CorrPerItem(iTime1, iTime2, iStim,iRepeat) = Dist;
                rvals_CorrPerItem(iTime1, iTime2, iStim,iRepeat)=r;
                
            end
        end
    end
end

 MeanCorrPerItem = squeeze(nanmean(CorrPerItem,4)); %avg over reps
 Mean_rvals_CorrPerItem = squeeze(nanmean(rvals_CorrPerItem,4));



RC_rVals_Early =squeeze(Mean_rvals_CorrPerItem(RefTimeStimEnd,RefTimeEarly,:));
RC_rVals_Late =squeeze(Mean_rvals_CorrPerItem(RefTimeStimEnd,RefTimeLate,:));

clear I1 I2
EarlyVsLateRatio_RC = RC_rVals_Early ./RC_rVals_Late;
[EarlyVsLateRatio_RC, I1] = rmoutliers(EarlyVsLateRatio_RC);


[h, pRC, ~, statsRC] = ttest2(RC_rVals_Early, RC_rVals_Late);



%% plot

% 1. plot ratio graph
X = categorical({'Amplitude','PV','RC'});
X = reordercats(X,{'Amplitude','PV','RC'});
DataToPlot = [mean(EarlyVsLateRatio_Amplitude), mean(EarlyVsLateRatio_PV),  mean(EarlyVsLateRatio_RC)];
figure;  set(gcf,'Position',[200 200 200 400])

b=bar (X, DataToPlot);
b.FaceColor = 'flat';
b.CData(1,:) = [178/255,34/255,34/255];
b.CData(2,:) = [28/255,99/255,124/25];
b.CData(3,:) = [62/255,157/255,129/255];
ylim([-1, 18])
hold on;
plot(repmat(X(1),1,numel(EarlyVsLateRatio_Amplitude)), EarlyVsLateRatio_Amplitude, '.', 'MarkerSize', 11, 'Color', [0.5, 0.5, 0.5]);
hold on;
plot(repmat(X(2),1,numel(EarlyVsLateRatio_PV)), EarlyVsLateRatio_PV, '.', 'MarkerSize', 11, 'Color', [0.5, 0.5, 0.5])
hold on;
plot(repmat(X(3),1,numel(EarlyVsLateRatio_RC)), EarlyVsLateRatio_RC, '.', 'MarkerSize', 11, 'Color', [0.5, 0.5, 0.5])

title('Early Vs Late Time Ratio');


 %% run anova
 DataForAnova = [EarlyVsLateRatio_Amplitude'; EarlyVsLateRatio_PV; EarlyVsLateRatio_RC];
 GroupLabels = [ones(numel(EarlyVsLateRatio_Amplitude),1); 2*ones(numel(EarlyVsLateRatio_PV),1); 3*ones(numel(EarlyVsLateRatio_RC),1)];
 
 
 [p, F, df1, df2] = wanova (DataForAnova,  GroupLabels);
 [~, p_AmpVsPV, ~, stats_AmpVsPV] = ttest2(EarlyVsLateRatio_Amplitude',EarlyVsLateRatio_PV );
  [~, p_AmpVsRC, ~, stats_AmpVsRC] = ttest2(EarlyVsLateRatio_Amplitude',EarlyVsLateRatio_RC );
  [~, p_PVVsRC, ~, stats_PVVsRC] = ttest2(EarlyVsLateRatio_PV ,EarlyVsLateRatio_RC );
  
  

    


