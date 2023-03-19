%calculate the population pattern stability across time, averaged across
%all stimuli, for the different ROI electrode groups (figure 4 a and figure 6 left panels- pattern stability)

clear all
close all
clc;
% set the relevant paths:
path_to_toolboxes = 'E:\Rotem\MATLAB_ToolBoxes\';
addpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
rmpath(genpath(fullfile(path_to_toolboxes,'fieldtrip-20161028')));
rand('seed',sum(100*clock));

ElecGroups = {'ContentSelectiveElecs', 'FaceSelectiveElecs', 'EVElecs'};
DataPath = ['E:\Rotem\Adaptation Project\Scripts\ScriptsAndData_BrodayDvir_etal2002\PreprocessedData\']; %edit accordingly

for iROI = 1:numel(ElecGroups)
%% Load datasets

clear CurROI ROI_Data TmpData ROI_Struct
CurROI = ElecGroups{iROI};
TmpData= struct2cell(load([DataPath, CurROI, '.mat' ]));
ROI_Struct = TmpData{1};

%%
%calculate the distances (1-correlation) between all time bins of each
%stimulus related pattern
%Correlations/distances calculated in a leave-1-out- manner across the
%individual repetitions of each stimulus
%store the results in a time*time (DS)*stimuli matrix
nRepeats=4;
nElec = numel(ROI_Struct.Elecs);
nTimeBins = numel(ROI_Struct.Time);
nStimuli = numel(ROI_Struct.Stimuli);
ROI_Data= ROI_Struct.Data;


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
AvgDistPerStimMat = nanmean(DistPerStimMat,4); %mean across 4 leave-1-out iterations
Avg_r_values = nanmean(r_values,4);
Avg_p_values = nanmean(p_values,4);

%% 3. shuffle across 1000 permutations

ShuffledAvgDistPerStimMat = nan(1000,nTimeBins, nTimeBins, nStimuli);
ShuffledAvg_r_values=nan(1000,nTimeBins, nTimeBins, nStimuli);
ElecInd = 1:nElec;
parfor iPer = 1:1000
ShuffledDistPerStimMat = nan(nTimeBins, nTimeBins, nStimuli, nRepeats);
Shuffled_r_values = nan(nTimeBins, nTimeBins, nStimuli, nRepeats);
Shuffled = ElecInd (randperm(nElec));
%loop over stimuli and repetitions
for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
for iStim = 1:nStimuli
  
        tMat1 = squeeze(ROI_Data(Shuffled , :, iStim, iRepeat));

            tMat2 = squeeze(nanmean(ROI_Data(:, :, iStim,AvgOver),4));
        %get 1-pearson correlation
        [r,p] = corr(tMat1, tMat2, 'Type', 'Pearson', 'rows', 'complete');
        Dist = 1-r;
        %save in matrix
        ShuffledDistPerStimMat(:, :, iStim, iRepeat) = Dist;
        Shuffled_r_values (:, :, iStim, iRepeat) = r;
      

end
end
ShuffledAvgDistPerStimMat(iPer,:,:,:) = nanmean(ShuffledDistPerStimMat,4);
ShuffledAvg_r_values(iPer,:,:,:) = nanmean(Shuffled_r_values,4);

end


 
%% fdr and cluster corrections 
Time = ROI_Struct.Time;
MeanDistMat = nanmean(AvgDistPerStimMat(:,:,:),3);
MeanShuffledAvgDistMat = nanmean(ShuffledAvgDistPerStimMat(:,:,:,:),4);

pMat = nan(numel(Time), numel(Time));
for iTime1 = 1:numel(Time)
    for iTime2 = 1:numel(Time)
       pMat(iTime1,iTime2) = 1-sum(MeanDistMat(iTime1,iTime2)<MeanShuffledAvgDistMat(:,iTime1,iTime2))/1000;
    end
end

[p_fdr, p_maskedAll] = fdr( pMat, 0.001);

%% cluster correction
clusterTh=1.96; % cluster-defining threshold in zscore units (p<0.05)
I = size(ShuffledAvgDistPerStimMat,1);
meanShuffledData = nanmean(ShuffledAvgDistPerStimMat(:));

v = nan(I,1);
for i = 1:I
    tmp = squeeze(ShuffledAvgDistPerStimMat(i,:,:));
    v(i) = nanvar(tmp(:));
end
MSD = sqrt(nanmean(v));   % mean variance computed across shuffled maps
TH = [meanShuffledData-(clusterTh*MSD) meanShuffledData+(clusterTh*MSD)]; 
shuffledDataThr = ShuffledAvgDistPerStimMat<TH(1)|ShuffledAvgDistPerStimMat>TH(2);
realDataThr = MeanDistMat<TH(1)|MeanDistMat>TH(2);

MaxClusterSize = [];
for k = 1:size(shuffledDataThr,1)
    tmp = squeeze(shuffledDataThr(k,:,:));
    MaxClusterSize(k) = sum(sum(bwareafilt(tmp,1))); % returns the size of the largest cluster
end

[ll,C] = bwlabel(realDataThr,8);
clusterSize = []; 
P_clusters = [];
for k = 1:C
    clusterSize(k) = sum(sum(ll==k));
    currentPvalue = (sum(MaxClusterSize > clusterSize(k)) + 1) / (length(MaxClusterSize)+1); % see Davison and Hinkley (1997) formula
    P_clusters(k) = currentPvalue;
    if currentPvalue < 0.05
        sig = ll==k;
 
    end
end

%% Plot results
figure; set(gcf,'Position',[400 400 700 700])

h3=imagesc(MeanDistMat);
% 
TmpMask = sig+p_maskedAll;
MyMask = TmpMask==2;
MyMask=MyMask+0.6;
set(h3, 'AlphaData',MyMask);
set(gca, 'XTick', 3:3:28);
set(gca, 'XTickLabels', Time(3:3:28));
set(gca, 'YTick', 3:3:28);
set(gca, 'YTickLabels', Time(3:3:28));
axis square xy;

title({'Average pattern stability across time', CurROI}) 
xlabel('Time (s)'); ylabel('Time (s)')
colormap(gca,'parula');
h2=colorbar() ;
caxis([0.3 1.1]);
ylabel(h2, '1-correlation')
set(gca, 'FontSize', 12);

end