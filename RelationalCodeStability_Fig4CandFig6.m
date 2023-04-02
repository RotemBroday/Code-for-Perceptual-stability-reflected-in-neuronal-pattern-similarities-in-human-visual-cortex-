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
% Calculate the RDM (stimuli-pairwise distances matrix) for each time bin
% seperately, Correlations/distances calculated in a leave-1-out- manner across the
%individual repetitions of each stimulus
%store the results in a stimuli*stimuli matrix, for each time bin and
%leave-1-out-iteration- resulting in a grand nstim*nstim*ntime*nrepeats
%matrix
nRepeats=4;
nElec = numel(ROI_Struct.Elecs);
nTimeBins = numel(ROI_Struct.Time);
nStimuli = numel(ROI_Struct.Stimuli);
ROI_Data= ROI_Struct.Data;
Time = ROI_Struct.Time;
PairwiseDistMat = nan(nStimuli, nStimuli, nTimeBins, nRepeats); % mat containing all stimuli pairwise distances
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


CorrPerItem = nan(nTimeBins, nTimeBins, nStimuli, iRepeat);
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

MeanCorrPerItem = nanmean(CorrPerItem,4); %avg over repeats
Mean_rvals_CorrPerItem = nanmean(rvals_CorrPerItem,4);


GrandMean_CorrPerItem = nanmean(MeanCorrPerItem,3);
GrandMean_rvals_CorrPerItem = nanmean(Mean_rvals_CorrPerItem,3);


%% Shuffling 
     


ShuffledPermDistMat =  nan(1000,nStimuli, nStimuli, nTimeBins, nRepeats);
parfor iPerm = 1:1000
    ShuffledPairwiseDistMat = nan(nStimuli, nStimuli, nTimeBins, nRepeats);
    Shuffled_r_values = nan(nStimuli, nStimuli, nTimeBins, nRepeats);
    Shuffled_p_values = nan(nStimuli, nStimuli, nTimeBins, nRepeats);
     Shuffled = randperm(nElec); 
 
for iTime = 1:numel(Time)
    for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
    
        
        StimVec1 = squeeze(ROI_Data(Shuffled, iTime, :, iRepeat));
        
           
            StimVec2 =squeeze(nanmean(ROI_Data(:, iTime,:, AvgOver),4));
            %get 1-pearson correlation
            [r,p] = corr(StimVec1, StimVec2, 'Type', 'Pearson', 'rows', 'complete');
            Dist = 1-r;
            %save in matrix
           ShuffledPairwiseDistMat (:,:,iTime,iRepeat) = Dist;
        
            
        end
 
end
ShuffledPermDistMat(iPerm,:,:,:,:) = ShuffledPairwiseDistMat;
end

%% shuffling continued
ShuffledPermCorrPerItem = nan(1000,nTimeBins, nTimeBins, nStimuli, nRepeats);
parfor iPerm=1:1000
ShuffledCorrPerItem = nan(nTimeBins, nTimeBins, nStimuli, nRepeats);
for iRepeat = 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
    for iStim = 1:nStimuli
    
                
                Vec1 = squeeze(ShuffledPermDistMat(iPerm,iStim,:,:, iRepeat));
                
                Vec2 = squeeze(nanmean(ShuffledPermDistMat(iPerm, iStim,:,:, AvgOver),5));
                [r,p] = corr(Vec1, Vec2, 'Type', 'Pearson', 'rows', 'complete');
                Dist = 1-r;
                ShuffledCorrPerItem(:, :, iStim,iRepeat) = Dist;
                
                

    end
end
ShuffledPermCorrPerItem(iPerm,:,:,:,:) = ShuffledCorrPerItem;
end
 ShuffledMeanCorrPerItem = nanmean(ShuffledPermCorrPerItem,5); %avg over repeats

ShuffledGrandMean_CorrPerItem = nanmean( ShuffledMeanCorrPerItem,4); %over stimuli

%% fdr and cluster corrections 

pMat = nan(numel(Time), numel(Time));
for iTime1 = 1:numel(Time)
    for iTime2 = 1:numel(Time)
       pMat(iTime1,iTime2) = 1-sum(GrandMean_CorrPerItem(iTime1,iTime2)<ShuffledGrandMean_CorrPerItem(:,iTime1,iTime2))/1000;
    end
end
[p_fdr, p_masked] = fdr( pMat, 0.005);

%% cluster correction
clusterTh=1.96; % cluster-defining threshold in zscore units 
I = size(ShuffledGrandMean_CorrPerItem,1);
meanShuffledData = nanmean(ShuffledGrandMean_CorrPerItem(:));

v = nan(I,1);
for i = 1:I
   tmp = squeeze(ShuffledGrandMean_CorrPerItem(i,:,:));
    v(i) = nanvar(tmp(:));
end
MSD = sqrt(nanmean(v));   % mean variance computed across shuffled maps
TH = [meanShuffledData-(clusterTh*MSD) meanShuffledData+(clusterTh*MSD)]; 
shuffledDataThr = ShuffledGrandMean_CorrPerItem<TH(1)|ShuffledGrandMean_CorrPerItem>TH(2);
realDataThr = GrandMean_CorrPerItem<TH(1)|GrandMean_CorrPerItem>TH(2);

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

h=imagesc(GrandMean_CorrPerItem);
% 
TmpMask = sig+p_masked;
MyMask = TmpMask==2;
MyMask=MyMask+0.6;
set(h, 'AlphaData',MyMask);
set(gca, 'XTick', 3:3:28);
set(gca, 'XTickLabels', Time(3:3:28));
set(gca, 'YTick', 3:3:28);
set(gca, 'YTickLabels', Time(3:3:28));
axis square xy;

title({'Relational code stability across time', CurROI}) 
xlabel('Time (s)'); ylabel('Time (s)')
colormap(gca,'parula');
h2=colorbar() ;
caxis([0.3 1]);
ylabel(h2, '1-correlation')
set(gca, 'FontSize', 12);




end