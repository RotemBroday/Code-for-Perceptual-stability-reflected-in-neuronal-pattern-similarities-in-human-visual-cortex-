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

DataPath = ['E:\Rotem\Adaptation Project\Scripts\ScriptsAndData_BrodayDvir_etal2002\PreprocessedData\']; %edit accordingly


%% Load dataset


CurROI = 'ContentSelectiveElecs';
TmpData= struct2cell(load([DataPath, CurROI, '.mat' ]));
ROI_Struct = TmpData{1};


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


TmpSig = sig+p_masked;
SigMat = TmpSig==2;

%% plot correlation time courses

  Selected_Tpoints = [-0.1, 0.3, 0.7, 1.1,1.5,  2];
 [I, i1, ~]=intersect(Time,Selected_Tpoints);
 figure;
 set(gcf,'Position',[100 100 1800 300])
 
 for iTC = 1:numel(I)
     R_TC = GrandMean_rvals_CorrPerItem (i1(iTC),:);
    
      pTC = SigMat(i1(iTC),:);
      
      R_TC_std =nanstd(squeeze(Mean_rvals_CorrPerItem(i1(iTC),:,:)));
   
     subplot(1,6, iTC); 
 hold on;
   h1=shadedErrorBar(Time,R_TC,R_TC_std./sqrt(28),{ 'color', [ 25/255,25/255,112/255],'LineWidth',4,'LineSmoothing','on'},1);

       hold on;
      plot(Time(find(pTC==1)), -0.08*ones(1,numel(find(pTC==1))), 'Color', [1 174/255, 86/255], 'LineWidth', 4)
     title(['Time bin: ' num2str(Selected_Tpoints(iTC)) 's'])
 xlim([-0.5 2.2]);
 ylim([-0.2 1])
 xticks([-0.5:0.5:2])
 yticks([-0.2:0.2: 1]);
 boxplotX = [-0.75,-0.01,0,1.5,1.51,2.3];
boxplotY = [-0.195,-0.195,-0.12,-0.12,-0.195,-0.195];

hold on; plot(boxplotX, boxplotY, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on; plot([Selected_Tpoints(iTC), Selected_Tpoints(iTC)], [-0.2, 1], 'k--');
xlabel('Time [sec]')

 ylabel('Correlation (r)')
 set(gca, 'FontSize', 14);
 box off
 end