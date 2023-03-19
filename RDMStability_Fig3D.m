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

ElecGroups = {'ContentSelectiveElecs'};
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

%% Divide RDMs to upper and lower half (later averaging across the results of both halves)
%Upper half
for iTime = 1:nTimeBins
    for iRepeat = 1:nRepeats
    tmpVec=[];
    c=0;
    for iRow = 1:nStimuli
        c=c+1;
        tmpVec = [tmpVec PairwiseDistMat(iRow, c:28, iTime,iRepeat)];
    end
    UpperDistVec(iTime, iRepeat,:) = tmpVec;
end
end
 %Lower half
for iTime = 1:nTimeBins
    for iRepeat = 1:nRepeats
    tmpVec=[];
    c=0;
    for iRow = 1:nStimuli
        c=c+1;
        tmpVec = [tmpVec ;PairwiseDistMat(c:28, iRow, iTime, iRepeat)];
    end
    LowerDistVec(iTime,iRepeat,:) = tmpVec';
    end
end

%% get distance matrices between RDMs from different time bins
UpperTimeDistMat = nan(nTimeBins, nTimeBins, nRepeats);
Upper_rvals_TimeDistMat = nan(nTimeBins, nTimeBins, nRepeats);
for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
    for iTime1 = 1:nTimeBins
        clear Vec1
        Vec1 = squeeze(UpperDistVec(iTime1, iRepeat,:));
        for iTime2 = 1:nTimeBins
            clear Vec2
            Vec2 = squeeze(nanmean(UpperDistVec(iTime2,AvgOver, :),2));
            [r,p] = corr(Vec1, Vec2);
            Dist = 1-r;
            UpperTimeDistMat(iTime1, iTime2, iRepeat) = Dist;
            Upper_rvals_TimeDistMat(iTime1, iTime2, iRepeat)=r;
        end
    end
end
UpperAvgTimeDistMat = nanmean(UpperTimeDistMat,3);
UpperAvg_rvals_TimeDistMat = nanmean(Upper_rvals_TimeDistMat,3);


LowerTimeDistMat = nan(nTimeBins, nTimeBins, nRepeats);
Lower_rvals_TimeDistMat = nan(nTimeBins, nTimeBins, nRepeats);
for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
    for iTime1 = 1:nTimeBins
        clear Vec1
        Vec1 = squeeze(LowerDistVec(iTime1, iRepeat,:));
        for iTime2 = 1:nTimeBins
            clear Vec2
            Vec2 = squeeze(nanmean(LowerDistVec(iTime2,AvgOver, :),2));
            [r,p] = corr(Vec1, Vec2);
            Dist = 1-r;
            LowerTimeDistMat(iTime1, iTime2, iRepeat) = Dist;
            Lower_rvals_TimeDistMat(iTime1, iTime2, iRepeat)=r;
        end
    end
end
LowerAvgTimeDistMat = nanmean(LowerTimeDistMat,3);
LowerAvg_rvals_TimeDistMat = nanmean(Lower_rvals_TimeDistMat,3);

Mean_rvals_TimeDistMat = (Lower_rvals_TimeDistMat+Upper_rvals_TimeDistMat)./2;

AvgTimeDistMat =(LowerAvgTimeDistMat+UpperAvgTimeDistMat)./2;
Avg_rvals_TimeDistMat = (LowerAvg_rvals_TimeDistMat+UpperAvg_rvals_TimeDistMat)./2;

%% shuffled permutations for statistical validation
ShuffledTimeDistMat =  nan(1000,nTimeBins, nTimeBins);
ShuffledAvg_rvals_TimeDistMat=  nan(1000,nTimeBins, nTimeBins);
parfor iPerm = 1:1000
    ShuffledPairwiseDistMat = nan(nStimuli, nStimuli, nTimeBins, nRepeats);
    Shuffled_r_values = nan(nStimuli, nStimuli, nTimeBins, nRepeats);
    Shuffled = randperm(nElec);

    for iTime = 1:nTimeBins
    for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];
     StimMat1 = squeeze(ROI_Data(Shuffled, iTime, :, iRepeat));
    
     StimMat2 =squeeze(nanmean(ROI_Data(:, iTime, :, AvgOver),4));

    [r,p] = corr(StimMat1, StimMat2, 'Type', 'Pearson', 'rows', 'complete');
    Dist = 1-r;
    ShuffledPairwiseDistMat(:,:,iTime,iRepeat) = Dist;
    Shuffled_r_values(:,:,iTime,iRepeat) = r;
    
    end
    end
     AllPerShuffledPairwiseDistMat(iPerm,:,:,:,:) = ShuffledPairwiseDistMat;
     AllPerShuffled_r_values(iPerm,:,:,:,:) = Shuffled_r_values;
   
end
%% Shuffling - Divide RDMs to upper and lower half (later averaging across the results of both halves)
%Upper half
for iPer=1:1000

for iTime = 1:nTimeBins
    for iRepeat = 1:nRepeats
    tmpVec=[];
    c=0;
    for iRow = 1:nStimuli
        c=c+1;
        tmpVec = [tmpVec ; squeeze(AllPerShuffledPairwiseDistMat(iPer,iRow, c:28, iTime,iRepeat))];
    end
    ShuffledUpperDistVec(iPer,iTime, iRepeat,:) = tmpVec;
end
end
%Lower half
for iTime = 1:nTimeBins
    for iRepeat = 1:nRepeats
    tmpVec=[];
    c=0;
    for iRow = 1:nStimuli
        c=c+1;
        tmpVec = [tmpVec  squeeze(AllPerShuffledPairwiseDistMat(iPer,c:28, iRow, iTime, iRepeat))];
    end
    ShuffledLowerDistVec(iPer,iTime,iRepeat,:) = tmpVec';
    end
end
end
%% shuffling- get distance matrices between shuffled RDMs from different time bins
UpperShuffledTimeDistMat = nan(1000,nTimeBins, nTimeBins, nRepeats);
UpperShuffledrvals_TimeDistMat = nan(1000,nTimeBins, nTimeBins, nRepeats);
parfor iPer=1:1000
for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];

        Mat1 = squeeze(ShuffledUpperDistVec(iPer,:, iRepeat,:));
      
            
            Mat2 = squeeze(nanmean(ShuffledUpperDistVec(iPer,:,AvgOver, :),3));
            [r,p] = corr(Mat1', Mat2');
            Dist = 1-r;
            UpperShuffledTimeDistMat(iPer,:, :, iRepeat) = Dist;
            UpperShuffledrvals_TimeDistMat(iPer,:,:, iRepeat)=r;
 
end
end
UpperAvgShuffledTimeDistMat = nanmean(UpperShuffledTimeDistMat,4);
UpperAvgShuffled_rvals_TimeDistMat = nanmean(UpperShuffledrvals_TimeDistMat,4);

LowerShuffledTimeDistMat = nan(1000,nTimeBins, nTimeBins, nRepeats);
LowerShuffledrvals_TimeDistMat = nan(1000,nTimeBins, nTimeBins, nRepeats);
parfor iPer=1:1000
for iRepeat= 1:nRepeats
    AvgOver=[1:4];
    AvgOver(AvgOver==iRepeat)=[];

        Mat1 = squeeze(ShuffledLowerDistVec(iPer,:, iRepeat,:));
      
            
            Mat2 = squeeze(nanmean(ShuffledLowerDistVec(iPer,:,AvgOver, :),3));
            [r,p] = corr(Mat1', Mat2');
            Dist = 1-r;
            LowerShuffledTimeDistMat(iPer,:, :, iRepeat) = Dist;
            LowerShuffledrvals_TimeDistMat(iPer,:,:, iRepeat)=r;
 
end
end
LowerAvgShuffledTimeDistMat = nanmean(LowerShuffledTimeDistMat,4);
LowerAvgShuffled_rvals_TimeDistMat = nanmean(LowerShuffledrvals_TimeDistMat,4);


AvgShuffledTimeDistMat = (LowerAvgShuffledTimeDistMat+UpperAvgShuffledTimeDistMat)./2;
AvgShuffled_rvals_TimeDistMat = (LowerAvgShuffled_rvals_TimeDistMat+UpperAvgShuffled_rvals_TimeDistMat)./2;



%% cluster correction
clusterTh=1.96; % cluster-defining threshold in zscore units 
I = size(AvgShuffledTimeDistMat,1);
meanShuffledData = nanmean(AvgShuffledTimeDistMat(:));

v = nan(I,1);
for i = 1:I
    tmp = squeeze(AvgShuffledTimeDistMat(i,:,:));
    v(i) = nanvar(tmp(:));
end
MSD = sqrt(nanmean(v));   % mean variance computed across shuffled maps
TH = [meanShuffledData-(clusterTh*MSD) meanShuffledData+(clusterTh*MSD)]; 
shuffledDataThr = AvgShuffledTimeDistMat<TH(1)|AvgShuffledTimeDistMat>TH(2);
realDataThr = AvgTimeDistMat<TH(1)|AvgTimeDistMat>TH(2);

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
 %% Plot R-value timecourses for a few time points
 % subplot, each plot shows mean r timecourse
Time = ROI_Struct.Time;
 Selected_Tpoints = round([-0.1, 0.3, 0.7, 1.1,1.5,  2],1);
 
 [I, i1, i2]=intersect(Time,Selected_Tpoints);
 figure;
 set(gcf,'Position',[100 100 1800 300])
 
 for iTC = 1:numel(I)
     R_TC = Avg_rvals_TimeDistMat (i1(iTC),:);
     
     pTC = sig(i1(iTC),:);
     
     R_TC_std =nanstd(squeeze(Mean_rvals_TimeDistMat(i1(iTC),:,:))');
     
     subplot(1,6, iTC);
     hold on;
     h1=shadedErrorBar(Time,R_TC,R_TC_std./sqrt(4),{ 'color', [ 25/255,25/255,112/255],'LineWidth',4,'LineSmoothing','on'},1);
     
     
     hold on;
     SigTimes = Time(find(pTC==1));
     
     plot(SigTimes, -0.08*ones(1,numel(find(pTC==1))), 'Color', [1 174/255, 86/255], 'LineWidth', 4)
     
     title(['Time bin: ' num2str(Selected_Tpoints(iTC)) 's'])
     xlim([-0.5 2.2]);
     ylim([-0.2 0.9])
     xticks([-0.5:0.5:2])
     yticks([-0.2:0.2: 0.8]);
     boxplotX = [-0.75,-0.01,0,1.5,1.51,2.3];
     boxplotY = [-0.195,-0.195,-0.12,-0.12,-0.195,-0.195];
     
     hold on; plot(boxplotX, boxplotY, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
     hold on; plot([Selected_Tpoints(iTC), Selected_Tpoints(iTC)], [-0.2, 1], 'k--');
     xlabel('Time [sec]')
     
     ylabel('Correlation (r)')
     set(gca, 'FontSize', 14);
     box off
 end


end