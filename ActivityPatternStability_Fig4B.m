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
TmpSig = sig+p_maskedAll;
SigMat = TmpSig==2;
%% plot correlation time courses

  Selected_Tpoints = [-0.1, 0.3, 0.7, 1.1,1.5,  2];
 [I, i1, ~]=intersect(Time,Selected_Tpoints);
 figure;
 set(gcf,'Position',[100 100 1800 300])
 
 for iTC = 1:numel(I)
     R_TC =squeeze(nanmean( Avg_r_values(i1(iTC),:,:),3)); %average across stimuli
    
      pTC = SigMat(i1(iTC),:);
      
     R_TC_std =nanstd(squeeze(Avg_r_values(i1(iTC),:,:)));
   
     subplot(1,6, iTC); 
 hold on;
   h1=shadedErrorBar(Time,R_TC,R_TC_std./sqrt(28),{ 'color', [ 25/255,25/255,112/255],'LineWidth',4,'LineSmoothing','on'},1);

       hold on;
      plot(Time(find(pTC==1)), -0.08*ones(1,numel(find(pTC==1))), 'Color', [1 174/255, 86/255], 'LineWidth', 4)
     title(['Time bin: ' num2str(Selected_Tpoints(iTC)) 's'])
 xlim([-0.5 2.2]);
 ylim([-0.2 0.8])
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