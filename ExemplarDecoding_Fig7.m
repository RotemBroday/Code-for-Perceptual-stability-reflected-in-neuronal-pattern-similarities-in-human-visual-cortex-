%% Exemplar decoding- train and test on different time bins
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
StimNames = {'Face1', 'Face2', 'Face3', 'Face4', 'Face5', 'Face6', 'Face7', 'Face11', 'Face12', 'Face13', 'Face14', 'Face15'...
    'Face16', 'Face17','Place1', 'Place2', 'Place3', 'Place4', 'Place5', 'Place6', 'Place7', 'Place11', 'Place12', 'Place13' ...
    'Place14', 'Place15', 'Place16', 'Place17'};
for iROI = 1:numel(ElecGroups)
%% Load datasets

clear CurROI ROI_Data TmpData ROI_Struct
CurROI = ElecGroups{iROI};
TmpData= struct2cell(load([DataPath, CurROI, '.mat' ]));
ROI_Struct = TmpData{1};
Time = ROI_Struct.Time;

nRepeats=4;
nElec = numel(ROI_Struct.Elecs);
nTimeBins = numel(ROI_Struct.Time);
nStimuli = numel(ROI_Struct.Stimuli);
ROI_Data= ROI_Struct.Data;



%% PCA dimensionality reduction
MeanRepsMat = nanmean(ROI_Data ,4);
MeanAcrossTime = squeeze(nanmean(MeanRepsMat(:,:,:),2));
 dim_num = floor(intrinsic_dim(MeanAcrossTime','MLE')); 

[mappedAvg, mapping] = compute_mapping(MeanAcrossTime','PCA',dim_num);

    nd = size(mappedAvg,2);
    [~, mapping_allpc] = pca(MeanAcrossTime', size(MeanAcrossTime', 2));
    lambda = mapping_allpc.lambda ./ sum(mapping_allpc.lambda);
    explained_var = sum(lambda(1:nd));
        % PLOT principle components:
    figure('name','All PCs','color','w','position',[0 0 200 150]);
    imagesc(1:28,1:size(mappedAvg,2),mappedAvg'); title(sprintf('PCA - viewing (%.1f%% explained var)',explained_var*100));
 caxis([-3 3]); axis xy; axis square
    ylabel('PC #'); xlabel('Item #');
    set(gca,'ytick',[0:5:20],'xtick',[1 15 28])
    cbar; pos=get(gca,'position'); title('score');
    
          %% map all data to pca space
    PCA_AllData = [];
    for iTime = 1:numel(Time)
        for iStim =1:nStimuli
            for iRepeat = 1:nRepeats
                clear tmp
                tmp = ROI_Data(:, iTime, iStim, iRepeat);
                tmp(isnan(tmp)) = 0; % replace NaN with zeros (i.e. the mean) prior to dimentionality reduction
               PCA_AllData(:, iTime, iStim, iRepeat) =  out_of_sample(tmp',mapping)'; 
            end
        end
    end

    %% decoding
PCAScoreVecExemplarSelecAllStim = [];
PerMeanScoreVec_ExemplarSelecAllStim = [];
for TrainingTimeInd = 1:numel(Time)
    for TestTimeInd = 1:numel(Time)
        parfor iIter = 1:1000
            %chose test data- random trial for each stimuli
            TestData = nan(dim_num, nStimuli);
            TrainData = nan(dim_num,nStimuli);
            TestTrialsPerStim = randi([1 nRepeats],nStimuli,1);
            TrainTrialsPerStim = [];
            for iStim = 1:nStimuli
                tmp = [1:4];
                tmp(TestTrialsPerStim(iStim)==tmp)=[];
                TrainTrialsPerStim(iStim, :)=tmp;
                TestData(:,iStim) = PCA_AllData(:,TestTimeInd,iStim,TestTrialsPerStim(iStim));
                TrainData(:,iStim) = nanmean(PCA_AllData(:,TrainingTimeInd,iStim,tmp),4);
            end
            TestLabels = StimNames;
            TrainLabels = StimNames;
            CorrectCounter = 0;
            
            for iStim = 1:nStimuli
                CorrMat = corr(TrainData, TestData, 'rows',  'complete');
                [M I] = max(CorrMat(:));
                [row,col] = ind2sub(size(CorrMat),I);
                if strcmpi(TrainLabels{row}, TestLabels{col})
                    CorrectCounter = CorrectCounter+1;
                end
                TrainData(:,row) = [];
                TestData(:, col)=[];
                TrainLabels(row)=[];
                TestLabels(col) = [];
            end
            PCAScoreVecExemplarSelecAllStim(iIter,TrainingTimeInd, TestTimeInd) = (CorrectCounter/28)*100;
        end
 disp(TrainingTimeInd); disp(TestTimeInd)       
     save('PCAScoreVecExemplarSelecAllStim', 'PCAScoreVecExemplarSelecAllStim')   
             %% shuffling- 1000 permutations
        
        parfor iPer = 1:1000
            Shuffled = randperm(nStimuli);
            ShuffledScoreVec = [];
            for iIter = 1:1000
                %chose test data- random trial for each stimuli
                TestData = nan(dim_num, nStimuli);
                TrainData = nan(dim_num,nStimuli);
                TestTrialsPerStim = randi([1 nRepeats],nStimuli,1);
                TrainTrialsPerStim = [];
                for iStim = 1:nStimuli
                    tmp = [1:4];
                    tmp(TestTrialsPerStim(iStim)==tmp)=[];
                    TrainTrialsPerStim(iStim, :)=tmp;
                    TestData(:,iStim) = PCA_AllData(:,TestTimeInd,iStim,TestTrialsPerStim(iStim));
                    TrainData(:,iStim) = nanmean(PCA_AllData(:,TrainingTimeInd,iStim,tmp),4);
                end
                TestLabels = StimNames(Shuffled);
                TrainLabels = StimNames;
                CorrectCounter = 0;
                
                for iStim = 1:nStimuli
                    CorrMat = corr(TrainData, TestData, 'rows',  'complete');
                    [M I] = max(CorrMat(:));
                    [row,col] = ind2sub(size(CorrMat),I);
                    if strcmpi(TrainLabels{row}, TestLabels{col})
                        CorrectCounter = CorrectCounter+1;
                    end
                    TrainData(:,row) = [];
                    TestData(:, col)=[];
                    TrainLabels(row)=[];
                    TestLabels(col) = [];
                end
                ShuffledScoreVec(iIter) = (CorrectCounter/28)*100;
            end
            
            PerMeanScoreVec_ExemplarSelecAllStim (iPer,TrainingTimeInd, TestTimeInd) = mean( ShuffledScoreVec);
            
        end
        save('PerMeanScoreVec_ExemplarSelecAllStim', 'PerMeanScoreVec_ExemplarSelecAllStim')
    end
end


% 
% load('PerMeanScoreVec_ExemplarSelecAllStim.mat');
% load('PCAScoreVecExemplarSelecAllStim.mat')
%% calculating stats from shuffling
pMat = nan(numel(Time), numel(Time));
for iTrainingTime = 1:numel(Time)
    for iTestTime = 1:numel(Time)
       pMat(iTrainingTime,iTestTime) = 1-sum(mean(PCAScoreVecExemplarSelecAllStim(:,iTrainingTime,iTestTime))>PerMeanScoreVec_ExemplarSelecAllStim (:,iTrainingTime,iTestTime))/1000;
    end
end
     
[p_fdr, p_masked] = fdr( pMat, 0.001);
if iROI==2
  [p_fdr, p_masked] = fdr( pMat, 0.005);  
end
MeanScore = squeeze(mean(PCAScoreVecExemplarSelecAllStim));
 %% cluster correction
clusterTh=1.96; % cluster-defining threshold in zscore units (for two sided ztest at p<0.05: clusterTh=1.96)
I = size(PerMeanScoreVec_ExemplarSelecAllStim,1);
meanShuffledData = nanmean(PerMeanScoreVec_ExemplarSelecAllStim(:));

v = nan(I,1);
for i = 1:I
    tmp = squeeze(PerMeanScoreVec_ExemplarSelecAllStim(i,:,:));
    v(i) = nanvar(tmp(:));
end
MSD = sqrt(nanmean(v));   % mean variance computed across shuffled maps
TH = [meanShuffledData-(clusterTh*MSD) meanShuffledData+(clusterTh*MSD)]; 
shuffledDataThr = PerMeanScoreVec_ExemplarSelecAllStim<TH(1)|PerMeanScoreVec_ExemplarSelecAllStim>TH(2);
realDataThr = MeanScore<TH(1)|MeanScore>TH(2);

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

%% Plot in matrix format time*time bins
figure;h= imagesc(squeeze(mean(PCAScoreVecExemplarSelecAllStim)));

axis square xy;
set(gca, 'XTick', 3:3:28);
set(gca, 'XTickLabels', Time(3:3:28));
set(gca, 'YTick', 3:3:28);
set(gca, 'YTickLabels', Time(3:3:28));
TmpMask = sig+p_masked;
MyMask = TmpMask==2;

MyMask = bwlabel(MyMask,8);
MyMask = MyMask==1;
MyMask=MyMask+0.5;
set(h, 'AlphaData',MyMask);

title({'Exemplar decoding across time', CurROI}) 
xlabel('Training time (s)'); ylabel('Test time (s)')
colormap(gca,'parula');

h2=colorbar() ;
if iROI ==1
    caxis([2 19])
elseif iROI==2
    caxis([2 13])
elseif iROI==3
    caxis([2 17])

end
ylabel(h2, ' Accuracy (% correct)')

end