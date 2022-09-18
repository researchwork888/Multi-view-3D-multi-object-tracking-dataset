dbstop if error
addpath(genpath('./'))
if strcmp(model.dataset, 'CMC1')
    
    resDir = './Results/EST_CMC1_WORLD_CENTROID.txt';
    gtDataDir =  '../CMC1/GT_CMC1_WORLD_CENTROID.txt';%'GT_J_1_1_WORLD_CENTROID.txt'; %gt
    
    benchmark  = ['_',dataset] ;
elseif strcmp(model.dataset,'CMC2')
    resDir = './Results/EST_CMC2_WORLD_CENTROID.txt';
    gtDataDir = '../CMC2/GT_CMC2_WORLD_CENTROID.txt'; %gt
    
    benchmark  = ['_',dataset];
    
elseif strcmp(model.dataset,'CMC3')
    resDir = './Results/EST_CMC3_WORLD_CENTROID.txt';
    gtDataDir = '../CMC3/GT_CMC3_WORLD_CENTROID.txt'; %gt
    
    benchmark  = ['_',dataset];
    
elseif strcmp(model.dataset,'CMC4')
    resDir = './Results/EST_CMC4_WORLD_CENTROID.txt';
    gtDataDir = '../CMC4/GT_CMC4_WORLD_CENTROID.txt'; %gt
    benchmark  = ['_',dataset];
    
elseif strcmp(model.dataset,'CMC5')
    resDir = './Results/EST_CMC5_WORLD_CENTROID.txt';
    
    gtDataDir = '../CMC5/GT_CMC5_WORLD_CENTROID.txt'; %gt
    
    benchmark  = ['_',dataset];
    
end

world = 1;
threshold = 0.5;
sequenceName = benchmark;
%%
gtFilename = gtDataDir;

gtdata = dlmread(gtFilename);
gtMat{1} = gtdata;

%%
resFilename = resDir;
resdata = dlmread(resFilename);
resMat{1} = resdata;
%%
frameIdPairs = resMat{1}(:,1:2);
[u,I,~] = unique(frameIdPairs, 'rows', 'first');
hasDuplicates = size(u,1) < size(frameIdPairs,1);
if hasDuplicates
    ixDupRows = setdiff(1:size(frameIdPairs,1), I);
    dupFrameIdExample = frameIdPairs(ixDupRows(1),:);
    rows = find(ismember(frameIdPairs, dupFrameIdExample, 'rows'));
    
    errorMessage = sprintf('Invalid submission: Found duplicate ID/Frame pairs in sequence %s.\nInstance:\n', sequenceName);
    errorMessage = [errorMessage, sprintf('%10.2f', resMat{ind}(rows(1),:)), sprintf('\n')];
    errorMessage = [errorMessage, sprintf('%10.2f', resMat{ind}(rows(2),:)), sprintf('\n')];
    assert(~hasDuplicates, errorMessage);
end
%%
%  [metsCLEAR, mInf, additionalInfo] = CLEAR_MOT_HUN(gtMat{1}, resMat{1}, threshold, world);
[metsCLEAR, mInf, additionalInfo] = CLEAR_MOT_HUN_gIoU(gtMat{1}, resMat{1}, threshold, world);

metsID = IDmeasures(gtMat{1}, resMat{1}, threshold, world);
mets = [metsID.IDF1, metsID.IDP, metsID.IDR, metsCLEAR];
allMets(1).name = sequenceName;
allMets(1).m    = mets;
allMets(1).IDmeasures = metsID;
allMets(1).additionalInfo = additionalInfo;
fprintf('%s\n', sequenceName);
printMetrics(mets);
fprintf('\n');
evalFile = fullfile(resDir, sprintf('eval_%s.txt',sequenceName));
%     dlmwrite(evalFile, mets);


%%