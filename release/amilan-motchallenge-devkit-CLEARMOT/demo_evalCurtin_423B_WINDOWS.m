clear all
% benchmarkDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/GT_J_5_2_Cam1.txt/';
% evaluateDetection('c9-train.txt', 'res/MOT17Det/DPM/data', benchmarkDir);
addpath(genpath('/Users/jonahong/Google Drive/PHD 2018/Codes/Visual MOT Challenge DevKit/amilan-motchallenge-devkit-aa05c63476c9'))
dataset = 'J_FJ_7'
base_distance = 'giou' %"giou" "norm" 
dte = {'yolov3';'faster_rcnn_vgg16'};
wi_wo_occ_model = { 'with_occ_model'; 'without_occ_model'};
wi_wo_cam_config = {'without_cam_config';'with_cam_config'};
world = 1; % must be 1 
threshold =0.75;%150;1;0.7 depending on unit if using base distance as norm
if 1 % 1 : GLMB-OC; 0 : MCD_GLMB
for de = 1 :size(dte,1)
    for wo = 1 : size(wi_wo_occ_model,1)
        if strcmp(wi_wo_occ_model{wo}, "with_occ_model")
            trig_config = [1 2];
        else
            trig_config = [1];
        end
        for cam_config =trig_config
            if strcmp(dataset, 'J_1_1')
                clear mets
                dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Efficient-MSGLMB (no mode WINDOWS)";
                resDir=fullfile(dir,dataset,dte{de},wi_wo_occ_model{wo},wi_wo_cam_config{cam_config},'EST_J_1_1_WORLD_CENTROID_EXTENT.txt');
%                 resDir = 'EST_J_1_1_WORLD_CENTROID.txt';
                gtDataDir = 'C:\Users\17655280\Documents\TPAMI Dataset\Curtin 423B Lab\Journal 1\GT_J_1_1_WORLD_CENTROID_EXTENT.txt'; %gt
%                 gtDataDir = 'GT_J_1_1_WORLD_CENTROID.txt'
                benchmark  = ['Curtin423B_',dataset] ;
            elseif strcmp(dataset, 'J_5_2')
                clear mets
                dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Efficient-MSGLMB (no mode WINDOWS)";
                resDir=fullfile(dir,dataset,dte{de},wi_wo_occ_model{wo},wi_wo_cam_config{cam_config},"EST_J_5_2_WORLD_CENTROID_EXTENT.txt");
                gtDataDir = 'C:\Users\17655280\Documents\TPAMI Dataset\Curtin 423B Lab\Journal 5\J_5_2\GT_J_5_2_WORLD_CENTROID_EXTENT.txt'; %gt
                benchmark  = ['Curtin423B_',dataset] ;
            elseif strcmp(dataset, 'J_6_1')
                clear mets
                dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Efficient-MSGLMB (no mode WINDOWS)";
                resDir=fullfile(dir,dataset,dte{de},wi_wo_occ_model{wo},wi_wo_cam_config{cam_config},"EST_J_6_1_WORLD_CENTROID_EXTENT.txt");
                gtDataDir = 'C:\Users\17655280\Documents\TPAMI Dataset\Curtin 423B Lab\Journal 6\J_6_1\GT_J_6_1_WORLD_CENTROID_EXTENT.txt'; %gt
                benchmark  = ['Curtin423B_',dataset] ;
            elseif strcmp(dataset, 'J_FJ_4')
                clear mets
                dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Efficient-MSGLMB (mode WINDOWS)";
                resDir=fullfile(dir,dataset,dte{de},wi_wo_occ_model{wo},wi_wo_cam_config{cam_config},"EST_J_FJ_4_WORLD_CENTROID_EXTENT.txt");
                gtDataDir = 'C:\Users\17655280\Documents\TPAMI Dataset\Curtin 423B Lab\Journal (Fall and Jump)\J_FJ_4\GT_J_FJ_4_WORLD_CENTROID_EXTENT.txt'; %gt
                benchmark  = ['Curtin423B_',dataset] ;
            elseif strcmp(dataset, 'J_FJ_5')
                clear mets
                dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Efficient-MSGLMB (mode WINDOWS)";
                resDir=fullfile(dir,dataset,dte{de},wi_wo_occ_model{wo},wi_wo_cam_config{cam_config},"EST_J_FJ_5_WORLD_CENTROID_EXTENT.txt");
                gtDataDir = 'C:\Users\17655280\Documents\TPAMI Dataset\Curtin 423B Lab\Journal (Fall and Jump)\J_FJ_5\GT_J_FJ_5_WORLD_CENTROID_EXTENT.txt'; %gt
                benchmark  = ['Curtin423B_',dataset] ;
            elseif strcmp(dataset, 'J_FJ_7')
                clear mets
                dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Efficient-MSGLMB (mode WINDOWS)";
                resDir=fullfile(dir,dataset,dte{de},wi_wo_occ_model{wo},wi_wo_cam_config{cam_config},"EST_J_FJ_7_WORLD_CENTROID_EXTENT.txt");
                gtDataDir = 'C:\Users\17655280\Documents\TPAMI Dataset\Curtin 423B Lab\Journal (Fall and Jump)\J_FJ_7\GT_J_FJ_7_WORLD_CENTROID_EXTENT.txt'; %gt
                benchmark  = ['Curtin423B_',dataset] ;
            elseif strcmp(dataset, 'WILDTRACKS')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% put int here %%%%%%%%%%%%%%%%%%%%%
                clear mets
%                 dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\WILDTRACKS\efficient_MSGLMB principled update";
                dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\WILDTRACKS\New folder";

                resDir=fullfile(dir,dte{de},wi_wo_occ_model{wo},"EST_WILDTRACKS_WORLD.txt");
                gtDataDir = 'C:\Users\17655280\Documents\TPAMI Dataset\WILDTRACKS\Wildtrack_dataset\GT_WILDTRACKS_WORLD.txt'; %gt
                benchmark  = ['Curtin423B_',dataset] ;
            end            
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
            if strcmp(base_distance,"norm")
                [metsCLEAR, mInf, additionalInfo] = CLEAR_MOT_HUN(gtMat{1}, resMat{1}, threshold, world);
            elseif  strcmp(base_distance,"giou")
                
                [metsCLEAR, mInf, additionalInfo] = CLEAR_MOT_HUN_gIoU(gtMat{1}, resMat{1}, threshold, world);
                
            end
            metsID = IDmeasures(gtMat{1}, resMat{1}, threshold, world);
            mets = [metsID.IDF1, metsID.IDP, metsID.IDR, metsCLEAR];
            allMets(1).name = sequenceName;
            allMets(1).m    = mets;
            allMets(1).IDmeasures = metsID;
            allMets(1).additionalInfo = additionalInfo;
            fprintf('%s\n', [sequenceName '(',dte{de}, ' | '  ,wi_wo_occ_model{wo},' | ',wi_wo_cam_config{cam_config},' | ',base_distance,')']); printMetrics(mets); fprintf('\n');
            evalFile = fullfile(resDir, sprintf('eval_%s.txt',sequenceName));
            %     dlmwrite(evalFile, mets);
            
            
            %%
        end
    end
end

else
            clear mets
%             dir = "C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\WILDTRACKS\efficient_MSGLMB principled update";
            dir = "C:\Users\17655280\Documents\TPAMI Dataset\WILDTRACKS\Wildtrack_dataset";

            resDir=fullfile(dir,"EST_WILDTRACKS_WORLD_MCD_GLMB.txt");
            gtDataDir = 'C:\Users\17655280\Documents\TPAMI Dataset\WILDTRACKS\Wildtrack_dataset\GT_WILDTRACKS_WORLD.txt'; %gt
            benchmark  = ['WILDTRACKS_MCD_GLMB',dataset] ;
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
            if strcmp(base_distance,"norm")
                [metsCLEAR, mInf, additionalInfo] = CLEAR_MOT_HUN(gtMat{1}, resMat{1}, threshold, world);
            elseif  strcmp(base_distance,"giou")
                
                [metsCLEAR, mInf, additionalInfo] = CLEAR_MOT_HUN_gIoU(gtMat{1}, resMat{1}, threshold, world);
                
            end
            metsID = IDmeasures(gtMat{1}, resMat{1}, threshold, world);
            mets = [metsID.IDF1, metsID.IDP, metsID.IDR, metsCLEAR];
            allMets(1).name = sequenceName;
            allMets(1).m    = mets;
            allMets(1).IDmeasures = metsID;
            allMets(1).additionalInfo = additionalInfo;
%             fprintf('%s\n', [sequenceName '(',dte{de}, ' | '  ,wi_wo_occ_model{wo},' | ',wi_wo_cam_config{cam_config},' | ',base_distance,')']); 
            printMetrics(mets); fprintf('\n');
            evalFile = fullfile(resDir, sprintf('eval_%s.txt',sequenceName));
            %     dlmwrite(evalFile, mets);
    
    
end
%%
if strcmp(dataset, 'J_1_1')
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/EST_J_1_1_WORLD_CENTROID.txt';
% %     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
%     gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
    
        resDir = 'EST_J_1_1_WORLD_CENTROID.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
    gtDataDir = 'GT_J_1_1_WORLD_CENTROID.txt'; %gt 

    benchmark  = ['Curtin423B_',dataset] ;
elseif strcmp(dataset,'J_FJ_4')
    % resDir = 'C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Annotation Backup\J_5_2\Results (estimates)\EST_J_5_2_WORLD_CENTROID.txt'; %ESTIMATES   
    resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/EST_J_FJ_4_WORLD_CENTROID.txt';
    % gtDataDir = 'C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Annotation Backup\J_5_2\GT_J_5_2_WORLD_CENTROID.txt'; %gt 
    gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/GT_J_FJ_4_WORLD_CENTROID.txt'; %gt 

    benchmark  = ['Curtin423B_',dataset];
    
elseif strcmp(dataset,'J_5_2')
%     % resDir = 'C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Annotation Backup\J_5_2\Results (estimates)\EST_J_5_2_WORLD_CENTROID.txt'; %ESTIMATES   
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/EST_J_5_2_WORLD_CENTROID.txt';
%     % gtDataDir = 'C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Annotation Backup\J_5_2\GT_J_5_2_WORLD_CENTROID.txt'; %gt 
%     gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/GT_J_5_2_WORLD_CENTROID.txt'; %gt 
        resDir = 'EST_J_5_2_WORLD_CENTROID.txt';
    gtDataDir = 'GT_J_5_2_WORLD_CENTROID.txt'; %gt 

    benchmark  = ['Curtin423B_',dataset];
    
elseif strcmp(dataset,'J_6_1')
%     % resDir = 'C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Annotation Backup\J_5_2\Results (estimates)\EST_J_5_2_WORLD_CENTROID.txt'; %ESTIMATES   
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/EST_J_6_1_WORLD_CENTROID.txt';
% %     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/GT_J_6_1_WORLD_CENTROID.txt';
% 
%     % gtDataDir = 'C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Annotation Backup\J_5_2\GT_J_5_2_WORLD_CENTROID.txt'; %gt 
%     gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/GT_J_6_1_WORLD_CENTROID.txt'; %gt 
        resDir = 'EST_J_6_1_WORLD_CENTROID.txt';
    gtDataDir = 'GT_J_6_1_WORLD_CENTROID.txt'; %gt 
    benchmark  = ['Curtin423B_',dataset];
    
elseif strcmp(dataset,'WILDTRACKS')
    % resDir = 'C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Annotation Backup\J_5_2\Results (estimates)\EST_J_5_2_WORLD_CENTROID.txt'; %ESTIMATES   
    resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/EST_WILDTRACKS_WORLD.txt';
%     resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_WORLD.txt';

    % gtDataDir = 'C:\Users\17655280\Google Drive\PHD 2018\Codes\Visual Tracking (Visual Task)\VGLMB_Jonah\Multiview_ppl_Tracking_GLMB\Curtin Acoustic Lab 423B\Annotation Backup\J_5_2\GT_J_5_2_WORLD_CENTROID.txt'; %gt 
    gtDataDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_WORLD.txt'; %gt 

    benchmark  = ['Curtin423B_',dataset];

end

