cd('/Users/jonahong/Google Drive/PHD 2018/Codes/Visual MOT Challenge DevKit/amilan-motchallenge-devkit-aa05c63476c9')
addpath(genpath('/Users/jonahong/Google Drive/PHD 2018/Codes/Visual MOT Challenge DevKit/amilan-motchallenge-devkit-aa05c63476c9'))

clear all

%%
% benchmarkDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/GT_J_5_2_Cam1.txt/';
% evaluateDetection('c9-train.txt', 'res/MOT17Det/DPM/data', benchmarkDir);
% det_method = {"yolov3";"faster_rcnn_vgg16"};
dataset = 'J_FJ_4'
cam_no =1


if strcmp(dataset, 'J_1_1')
    switch cam_no 
        case 1
                resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/Cam_1/Detection/MOT format/DET_J_1_1_Cam1.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_Cam_1.txt'; %gt 
%                   resDir = gtDataDir; 
        case 2
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/Cam_2/Detection/MOT format/DET_J_1_1_Cam2.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_Cam_2.txt'; %gt 


        case 3 
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/Cam_3/Detection/MOT format/DET_J_1_1_Cam3.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_Cam_3.txt'; %gt 
        case 4
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/Cam_4/Detection/MOT format/DET_J_1_1_Cam4.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_Cam_4.txt'; %gt             
    end

    seqName  = ['Curtin423B_',dataset,'_cam',num2str(cam_no)] ;
elseif strcmp(dataset,'J_FJ_4')
    switch cam_no 
        case 1
                resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/Cam_1/Detection/MOT format/DET_J_FJ_4_Cam1.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/GT_J_FJ_4_Cam1.txt'; %gt 
%                   resDir = gtDataDir; 
        case 2
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/Cam_2/Detection/MOT format/DET_J_FJ_4_Cam2.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/GT_J_FJ_4_Cam2.txt'; %gt 


        case 3 
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/Cam_3/Detection/MOT format/DET_J_FJ_4_Cam3.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/GT_J_FJ_4_Cam3.txt'; %gt 
        case 4
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/Cam_4/Detection/MOT format/DET_J_FJ_4_Cam4.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal (Fall and Jump)/J_FJ_4/GT_J_FJ_4_Cam4.txt'; %gt             
    end

    seqName  = ['Curtin423B_',dataset,'_cam',num2str(cam_no)] ;
    
elseif strcmp(dataset,'J_5_2')
    switch cam_no 
        case 1
                resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/Cam_1/Detection/MOT format/DET_J_5_2_Cam1.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/GT_J_5_2_Cam1.txt'; %gt 
%                   resDir = gtDataDir; 
        case 2
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/Cam_2/Detection/MOT format/DET_J_5_2_Cam2.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/GT_J_5_2_Cam2.txt'; %gt 


        case 3 
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/Cam_3/Detection/MOT format/DET_J_5_2_Cam3.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/GT_J_5_2_Cam3.txt'; %gt 
        case 4
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/Cam_4/Detection/MOT format/DET_J_5_2_Cam4.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 5/J_5_2/GT_J_5_2_Cam4.txt'; %gt             
    end

    seqName  = ['Curtin423B_',dataset,'_cam',num2str(cam_no)] ;
elseif strcmp(dataset,'J_6_1')
    switch cam_no 
        case 1
                resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/Cam_1/Detection/MOT format/DET_J_6_1_Cam1.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/GT_J_6_1_Cam1.txt'; %gt 
%                   resDir = gtDataDir; 
        case 2
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/Cam_2/Detection/MOT format/DET_J_6_1_Cam2.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/GT_J_6_1_Cam2.txt'; %gt 


        case 3 
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/Cam_3/Detection/MOT format/DET_J_6_1_Cam3.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/GT_J_6_1_Cam3.txt'; %gt 
        case 4
            resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/Cam_4/Detection/MOT format/DET_J_6_1_Cam4.txt';
            gtDataDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 6/J_6_1/GT_J_6_1_Cam4.txt'; %gt             
    end

    seqName  = ['Curtin423B_',dataset,'_cam',num2str(cam_no)] ;

elseif strcmp(dataset,'WILDTRACKS')
     switch cam_no 
        case 1
                resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/Image_subsets/C1/Detection/MOT format/DET_WILDTRACKS_Cam1.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_Cam_1.txt'; %gt 
%                   resDir = gtDataDir; 
        case 2
                resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/Image_subsets/C2/Detection/MOT format/DET_WILDTRACKS_Cam2.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_Cam_2.txt'; %gt 
%                   resDir = gtDataDir;
        case 3
                resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/Image_subsets/C3/Detection/MOT format/DET_WILDTRACKS_Cam3.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_Cam_3.txt'; %gt 
%                   resDir = gtDataDir;
        case 4
                resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/Image_subsets/C4/Detection/MOT format/DET_WILDTRACKS_Cam4.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_Cam_4.txt'; %gt 
%                   resDir = gtDataDir;
        case 5
                resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/Image_subsets/C5/Detection/MOT format/DET_WILDTRACKS_Cam5.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_Cam_5.txt'; %gt 
%                   resDir = gtDataDir;
        case 6
                resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/Image_subsets/C6/Detection/MOT format/DET_WILDTRACKS_Cam6.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_Cam_6.txt'; %gt 
%                   resDir = gtDataDir;
        case 7
                resDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/Image_subsets/C7/Detection/MOT format/DET_WILDTRACKS_Cam7.txt';
%     resDir = '/Users/jonahong/Documents/MATLAB/Curtin 423B Lab/Journal 1/GT_J_1_1_WORLD_CENTROID.txt'; %gt 
                gtDataDir = '/Users/jonahong/Documents/WILDTRACKS /Wildtrack_dataset/GT_WILDTRACKS_Cam_7.txt'; %gt 
%                   resDir = gtDataDir;
     end
        seqName  = ['WILDTRACKS_',dataset,'_cam',num2str(cam_no)] ;

end

cls = [2,7,8,12]; %% ambiguous classes
minvis = 0.5;
ref=0:.025:1;
ref=0:.1:1;
showEachRef=1;
% Find out the length of each sequence
% and concatenate ground truth
gtInfoSingle=[];
gtAll={};
detInfoSingle=[];
detAll={};
seqCnt=0;
allFrCnt=0;
evalMethod=1;
gtAllMatrix=zeros(0,6);
detAllMatrix=zeros(0,7);

%% 
for s=1:1
    seqCnt=seqCnt+1;
        mcnt=1;

        detRaw=dlmread(resDir);
        gtRaw = dlmread(gtDataDir);

    
    % 
    gtOne= {};
    detOne = {};
    F =gtRaw(end,1);
    for t=1:F
        allFrCnt=allFrCnt+1;
        
        % keep pedestrians only and vis >= minvis
        exgt=find(gtRaw(:,1)==t & gtRaw(:,8)==1 & gtRaw(:,9)>=minvis);
        gtAll{allFrCnt}=[gtRaw(exgt,3:6) zeros(length(exgt),1)];
        gtOne{t}=[gtRaw(exgt,3:6) zeros(length(exgt),1)];
        
        ng = length(exgt);
        oneFrame=[allFrCnt*ones(ng,1), (1:ng)', gtRaw(exgt,3:6)]; % set IDs to 1..ng
        gtAllMatrix=[gtAllMatrix; oneFrame];

        exdet=find(detRaw(:,1)==t);
        bbox=detRaw(exdet,3:7);
        detAll{allFrCnt}=bbox;
        detOne{t}=bbox;
        
        ng = length(exdet);
        oneFrame=[allFrCnt*ones(ng,1), (1:ng)', detRaw(exdet,3:7)]; % set IDs to 1..ng
        detAllMatrix=[detAllMatrix; oneFrame];        
    end
    

    allFgt(seqCnt) = F;    
    gtInfoSingle(seqCnt).gt=gtOne;
    gtInfoSingle(seqCnt).gtMat=gtRaw(find(gtRaw(:,8)==1 & gtRaw(:,9)>=minvis),1:6);
    detInfoSingle(seqCnt).det = detOne;    
    detInfoSingle(seqCnt).detMat = detRaw;
    
    
        gt0=gtInfoSingle(seqCnt).gt;
        dt0=detInfoSingle(seqCnt).det;
        [gt,dt]=bbGt('evalRes',gt0,dt0);
        [rc,pr,scores,refprcn] = bbGt('compRoc',gt,dt,0,ref);
        
        
%         rc
%         pr
%         score
%         refprcn
        AP = mean(refprcn);
%         pause
        
        detResults(mcnt).mets(seqCnt).rc=rc;
        detResults(mcnt).mets(seqCnt).pr=pr;
        detResults(mcnt).mets(seqCnt).ref=refprcn;
        detResults(mcnt).mets(seqCnt).AP=AP;
        detResults(mcnt).mets(seqCnt).name=seqName;

        gtRawPed = gtInfoSingle(seqCnt).gtMat;
        detRawPed = detInfoSingle(seqCnt).detMat;
        [detMets, detMetsInfo, detMetsAddInfo]=CLEAR_MOD_HUN(gtRawPed,detRawPed);
%         printMetrics(detMets);
        detResults(mcnt).mets(seqCnt).detMets = detMets;
        detResults(mcnt).mets(seqCnt).detMetsInfo = detMetsInfo;
        detResults(mcnt).mets(seqCnt).detMetsAddInfo = detMetsAddInfo;

        
        
        refprstr = '';
        for r=1:length(refprcn)
            refprstr=[refprstr,sprintf('%.4f',refprcn(r))];
            if r<length(refprcn), refprstr=[refprstr,',']; end
        end
        
            refprcn = detResults(mcnt).mets(seqCnt).ref;
            AP = detResults(mcnt).mets(seqCnt).AP;
            detMets = detResults(mcnt).mets(seqCnt).detMets;
            
            fprintf('\t... %s\n',seqName);
            fprintf('Recall:    ')
            for r=1:showEachRef:length(ref)
                fprintf('%6.3f',ref(r)); 
            end
            fprintf('\n')
            fprintf('Precision: ')
            for r=1:showEachRef:length(ref)
                fprintf('%6.3f',refprcn(r)); 
            end
            fprintf('\n');
            fprintf('Average Precision: %.4f\n',AP);      
            printMetrics(detMets);

            fprintf('\n');
        
end


%%
