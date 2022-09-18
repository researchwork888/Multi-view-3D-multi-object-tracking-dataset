function meas = gen_meas(model,dataset)

if strcmp(dataset,'CMC1')
    cam1_img_direc = fullfile('../','CMC1','Cam_1','/');
    cam1_detection_direc = fullfile('../','CMC1','Cam_1','Detection','/');
    cam2_img_direc = fullfile('../','CMC1','Cam_2','/');
    cam2_detection_direc = fullfile('../','CMC1','Cam_2','Detection','/');
    cam3_img_direc = fullfile('../','CMC1','Cam_3','/');
    cam3_detection_direc = fullfile('../','CMC1','Cam_3','Detection','/');
    cam4_img_direc = fullfile('../','CMC1','Cam_4','/');
    cam4_detection_direc = fullfile('../','CMC1','Cam_4','Detection','/');
    
    cam1_img_direcs = dir([cam1_img_direc, '*.png']);
    cam2_img_direcs = dir([cam2_img_direc, '*.png']);
    cam3_img_direcs = dir([cam3_img_direc, '*.png']);
    cam4_img_direcs = dir([cam4_img_direc, '*.png']);
    
    cam1_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC1/Cam_1/';
    cam2_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC1/Cam_2/';
    cam3_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC1/Cam_3/';
    cam4_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC1/Cam_4/';
    
    meas.K = size(cam1_img_direcs,1);
    meas.start = 1;
    meas.end = meas.K;
elseif strcmp(dataset,'CMC2')
    cam1_img_direc = fullfile('../','CMC2','Cam_1','/');
    cam1_detection_direc = fullfile('../','CMC2','Cam_1','Detection','/');
    cam2_img_direc = fullfile('../','CMC2','Cam_2','/');
    cam2_detection_direc = fullfile('../','CMC2','Cam_2','Detection','/');
    cam3_img_direc = fullfile('../','CMC2','Cam_3','/');
    cam3_detection_direc = fullfile('../','CMC2','Cam_3','Detection','/');
    cam4_img_direc = fullfile('../','CMC2','Cam_4','/');
    cam4_detection_direc = fullfile('../','CMC2','Cam_4','Detection','/');
    
    cam1_img_direcs = dir([cam1_img_direc, '*.png']);
    cam2_img_direcs = dir([cam2_img_direc, '*.png']);
    cam3_img_direcs = dir([cam3_img_direc, '*.png']);
    cam4_img_direcs = dir([cam4_img_direc, '*.png']);
    
    cam1_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC2/Cam_1/';
    cam2_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC2/Cam_2/';
    cam3_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC2/Cam_3/';
    cam4_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC2/Cam_4/';
    
    meas.K = size(cam1_img_direcs,1);
    meas.start = 1;
    meas.end = meas.K;
elseif strcmp(dataset,'CMC3')
    cam1_img_direc = fullfile('../','CMC3','Cam_1','/');
    cam1_detection_direc = fullfile('../','CMC3','Cam_1','Detection','/');
    cam2_img_direc = fullfile('../','CMC3','Cam_2','/');
    cam2_detection_direc = fullfile('../','CMC3','Cam_2','Detection','/');
    cam3_img_direc = fullfile('../','CMC3','Cam_3','/');
    cam3_detection_direc = fullfile('../','CMC3','Cam_3','Detection','/');
    cam4_img_direc = fullfile('../','CMC3','Cam_4','/');
    cam4_detection_direc = fullfile('../','CMC3','Cam_4','Detection','/');
    
    cam1_img_direcs = dir([cam1_img_direc, '*.png']);
    cam2_img_direcs = dir([cam2_img_direc, '*.png']);
    cam3_img_direcs = dir([cam3_img_direc, '*.png']);
    cam4_img_direcs = dir([cam4_img_direc, '*.png']);
    
    cam1_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC3/Cam_1/';
    cam2_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC3/Cam_2/';
    cam3_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC3/Cam_3/';
    cam4_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC3/Cam_4/';
    
    meas.K = size(cam1_img_direcs,1);
    meas.start = 1;
    meas.end = meas.K;
elseif strcmp(dataset,'CMC4')
    cam1_img_direc = fullfile('../','CMC4','Cam_1','/');
    cam1_detection_direc = fullfile('../','CMC4','Cam_1','Detection','/');
    cam2_img_direc = fullfile('../','CMC4','Cam_2','/');
    cam2_detection_direc = fullfile('../','CMC4','Cam_2','Detection','/');
    cam3_img_direc = fullfile('../','CMC4','Cam_3','/');
    cam3_detection_direc = fullfile('../','CMC4','Cam_3','Detection','/');
    cam4_img_direc = fullfile('../','CMC4','Cam_4','/');
    cam4_detection_direc = fullfile('../','CMC4','Cam_4','Detection','/');
    
    cam1_img_direcs = dir([cam1_img_direc, '*.png']);
    cam2_img_direcs = dir([cam2_img_direc, '*.png']);
    cam3_img_direcs = dir([cam3_img_direc, '*.png']);
    cam4_img_direcs = dir([cam4_img_direc, '*.png']);
    
    cam1_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC4/Cam_1/';
    cam2_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC4/Cam_2/';
    cam3_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC4/Cam_3/';
    cam4_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC4/Cam_4/';
    
    meas.K = size(cam1_img_direcs,1);
    meas.start = 1;
    meas.end = meas.K;
elseif strcmp(dataset,'CMC5')
    cam1_img_direc = fullfile('../','CMC5','Cam_1','/');
    cam1_detection_direc = fullfile('../','CMC5','Cam_1','Detection','/');
    cam2_img_direc = fullfile('../','CMC5','Cam_2','/');
    cam2_detection_direc = fullfile('../','CMC5','Cam_2','Detection','/');
    cam3_img_direc = fullfile('../','CMC5','Cam_3','/');
    cam3_detection_direc = fullfile('../','CMC5','Cam_3','Detection','/');
    cam4_img_direc = fullfile('../','CMC5','Cam_4','/');
    cam4_detection_direc = fullfile('../','CMC5','Cam_4','Detection','/');
    
    cam1_img_direcs = dir([cam1_img_direc, '*.png']);
    cam2_img_direcs = dir([cam2_img_direc, '*.png']);
    cam3_img_direcs = dir([cam3_img_direc, '*.png']);
    cam4_img_direcs = dir([cam4_img_direc, '*.png']);
    

    cam1_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC5/Cam_1/';
    cam2_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC5/Cam_2/';
    cam3_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC5/Cam_3/';
    cam4_frcnn_vgg16 = '../Detections/results/faster_rcnn_resnet50_coco_2018_01_28/CMC5/Cam_4/';
    meas.K = size(cam1_img_direcs,1);
    meas.start = 1;
    meas.end = meas.K;
end


for frame_no = meas.start : meas.end
    if strcmp(model.detector,"yolov3")
        % for camera 1
        try
            i = 1;
            fid1 = fopen([cam1_detection_direc,'detection_',num2str(frame_no-1),'.txt'],'r');
            while true
                tline1 = fgetl(fid1);
                detections1{frame_no}(:,i) = sscanf(tline1,['Bounding Box: Left=', '%d',', Top=','%d',', Right=','%d',', Bottom=','%d']);
                i = i+1;
            end
        catch
            fclose(fid1);
        end
        % for camera 2
        try
            i = 1;
            fid2 = fopen([cam2_detection_direc,'detection_',num2str(frame_no-1),'.txt'],'r');
            while true
                tline2 = fgetl(fid2);
                detections2{frame_no}(:,i) = sscanf(tline2,['Bounding Box: Left=', '%d',', Top=','%d',', Right=','%d',', Bottom=','%d']);
                i = i+1;
            end
        catch
            fclose(fid2);
        end
        % for camera 3
        try
            i = 1;
            fid3 = fopen([cam3_detection_direc,'detection_',num2str(frame_no-1),'.txt'],'r');
            while true
                tline3 = fgetl(fid3);
                detections3{frame_no}(:,i) = sscanf(tline3,['Bounding Box: Left=', '%d',', Top=','%d',', Right=','%d',', Bottom=','%d']);
                i = i+1;
            end
        catch
            fclose(fid3);
        end
        % for camera 4
        try
            i = 1;
            fid4 = fopen([cam4_detection_direc,'detection_',num2str(frame_no-1),'.txt'],'r');
            while true
                tline4 = fgetl(fid4);
                detections4{frame_no}(:,i) = sscanf(tline4,['Bounding Box: Left=', '%d',', Top=','%d',', Right=','%d',', Bottom=','%d']);
                i = i+1;
            end
        catch
            fclose(fid4);
        end
        
    elseif strcmp(model.detector, "faster_rcnn_vgg16")
        % cam 1
        dir_cam1_frcnn_vgg16 = dir([cam1_frcnn_vgg16, '*.csv']);
        cam1_img_file{frame_no} = fullfile(cam1_img_direc, cam1_img_direcs(frame_no).name);
        img_cam1{frame_no} = imread( cam1_img_file{frame_no});
        [img_y,img_x,~] = size( img_cam1{frame_no});
        read_temp = csvread(fullfile(cam1_frcnn_vgg16, dir_cam1_frcnn_vgg16(frame_no).name));
        read_temp = read_temp(round(read_temp(:,2) )== 1,:); % peel out person label
        read_temp = read_temp(read_temp(:,end)>=0.3,:); % only pick iou >= 0.3
        detections1{frame_no}  = zeros(4,size(read_temp,1));
        detections1{frame_no}(1,:) = read_temp(:,4)*img_x;
        detections1{frame_no}(2,:) = read_temp(:,3)*img_y;
        detections1{frame_no}(3,:)  = read_temp(:,6)*img_x;
        detections1{frame_no}(4,:)  = read_temp(:,5)*img_y;
        
        % for camera 2
        dir_cam2_frcnn_vgg16 = dir([cam2_frcnn_vgg16, '*.csv']);
        cam2_img_file{frame_no} = fullfile(cam2_img_direc, cam2_img_direcs(frame_no).name);
        img_cam2{frame_no} = imread( cam2_img_file{frame_no});
        [img_y,img_x,~] = size( img_cam2{frame_no});
        read_temp = csvread(fullfile(cam2_frcnn_vgg16, dir_cam2_frcnn_vgg16(frame_no).name));
        read_temp = read_temp(round(read_temp(:,2) )== 1,:); % peel out person label
        read_temp = read_temp(read_temp(:,end)>=0.3,:); % only pick iou >= 0.3
        detections2{frame_no}  = zeros(4,size(read_temp,1));
        detections2{frame_no}(1,:) = read_temp(:,4)*img_x;
        detections2{frame_no}(2,:) = read_temp(:,3)*img_y;
        detections2{frame_no}(3,:)  = read_temp(:,6)*img_x;
        detections2{frame_no}(4,:)  = read_temp(:,5)*img_y;
        % for camera 3
        dir_cam3_frcnn_vgg16 = dir([cam3_frcnn_vgg16, '*.csv']);
        cam3_img_file{frame_no} = fullfile(cam3_img_direc, cam3_img_direcs(frame_no).name);
        img_cam3{frame_no} = imread( cam3_img_file{frame_no});
        [img_y,img_x,~] = size( img_cam3{frame_no});
        read_temp = csvread(fullfile(cam3_frcnn_vgg16, dir_cam3_frcnn_vgg16(frame_no).name));
        read_temp = read_temp(round(read_temp(:,2) )== 1,:); % peel out person label
        read_temp = read_temp(read_temp(:,end)>=0.3,:); % only pick iou >= 0.3
        detections3{frame_no}  = zeros(4,size(read_temp,1));
        detections3{frame_no}(1,:) = read_temp(:,4)*img_x;
        detections3{frame_no}(2,:) = read_temp(:,3)*img_y;
        detections3{frame_no}(3,:)  = read_temp(:,6)*img_x;
        detections3{frame_no}(4,:)  = read_temp(:,5)*img_y;
        
        % for camera 4
        dir_cam4_frcnn_vgg16 = dir([cam4_frcnn_vgg16, '*.csv']);
        cam4_img_file{frame_no} = fullfile(cam4_img_direc, cam4_img_direcs(frame_no).name);
        img_cam4{frame_no} = imread( cam4_img_file{frame_no});
        [img_y,img_x,~] = size( img_cam4{frame_no});
        read_temp = csvread(fullfile(cam4_frcnn_vgg16, dir_cam4_frcnn_vgg16(frame_no).name));
        read_temp = read_temp(round(read_temp(:,2) )== 1,:); % peel out person label
        read_temp = read_temp(read_temp(:,end)>=0.3,:); % only pick iou >= 0.3
        detections4{frame_no}  = zeros(4,size(read_temp,1));
        detections4{frame_no}(1,:) = read_temp(:,4)*img_x;
        detections4{frame_no}(2,:) = read_temp(:,3)*img_y;
        detections4{frame_no}(3,:)  = read_temp(:,6)*img_x;
        detections4{frame_no}(4,:)  = read_temp(:,5)*img_y;
        
    end
    %
    cam1_img_file{frame_no} = fullfile(cam1_img_direc, cam1_img_direcs(frame_no).name);
    img_cam1{frame_no} = imread( cam1_img_file{frame_no});
    %     bbs1_pre{frame_no}=acfDetect(img_cam1{frame_no},detector)';   uncomment to run pitor detector
    bbs1_pre{frame_no} = convert_bbox_to_bbs(detections1{frame_no}); % using YOLO detector
    feet1{frame_no} = convert_bbs_to_feet(   bbs1_pre{frame_no});
    feet1_gp{frame_no} = homtrans(inv(model.cam1_homo), feet1{frame_no});
    indx1 =  find(feet1_gp{frame_no}(2,:) <= model.YMAX(2) & feet1_gp{frame_no}(2,:) >=  model.YMAX(1) ...
        & feet1_gp{frame_no}(1,:) <= model.XMAX(2) & feet1_gp{frame_no}(1,:) >= model.XMAX(1) ) ;
    bbs1{frame_no} =  bbs1_pre{frame_no}(:,indx1);
    meas.bbs1{frame_no} = bbs1{frame_no} ;
    meas.img1{frame_no} = img_cam1{frame_no};
    
    %   %
    cam2_img_file{frame_no} = fullfile(cam2_img_direc, cam2_img_direcs(frame_no).name);
    img_cam2{frame_no} = imread( cam2_img_file{frame_no});
    %   bbs2_pre{frame_no}=acfDetect(img_cam2{frame_no},detector)';uncomment to run pitor detector
    bbs2_pre{frame_no} = convert_bbox_to_bbs(detections2{frame_no}); % using YOLO detector
    feet2{frame_no} = convert_bbs_to_feet(   bbs2_pre{frame_no});
    feet2_gp{frame_no} = homtrans(inv(model.cam2_homo), feet2{frame_no});
    indx2 =  find(feet2_gp{frame_no}(2,:) <= model.YMAX(2) & feet2_gp{frame_no}(2,:) >=  model.YMAX(1) ...
        & feet2_gp{frame_no}(1,:) <= model.XMAX(2) & feet2_gp{frame_no}(1,:) >= model.XMAX(1) ) ;
    bbs2{frame_no} =  bbs2_pre{frame_no}(:,indx2);
    meas.bbs2{frame_no} =  bbs2{frame_no};
    meas.img2{frame_no} =  img_cam2{frame_no};
    %   %
    cam3_img_file{frame_no} = fullfile(cam3_img_direc, cam3_img_direcs(frame_no).name);
    img_cam3{frame_no} = imread( cam3_img_file{frame_no});
    %   bbs3_pre{frame_no}=acfDetect(img_cam3{frame_no},detector)';   uncomment to run pitor detector
    bbs3_pre{frame_no} = convert_bbox_to_bbs(detections3{frame_no}); % using YOLO detector
    feet3{frame_no} = convert_bbs_to_feet(   bbs3_pre{frame_no});
    feet3_gp{frame_no} = homtrans(inv(model.cam3_homo), feet3{frame_no});
    indx3 =  find(feet3_gp{frame_no}(2,:) <= model.YMAX(2) & feet3_gp{frame_no}(2,:) >=  model.YMAX(1) ...
        & feet3_gp{frame_no}(1,:) <= model.XMAX(2) & feet3_gp{frame_no}(1,:) >= model.XMAX(1)) ;
    bbs3{frame_no} =  bbs3_pre{frame_no}(:,indx3);
    meas.bbs3{frame_no} = bbs3{frame_no};
    meas.img3{frame_no} =  img_cam3{frame_no};
    %
    %   %
    cam4_img_file{frame_no} = fullfile(cam4_img_direc, cam4_img_direcs(frame_no).name);
    img_cam4{frame_no} = imread( cam4_img_file{frame_no});
    %   bbs4_pre{frame_no}=acfDetect(img_cam4{frame_no},detector)';   uncomment to run pitor detector
    bbs4_pre{frame_no} = convert_bbox_to_bbs(detections4{frame_no}); % using YOLO detector
    feet4{frame_no} = convert_bbs_to_feet(   bbs4_pre{frame_no});
    feet4_gp{frame_no} = homtrans(inv(model.cam4_homo), feet4{frame_no});
    indx4 =  find(feet4_gp{frame_no}(2,:) <= model.YMAX(2) & feet4_gp{frame_no}(2,:) >=  model.YMAX(1) ...
        & feet4_gp{frame_no}(1,:) <= model.XMAX(2) & feet4_gp{frame_no}(1,:) >= model.XMAX(1)) ;
    bbs4{frame_no} =  bbs4_pre{frame_no}(:,indx4);
    meas.bbs4{frame_no} = bbs4{frame_no};
    meas.img4{frame_no} =  img_cam4{frame_no};
end



end