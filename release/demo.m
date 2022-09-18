clear all 
close all 
dbstop if error
restoredefaultpath
globalStream = RandStream('mcg16807','Seed',10);
RandStream.setGlobalStream(globalStream);
set(0,'Defaultfigurewindowstyle','normal')
cd('./')
addpath(genpath('./'));
do_print_videos = false; %true - print  or false - no print
%% Settings
% name of dataset 
dataset = 'CMC1'; 
% random camera configuration (1 - activate, 0 - normal 4-camera configuration)
random_sensor =0; 
% occlusion model (1 - activate, 0 - deactivate)
occ_model_on = 1;
% detector (yolov3 or faster_rcnn_vgg16)
detector  = "yolov3" ; %"yolov3" or "faster_rcnn_vgg16" 

% evaluate using 3D positions only or 3D with extent
performance_eval = "3D_pos_with_extent"; % "3D_pos" or "3D_pos_with_extent"
%% Run filter
model = gen_model(dataset,random_sensor,occ_model_on,detector,performance_eval); 
meas=  gen_meas(model,dataset);
truth= gen_truth(model,meas,dataset); % truth.X format : [centroid (xyz),width]

[est,prior_estimator] =   run_filter(model,meas); 
[handle, estimate_handle]= plot_results_3D(model,meas,est,do_print_videos);
figure
plot_gt_OSPA(model,truth,meas,est,estimate_handle,performance_eval);

%% Generate video
if do_print_videos
    video = VideoWriter(['./Results/',dataset,'_3d_Cam_1'], 'MPEG-4'); %#ok<UNRCH>
    video.FrameRate = 4;
    open(video);
    writeVideo(video, handle.img1_click);
    close(video);
    
    video = VideoWriter(['./Results/',dataset,'_3d_Cam_2'], 'MPEG-4');
    video.FrameRate = 4;
    open(video);
    writeVideo(video, handle.img2_click);
    close(video);
    
    video = VideoWriter(['./Results/',dataset,'_3d_Cam_3'], 'MPEG-4');
    video.FrameRate = 4;
    open(video);
    writeVideo(video, handle.img3_click);
    close(video);
    
    video = VideoWriter(['./Results/',dataset,'_3d_Cam_4'], 'MPEG-4');
    video.FrameRate = 4;
    open(video);
    writeVideo(video, handle.img4_click);
    close(video);
    
    video = VideoWriter(['./Results/',dataset,'_3d_plot'], 'MPEG-4');
    video.FrameRate = 4;
    open(video);
    writeVideo(video, handle.threedplot);
    close(video);
    
end















