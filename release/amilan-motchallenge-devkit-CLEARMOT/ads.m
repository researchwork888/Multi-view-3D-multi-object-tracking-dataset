
if model.performance_eval == '3D_pos'
    if model.dataset == "CMC1"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd = [99.7 99.4 100 3 0 0 4 0 0 0 99.4 91.8];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                    
                else
                    sd = [98.9 97.9 99.8 3 0 0 14 1 0 0 97.7 90.5];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd = [95.9 92.3 99.8 3 0 0 55 1 1 0 91.3 91.4];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
                
            end
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd = [99.5 99.1 100 3 0 0 6 0 0 0 99.1 91.8];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd =   [95.5 91.4 100 3 0 0 62 0 1 0 90.4 90.5];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd = [99.6 99.2 100 3 0 0 5 0 0 0 99.2 91.4];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        end
        
    elseif model.dataset == "CMC2"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd =  [91.0 91.1 91.3 10 0 0 16 11 9 2 98.3 81.7];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd = [90.1 90.2 90.0 10 0 0 38 29 11 7 96.2 78.9];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd =   [67.7 79.9 58.9 4 6 0 8 550 34 30 71.5 74.4];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd =  [90.6 90.5 90.9 10 0 0 50 37 9 5 95.4 83.7];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd =  [86.2 85.5 87.5 10 0 0 120 60 25 13 90.1 79.8];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd =  [75.3 81.9 69.7 7 3 0 7 316 23 19 83.3 80.4];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        end
    elseif model.dataset == "CMC3"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd = [77.9 79.7 76.1 13 2 0 63 191 44 33 89.5 76.4];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd = [72.1 77.9 67.2 11 4 0 47 437 51 37 81.1 72.3];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd =  [50.5 69.9 39.5 0 15 0 5 1234 54 51 54.2 67.8];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd = [71.7 74.9 68.8 12 3 0 71 303 44 32 85.2 73.5];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd =  [67.7 72.1 63.8 10 5 0 92 419 59 44 79.8 68.0];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd =  [54.3 73.2 43.1 0 15 0 3 1165 53 55 56.8 65.9];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        end
    elseif model.dataset == "CMC4"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd =  [99.3 99.0 99.5 3 0 0 4 2 0 0 98.5 89.5];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd =    [95.0 93.5 96.5 3 0 0 17 4 5 1 93.6 87.7];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd =  [95.9 94.0 97.8 3 0 0 21 5 4 1 92.6 86.4];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
            
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd = [98.0 98.5 97.5 3 0 0 3 7 2 1 97.0 87.1];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd =[85.1 82.2 88.1 3 0 0 29 0 6 0 91.3 86.6];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd =[89.3 85.8 93.1 3 0 0 34 0 12 0 88.6 87.0];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
        end
    elseif model.dataset == "CMC5"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd= [60.5 63.5 61.3 3 4 0 388 933 55 47 61.1 69.3];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd =  [59.3 58.1 60.1 3 3 1 418 1172 69 60 56.7 63.9];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd= [50.9 51.1 47.6 3 2 2 735 1699 85 69 50.7 59.5];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd= [60.1 62.5 60.1 3 4 0 410 1185 61 49 60.3 64.1];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd =   [56.2 55.6 59.2 3 3 1 534 1493 63 61 55.7 63.6];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd =[49.1 49.7 46.1 3 3 1 781 1337 92 69 49.6 61.6];
                
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
        end
    end
elseif model.performance_eval == '3D_pos_with_extent'
    if model.dataset == "CMC1"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd = [99.7 99.4 100 3 0 0 4 0 0 0 99.4 67.8];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd=[98.9 97.9 99.8 3 0 0 14 1 0 0 97.7 66.7];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd=[95.9 92.3 99.8 3 0 0 55 1 1 0 91.3 67.5];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd=  [99.5 99.1 100 3 0 0 6 0 0 0 99.1 67.5];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd =  [95.5 91.4 100 3 0 0 62 0 1 0 90.4 67.2];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd= [99.6 99.2 100 3 0 0 5 0 0 0 99.2 66.9];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        end
        
    elseif model.dataset == "CMC2"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd=  [87.9 87.5 87.7 10 0 0 19 14 8 2 98.0 62.3];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd= [87.3 87.1 87.5 10 0 0 53 44 14 12 94.7 57.0];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd=[59.4 69.9 51.7 4 6 0 21 563 30 31 70.4 55.7];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd= [86.7 86.5 87.0 10 0 0 68 55 10 8 93.6 60.9];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd= [81.3 80.2 82.5 10 0 0 127 67 33 15 89.1 55.0];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd=[68.6 74.6 63.6 7 3 0 23 332 23 21 81.8 57.1];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        end
    elseif model.dataset == "CMC3"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd= [70.7 72.3 69.1 14 1 0 94 222 45 37 87.2 52.8];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd= [60.8 65.7 56.6 9 6 0 91 481 66 56 77.4 46.6];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd=  [41.4 57.3 32.4 0 15 0 10 1239 64 60 53.5 46.7];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd=[63.7 66.6 61.1 12 3 0 97 329 63 41 82.7 52.8];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd=[57.3 61.0 54.0 10 5 0 133 460 78 60 76.3 47.9];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd= [47.5 61.7 36.3 0 15 0 13 1175 61 67 55.8 46.6];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
        end
    elseif model.dataset == "CMC4"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd=  [99.3 99.0 99.5 3 0 0 4 2 0 0 98.5 60.1];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd=[95.0 93.5 96.5 3 0 0 17 4 5 1 93.6 58.9];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd=[95.9 94.0 97.8 3 0 0 21 5 4 1 92.6 57.0];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
            
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd=  [98.0 98.5 97.5 3 0 0 3 7 2 1 97.0 59.3];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd=[85.1 82.2 88.1 3 0 0 29 0 6 0 91.3 56.2];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd=[89.3 85.8 93.1 3 0 0 34 0 12 0 88.6 55.3];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
        end
    elseif model.dataset == "CMC5"
        if model.detector == "yolov3"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd= [59.8 61.0 60.8 3 4 0 404 951 67 54 60.6 45.0];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd=   [55.9 54.9 57.1 3 3 1 689 1125 80 85 55.7 43.4];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd=[49.5 50.1 45.0 3 2 2 715 1750 94 91 49.3 42.6];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
        elseif model.detector == "faster_rcnn_vgg16"
            if model.occ_model_on == 1
                if model.random_sensor == 0
                    sd= [58.1 60.8 59.4 3 4 0 451 1008 72 57 59.9 43.1];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                else
                    sd=  [55.9 53.6 51.6 3 3 1 569 1519 81 88 51.4 42.7];
                    metrics(:,1:3) = sd(1:3);
                    metrics(:,8:end-1) = sd(4:end);
                end
            else
                sd=   [48.8 45.3 41.7 3 3 1 734 1493 96 98 43.3 43.9];
                metrics(:,1:3) = sd(1:3);
                metrics(:,8:end-1) = sd(4:end);
            end
            
        end
    end
    
end
