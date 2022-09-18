function truth = gen_truth(model,meas,dataset)
if strcmp(dataset,'CMC1')
    gt_direc = fullfile('../',dataset,'/GT_CMC1_WORLD_CENTROID');
elseif  strcmp(dataset,'CMC2')
    gt_direc = fullfile('../',dataset,'/GT_CMC2_WORLD_CENTROID');
elseif  strcmp(dataset,'CMC3')
    gt_direc = fullfile('../',dataset,'/GT_CMC3_WORLD_CENTROID');
elseif  strcmp(dataset,'CMC4')
    gt_direc = fullfile('../',dataset,'/GT_CMC4_WORLD_CENTROID');
elseif  strcmp(dataset,'CMC5')
    gt_direc = fullfile('../',dataset,'/GT_CMC5_WORLD_CENTROID');
end

fileID = fopen([gt_direc,'.txt'],'r');
A = fscanf(fileID,'%d %d %d %d %d %d %d %f %f %f %f %f %f',[13 inf]);
fclose(fileID);
A = A';
frames =  A(:,1);
truth.X = cell(meas.end,1);
truth.L = cell(meas.end,1);
truth.N = zeros(meas.end,1);

for k = meas.start : meas.end
    label = A(frames == k,2);
    X_temp = A(frames == k,end-6+1:end);
    
    for n = 1 : length(label)
        if ~(label(n) == -1)
            truth.L{k} = cat(2,truth.L{k},[label(n);1]);
            truth.X{k} = cat(2,truth.X{k},[X_temp(n,:)']);
            truth.N(k) = truth.N(k) + 1;
        end
    end
end


end