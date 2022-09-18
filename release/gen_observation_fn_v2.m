function Z= gen_observation_fn_v2(model,X,W,s,mode)
cam_mat = model.cam_matrix(:,:,s);
%linear observation equation (position components only)
if ~isnumeric(W)
    if strcmp(W,'noise')
        vec_D = diag(model.D);
        WW(1:2,:) = diag(vec_D(1:2))*randn(size(vec_D,1)/2,size(X,2));
        WW(3,:) = lognrnd(0,vec_D(3),[1 size(X,2)]); % width
        WW(4,:) = lognrnd(0,vec_D(4),[1 size(X,2)]); % width
    elseif strcmp(W,'noiseless')
        WW= zeros(size(model.D,1),size(X,2));
    end
else
    WW = W;
end

if isempty(X)
    Z= [];
else
    X([7 8 9],:) = exp((X([7 8 9],:))/model.scale);
    
    for i = 1 :size(X,2)
        
        temp =X([1 3 5 7 8 9],i);
        
        % pnt1 = [temp(1)+(temp(4));temp(2);temp(3)/2]; % - right
        % pnt2 = [temp(1)-(temp(4));temp(2);temp(3)/2]; %  - left
        % pnt3 = [temp(1);temp(2)+(temp(5));temp(3)/2]; % | right
        % pnt4 = [temp(1);temp(2)-(temp(5));temp(3)/2]; % | left
        % pnt5 = [temp(1);temp(2);temp(3)];
        % pnt6 = [temp(1);temp(2);0];
        pnt1 = [temp(1)+(temp(4));temp(2);temp(3)]; % - right
        pnt2 = [temp(1)-(temp(4));temp(2);temp(3)]; %  - left
        pnt3 = [temp(1);temp(2)+(temp(5));temp(3)]; % | right
        pnt4 = [temp(1);temp(2)-(temp(5));temp(3)]; % | left
        pnt5 = [temp(1);temp(2);temp(3)+temp(6)];
        pnt6 = [temp(1);temp(2);temp(3)-temp(6)];
        vet2 = [pnt1 pnt2 pnt3 pnt4 pnt5 pnt6];
        
        
        
        temp = cam_mat*[vet2; ones(1,size(vet2,2))];
        vertices =  temp([1 2],:) ./ temp(3,:)   ;
        % patch('Faces',facs,'Vertices',vertices','FaceColor',[0.1 0.1 0.1],'FaceAlpha',0.2); hold on;
        
        x_2 = max(vertices(1,:));
        x_1 = min(vertices(1,:));
        y_2 = max(vertices(2,:));
        y_1 = min(vertices(2,:));
        %             bbs(:,i) = [x_1 y_1 log(x_2-x_1) log(y_2-y_1)]';
        bbs(:,i) = [x_1 y_1 log(x_2-x_1)+model.meas_n_mu(1,mode) log(y_2-y_1)+model.meas_n_mu(2,mode) ]';
        
        
    end
    
    
    %     Z = bbs + [WW(1:2,:) ; zeros(2,size(X,2))];
    Z = bbs + WW;
    
    
    
end