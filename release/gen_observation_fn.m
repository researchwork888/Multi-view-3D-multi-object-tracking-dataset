function Z= gen_observation_fn(model,X,W,s,mode)
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
%     X([7,8],:) = exp((X([7,8],:))/model.scale); 
    X([7 8 9],:) =  exp((X([7,8 9],:))); 

    for i = 1 :size(X,2)
        xc = X(1,i) ;
        yc = X(3,i);
        zc = X(5,i);
        ellipsoid_c = [xc,yc,zc];
%         r = X(7,i) ;% half length radius (xy)
%         hh = X(8,i) ;% half length height (z)
        rx = X(7,i) ;% half length radius (x)
        ry = X(8,i) ;% half length radius (x)
        hh = X(9,i) ;% half length height (z)
        % Quadric general equation
        % 0 = Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J
        % Q = [A D/2 E/2 G/2;
        %      D/2 B F/2 H/2;
        %        E/2 F/2 C I/2;
        %        G/2 H/2 I/2 J];
       
        % calculations for A, B, C ...
%         A = 1/(r^2); B = 1/(r^2); C = 1/(hh^2); % equal x-y , C from hh
        A = 1/(rx^2); B = 1/(ry^2); C = 1/(hh^2); % A from rx, B from ry , C from hh

        % calculations for D, E, F ... 
        % no rotation (axis-aligned) means D, E,F = 0 
        D = 0; E = 0;  F=0;
        
        % calculations for G,H,I,J ... 
        PSD = diag([A,B,C]); 
        [right_eig, eig_vals] = eig(PSD,'nobalance','matrix');
        eig_vals_vec = diag(eig_vals);
        temp_ellip_c = right_eig'*ellipsoid_c' ;
        ggs = [-2*temp_ellip_c.*eig_vals_vec]';
        desired = ggs*right_eig';
        G = desired(1); H = desired(2); I = desired(3);
       
        J_temp = sum([ggs.^2./(4*eig_vals_vec)']) ;% or sum(eig_vals_vec.*(-temp_ellip_c).^2)
        J = -1+(J_temp);
        
        % hence, 
        Q = [A D/2 E/2 G/2;
             D/2 B F/2 H/2;
               E/2 F/2 C I/2;
               G/2 H/2 I/2 J]; % 4x4 matrix
        
        C_t = cam_mat*inv(Q)*cam_mat';
%         C_t = cam_mat*(Q)^(-1)*cam_mat';

        C = inv(C_t);  % 3x3 matrix
        
        % Conic general equation
        % 0 = Ax^2 + By^2 + Cxy + Dx + Ey + F
        
        % C = [A C/2 D/2;
        %      C/2 B E/2;
        %      D/2 E/2 F ];
       
        C_strip = [C(1,1:2);C(2,1:2)]; % 2x2 matrix 
%         clear right_eig eig_vals eig_vals_vec
        [right_eig, eig_vals]= eig(C_strip);
        eig_vals_vec = diag(eig_vals);
        x_and_y_vec = [2.*C(3);2.*C(6)]; %extrack D and E 
        x_and_y_vec_transformed = (x_and_y_vec')*right_eig;
        h_temp = (x_and_y_vec_transformed./eig_vals_vec').*(1/2);
        h_temp_squared = eig_vals_vec'.*(h_temp.^2);
        
        h = -1*h_temp;
        ellipse_c = right_eig*h';
        
        offset = -sum(h_temp_squared) + C(9);

        uu = right_eig(:,1)*sqrt(-offset/eig_vals_vec(1));
        vv = right_eig(:,2)*sqrt(-offset/eig_vals_vec(2));
        if isreal(uu) && isreal(vv)
        
            e = sqrt(uu.*uu + vv.*vv );
            bbox = [ellipse_c - e , ellipse_c + e];
            top_left(1) = min(bbox(1,:));
            top_left(2) = min(bbox(2,:));
            bottm_right(1) = max(bbox(1,:));
            bottm_right(2) = max(bbox(2,:));
            
%             if top_left(1) < 0
%                 top_left(1) = 0;
%             elseif top_left(1) > model.imagesize(1)
%                 top_left(1) = model.imagesize(1);
%             elseif top_left(2) > model.imagesize(2)
%                 top_left(2) = model.imagesize(2);
%             elseif  top_left(2) < 0
%                 top_left(2) =0;
%             end
%             
%             if bottm_right(1) < 0
%                 bottm_right(1) = 0;
%             elseif bottm_right(1) > model.imagesize(1)
%                 bottm_right(1) = model.imagesize(1);
%             elseif bottm_right(2) > model.imagesize(2)
%                 bottm_right(2) = model.imagesize(2);
%             elseif  bottm_right(2) < 0
%                 bottm_right(2) =0;
%             end
            
            bbs_noiseless(:,i) = [top_left(1) , top_left(2), log((bottm_right(1)-top_left(1)))+model.meas_n_mu(1,mode), log((bottm_right(2)-top_left(2)))+model.meas_n_mu(2,mode) ]';
        else 
            top_left = [1 1]; 
            bottm_right = [model.imagesize(1) model.imagesize(2)] ;
            bbs_noiseless(:,i) = [top_left(1) , top_left(2), log((bottm_right(1)-top_left(1)))+model.meas_n_mu(1,mode), log((bottm_right(2)-top_left(2)))+model.meas_n_mu(2,mode) ]';

        end




        %% check projection
%         img = eval(['meas.img',num2str(q),'{',num2str(k),'}']);
%         imshow(img,'border','tight');hold on;
%         fimplicit(@(x,y) C(1).*x.^2 + C(5).*y.^2 + 2.*x.*y.*C(2) + 2.*C(3).*x + 2.*C(6).*y + C(9), [1 1920 ],'linewidth',2)
%         hold on;
%         W = bbs_noiseless(3,end); H = bbs_noiseless(4,end);
%         X = bbs_noiseless(1,end) + floor(W/2);
%         Y = bbs_noiseless(2,end) + floor(H/2);
%         plot(X,Y,'o','linewidth',2); % Center
%         %     X = state(1,i); Y= state(2,i);
%         %     W = state(3,i); H = state(4,i);
%         R1 = [X + floor(W/2);Y + floor(H/2)]; %center
%         R2 = [X - floor(W/2);Y + floor(H/2)]; % bottom left
%         R3 = [X - floor(W/2);Y - floor(H/2)]; % top left
%         R4 = [X + floor(W/2);Y - floor(H/2)]; % top right
%         R5 = [X + floor(W/2);Y + floor(H/2)]; % bottom right
%         %     R1 = [X ; Y];
%         %     R2 = [X + floor(W/2);Y];
%         %     R3 = [X ;Y + floor(H/2)];
%         %     R4 = [X - floor(W/2);Y + floor(H/2)];
%         %     R5 = [X - floor(W/2);Y - floor(H/2)];
%         BOX = [R1 R2 R3 R4 R5];
%         plot(BOX(1,:),BOX(2,:),'linewidth',3,'Color', [1,0,0]);
% 
%%
%         img = eval(['meas.img',num2str(q),'{',num2str(k),'}']);
%         imshow(img,'border','tight');hold on;
% %         fimplicit(@(x,y) C(1).*x.^2 + C(5).*y.^2 + 2.*x.*y.*C(2) + 2.*C(3).*x + 2.*C(6).*y + C(9), [1 1920 ],'linewidth',2)
% %         hold on;
%  for i = 1 :size(X,2)
% 
%         W = bbs(3,i); H = bbs(4,i);
%         XX = bbs(1,i) + floor(W/2);
%         Y = bbs(2,i) + floor(H/2);
%         plot(XX,Y,'o','linewidth',2); % Center
%         %     X = state(1,i); Y= state(2,i);
%         %     W = state(3,i); H = state(4,i);
%         R1 = [XX + floor(W/2);Y + floor(H/2)]; %center
%         R2 = [XX - floor(W/2);Y + floor(H/2)]; % bottom left
%         R3 = [XX - floor(W/2);Y - floor(H/2)]; % top left
%         R4 = [XX + floor(W/2);Y - floor(H/2)]; % top right
%         R5 = [XX + floor(W/2);Y + floor(H/2)]; % bottom right
%         %     R1 = [X ; Y];
%         %     R2 = [X + floor(W/2);Y];
%         %     R3 = [X ;Y + floor(H/2)];
%         %     R4 = [X - floor(W/2);Y + floor(H/2)];
%         %     R5 = [X - floor(W/2);Y - floor(H/2)];
%         BOX = [R1 R2 R3 R4 R5];
%         plot(BOX(1,:),BOX(2,:),'linewidth',3,'Color', [1,0,0]);
% end

 

   
    end
    bbs = bbs_noiseless + WW;
    Z = bbs ;
end