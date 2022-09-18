function pD_test = Detection_aka_occlusion_Model(extract,model,q,pidx,Ini_model_check,filter)

sensor_pos  = model.sensor_pos{q}; 
cam_mat = model.cam_matrix(:,:,q);
num_of_objects = length(extract); % no. of objects in this hypothesis
if num_of_objects == 0
    pD_test = [];
    return
end
for a = 1 : num_of_objects % Evaluate MAP of GM
    if length(extract{a}.mode_pd) == 1
        [~,indx]=max(extract{a}.mode_pd);
        [~,indxx]=max(extract{a}.w_pd{indx});
        X(:,a) =  extract{a}.m_pd{:,indxx};
        
        P(:,:,a) = extract{a}.P_pd(:,:,indxx); % not used
    else
        [~,indx]=max(extract{a}.mode_pd);
        [~,indxx]=max(extract{a}.w_pd{:,indx});
        X(:,a) =  extract{a}.m_pd{:,indx}(:,indxx);
        
        P(:,:,a) = extract{a}.P_pd{:,indx}(:,:,indxx); % not used
    end
end
% X([7 8],:) =  exp((X([7,8],:))/model.scale); 
X([7 8 9],:) =  exp((X([7,8 9],:))); 


% Ini_model_check = round((filter.H_upd)*sqrt(glmb.w(pidx))/sum(sqrt(glmb.w))) > 0 ;

if Ini_model_check > 0   && model.occ_model_on %round(filter.pd_model*sqrt(glmb.w(pidx))/sum(sqrt(glmb.w))) > 0 
   


    P_3 = [X(1,:) ; X(3,:); X(5,:)];
    for i = 1: size(P_3,2)
        dist(:,i) = norm(P_3(:,i)-sensor_pos, 2);
    end
    [dist_sorted, indx] = sort(dist); 
    check_pool = ones(size(indx))'; % give every object a flag
    
    for  j = 1 : size(indx,2) % start looping through objects in the hypothesis
        if check_pool(j,:) % check which object hasnt been evaluated
            temp = cam_mat*[P_3(:,indx(j)); 1];
            test_point_Img_plane =  temp([1 2],:) ./ temp(3,:)   ;
            f1 =  test_point_Img_plane(1) <= 0  ;
            f2 =  test_point_Img_plane(1) >= model.imagesize(1)  ;
            f3 =  test_point_Img_plane(2) <= 0   ;
            f4 =  test_point_Img_plane(2) >= model.imagesize(2)   ;
            if f1 || f2 || f3 || f4
                pD_test(indx(j),:) = model.Q_D; % if object is not in the image, then assign low pd.
            else
                pD_test(indx(j),:) = model.P_D; % object is in the image. high pd.
            end
            check_pool(j,:) = 0; % unchecked the object when pd has been assigned 
            curr_object = X(:,indx(j))   ; 
            curr_object_centroid = curr_object([1 3 5],:);
%             curr_object_r = curr_object(7,:);
            curr_object_rx = curr_object(7,:);
            curr_object_ry = curr_object(8,:);
                        
%             curr_object_h = curr_object(8,:);
            curr_object_h = curr_object(9,:);

            % Ellipsoid can be represented as 
            % v'*A*v + b*x + c = 0
            % the idea is to represent/substitute v as a line emanating from
            % sensor_pos. Equation of line in 3D is v(t) = e + t*d where t>=0 and d is a directional unit vector
            % This yields: 
            % \alpha*k^2 + \beta*k + \gamma
            % where \alpha = d'*A*d 
            % \beta = b*d + 2*d'*A*e
            % \gamma = e'*A*e + b*e + c 
            % (see: https://www.geometrictools.com/Documentation/PerspectiveProjectionEllipsoid.pdf )
            % if the above has two distinct real roots, then the ray
            % between the sensor and object to be evaluated, 
            % intersects the chosen/current object twice. This confirms
            % object being evaluated is occluded by chosen/current
            % object. 
            
            
%             A  =  [1/(curr_object_r)^2 0 0;
%                     0 1/(curr_object_r)^2 0 ;
%                     0 0 1/(curr_object_h)^2]; % only for axis-aligned, we have zeros in non-diagonal elements.
            A  =  [1/(curr_object_rx)^2 0 0;
                    0 1/(curr_object_ry)^2 0 ;
                    0 0 1/(curr_object_h)^2];
                
                
                % calculating b and c
            [right_eig, eig_vals] = eig(A,'nobalance','matrix');
            eig_vals_vec = diag(eig_vals);
            temp_ellip_c = right_eig'*curr_object_centroid ;
            ggs = [-2*temp_ellip_c.*eig_vals_vec]';
            b_desired = ggs*right_eig';% this is b
            J_const = sum([ggs.^2./(4*eig_vals_vec)']);  
             
            temp_indx = indx.*check_pool'; % indexes for others except for the one chosen above
            for i = 1 : size(temp_indx,2)
                if temp_indx(i) == 0 % skip if other object(s) has been to be evaluated
                    % nothing
                else
                    evaluate =  P_3(:,temp_indx(i)) ;
                    vector_evaluate =  evaluate - sensor_pos;
                    vector_evaluate = vector_evaluate/norm(vector_evaluate);
                    
                    alpha = vector_evaluate'*A*vector_evaluate;
                    beta = b_desired*vector_evaluate + 2*vector_evaluate'*A*sensor_pos;
                    gamma = sensor_pos'*A*sensor_pos + b_desired*sensor_pos + -1+(J_const);
                    
                    root = roots([alpha,beta,gamma]);
                    
                    
                    %plot and see ~ not used 
                    
                    if isreal(root) 
                        check_pool(i,:) = 0;
                        pD_test(indx(i),:) = model.Q_D;
                    end
                    
                end
            end


        end
    end
            
    if any(check_pool), error('Not all states are checked - P_d are not assigned entirely!!!'); end
          
else 
    pD_test = model.P_D*ones(size(X,2),1);

end
%   pD_test = 0.9;
end