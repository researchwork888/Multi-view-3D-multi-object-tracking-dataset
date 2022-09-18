% usage: linear system
function est = estimate_tracjectories_UKF(model,prior_estimator,meas,sensor_index)
    cut_track_threshold = 50  ; 
    alpha = 1 ; 
    kappa = 2 ;
    beta = 2 ;
    est.X= cell(meas.K,1);
    est.N= zeros(meas.K,1);
    est.L= cell(meas.K,1) ; 
    % Cut short term tracjectories
    cut_indc = [] ;
    for trjidx = 1 : length(prior_estimator.labels)
        if size(prior_estimator.thetas{trjidx},2) < cut_track_threshold
            cut_indc = [cut_indc , trjidx] ; 

        end
%         asd = reshape(prior_estimator.thetas{trjidx},[prod(size(prior_estimator.thetas{trjidx})),1]);
%         prior_estimator.thetas{trjidx} = reshape(asd,[model.sensors,size(prior_estimator.thetas{trjidx},1)])';
    end
    prior_estimator.labels(cut_indc) = [] ; 
    prior_estimator.thetas(cut_indc) = [] ;
    prior_estimator.inits(cut_indc) = [] ; 
    % Loop through each trajectory
    for trjidx = 1 : length(prior_estimator.labels)
        start_time = prior_estimator.labels{trjidx}(end-1) ; 
%         end_time = start_time + length(prior_estimator.thetas{trjidx}) - 1 ; 
%         end_time = start_time + size(prior_estimator.thetas{trjidx},1)  - 1 ; 
        end_time = start_time+ size(prior_estimator.thetas{trjidx},2) -1 ;

        forward = cell(end_time-start_time+1,1) ;
        forward{1}.w = prior_estimator.inits{trjidx}.w{1} ; 
        forward{1}.m = prior_estimator.inits{trjidx}.m{1} ; 
        forward{1}.P = prior_estimator.inits{trjidx}.P{1} ;
        for q = 1: length(sensor_index{start_time})
            
            z_idx = prior_estimator.thetas{trjidx}(sensor_index{start_time}(q),1);
             if z_idx ~= 0
            
                 meas.Z_bbs = eval(['meas.bbs',num2str(sensor_index{start_time}(q)),'{',num2str(start_time),'}']);
                 meas.Z{start_time} = meas.Z_bbs;
                 model.cam_matrix(:,:,q) = eval(['model.cam',num2str(sensor_index{start_time}(q)),'_cam_mat']);
                 [forward{1}.w, forward{1}.m , forward{1}.P] = forward_update(meas.Z{start_time}(:,z_idx) , model,forward{1}.w , forward{1}.m,forward{1}.P,alpha,beta,kappa,q,meas) ;
            
            end
        end
        for k = 2 : length(forward)
         
            [forward{k}.w, forward{k}.m , forward{k}.P] = forward_predict(model,forward{k-1}.w , forward{k-1}.m,forward{k-1}.P,alpha,beta,kappa) ; 
           for q = 1: length(sensor_index{start_time + k -1})
                z_idx = prior_estimator.thetas{trjidx}(sensor_index{start_time + k -1}(q),k); 
               if z_idx ~= 0 % if not mis-detected then run measurement update (else keep the same as prediction)
                            
                                meas.Z_bbs = eval(['meas.bbs',num2str(sensor_index{start_time + k -1}(q)),'{',num2str(start_time + k -1),'}']);
                                meas.Z{start_time + k -1} = meas.Z_bbs;
                                 model.cam_matrix(:,:,q) = eval(['model.cam',num2str(sensor_index{start_time + k -1}(q)),'_cam_mat']);
                                [forward{k}.w, forward{k}.m , forward{k}.P] = forward_update(meas.Z{start_time + k -1}(:,z_idx) , model, forward{k}.w , forward{k}.m,forward{k}.P,alpha,beta,kappa,q,meas) ;
                            
               end
            end
        end

        est.N(end_time) = est.N(end_time) + 1 ; 
        est.X{end_time} = [est.X{end_time},forward{end}.m] ; 
        est.L{end_time} = [est.L{end_time} , prior_estimator.labels{trjidx}] ; 
        k_offset = 0 ;
        for k = length(forward)-1 : -1 : 1
            k_offset = k_offset + 1 ; 
%             [forward{k}.w , forward{k}.m , forward{k}.P] = backward_smooth_nonlinear(model,forward{k}.w,forward{k}.m,forward{k}.P, ...
%                 forward{k+1}.m,forward{k+1}.P,alpha,beta,kappa) ; 
                [forward{k}.w , forward{k}.m , forward{k}.P] = backward_smooth_linear(model,forward{k}.w,forward{k}.m,forward{k}.P, ...
            forward{k+1}.m,forward{k+1}.P) ;

            est.N(end_time-k_offset) = est.N(end_time-k_offset) +1 ; 
            est.X{end_time-k_offset} = [est.X{end_time-k_offset} , forward{k}.m] ; 
            est.L{end_time-k_offset} = [est.L{end_time-k_offset} , prior_estimator.labels{trjidx}] ;         
        end
    end
    

end


function [w_predict ,m_predict,P_predict] = forward_predict(model,w,m,P,alpha,kappa,beta)      
    w_predict = w ;
    plength= size(m,2);

    
    m_predict = zeros(size(m));
    P_predict = zeros(size(P));

    for idxp=1:plength
%         [m_temp,P_temp] = ukf_predict_single(model,m(:,idxp),P(:,:,idxp),alpha,kappa,beta);
        [m_temp,P_temp] = kalman_predict_single(model.F,model.Q(:,:,1),m(:,idxp),P(:,:,idxp));

        m_predict(:,idxp) = m_temp;
        P_predict(:,:,idxp) = P_temp;
    end

    function [m_predict,P_predict] = ukf_predict_single(model,m,P,alpha,kappa,beta)

        [X_ukf,u]= ut( [m; zeros(model.v_dim,1) ], blkdiag(P,model.Q), alpha, kappa );
        X_pred= gen_newstate_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.v_dim),:) );

        m_predict = X_pred*u(:);
        X_temp= X_pred- repmat(m_predict,[1 length(u)]);
        u(1)= u(1)+(1-alpha^2+beta);
        P_predict= X_temp*diag(u)*X_temp';
    end
    function [m_predict,P_predict] = kalman_predict_single(F,Q,m,P)

        m_predict = F*m;
        P_predict = Q+F*P*F'; 
    end
end

function [w_update,m_update,P_update] = forward_update(z,model,w,m,P,alpha,kappa,beta,q,meas)
    k = 1;
    w_update = w ;
    plength= size(m,2);
    zlength= size(z,2);
    

    
    m_update = zeros(model.x_dim,plength,zlength);
    P_update = zeros(model.x_dim,model.x_dim,plength);

    for idxp=1:plength
            [m_temp,P_temp] = ukf_update_single(z,model,m(:,idxp),P(:,:,idxp),alpha,kappa,beta,q,meas,k);
            m_update(:,idxp,:) = m_temp;
            P_update(:,:,idxp) = P_temp;
    end

function [m_temp,P_temp] = ukf_update_single(z,model,m,P,alpha,kappa,beta,q,meas,k)
    z([3 4],:)= log(z([3 4],:));
    
    [X_ukf,u]= ut( [m; zeros(model.z_dim(q),1) ], blkdiag(P,model.R(:,:,1)), alpha, kappa );
%     [X_ukf,u]= ut( m, P, alpha, kappa );
    
% Z_pred= gen_observation_fn_v2( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.w_dim),:) ,q,meas,k,z,1);
Z_pred= gen_observation_fn_v2( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.z_dim(q)),:) ,q,1);


% Z_pred= gen_observation_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.w_dim),:) ,q,meas,k,z);


        eta= Z_pred*u(:);

        S_temp= Z_pred- repmat(eta,[1 length(u)]);
        u(1)= u(1)+(1-alpha^2+beta);
        S= S_temp*diag(u)*S_temp';
        Vs= chol(S); inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

        G_temp= X_ukf(1:model.x_dim,:)- repmat(m,[1 length(u)]);
        G= G_temp*diag(u)*S_temp';
        K  = G*iS;

        m_temp = repmat(m,[1 size(z,2)]) + K*(z-repmat(eta,[1 size(z,2)]));
        P_temp = P- G*iS*G';
    end
end

function [w_s , m_s , P_s] = backward_smooth_nonlinear(model , old_w , old_m , old_P , m_s_p , P_s_p,alpha,beta,kappa)
    w_s = old_w; 
    % ukf prediction
    [X_ukf,u]= ut( [old_m; zeros(model.v_dim,1) ], blkdiag(old_P,model.Q),alpha,kappa);
    X_pred= gen_newstate_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.v_dim),:) ); 
    m_predict = X_pred*u(:);
    X_temp_covariance_p= X_pred- repmat(m_predict,[1 length(u)]); % predicted X_temp
    u(1)= u(1)+(1-alpha^2+beta);
    P_n_p= X_temp_covariance_p*diag(u)*X_temp_covariance_p'; % not smoothed predicted covriance
    X_temp_covariance = X_ukf(1:model.x_dim,:) - repmat(old_m,[1 length(u)]) ; % Current time X_temp
    C_n_p= X_temp_covariance*diag(u)*X_temp_covariance_p' ; % cross-covariance between 2 time steps
    % Smoothing stuff
    D = C_n_p /P_n_p ; % smooth gain
    m_s = old_m + D*(m_s_p - m_predict) ; % smoothed mean
    P_s = old_P + D*(P_s_p - P_n_p)*D' ; % smoothed covariance
end

function [w_s , m_s , P_s] = backward_smooth_linear(model , old_w , old_m , old_P , m_s_p , P_s_p)
   
    w_s = old_w ; 
    plength= size(old_m,2);
    
   
    for idxp = 1 : plength
        %Prediction
        m_predict = model.F * old_m(:,idxp) ;
        P_predict = model.F * old_P * model.F' + model.Q ;
        %Linear smoothing
        D = old_P * model.F' / P_predict ;
        P_s = old_P - D*(P_predict - P_s_p)*D' ;
        m_s = old_m + D*(m_s_p - m_predict) ;
    end
end

function X= gen_newstate_fn_CT_smooth(model,Xd,V)

    %nonlinear state space equation (CT model)

    if ~isnumeric(V)
        if strcmp(V,'noise')
            V= model.B*randn(size(model.B,2),size(Xd,2));
        elseif strcmp(V,'noiseless')
            V= zeros(size(model.B,1),size(Xd,2));
        end
    end

    if isempty(Xd)
        X= [];
    else %modify below here for user specified transition model
        X= zeros(size(Xd));
        %-- short hand
        L= size(Xd,2);
        T= model.T;
        omega= Xd(5,:);
        tol= 1e-10;
        %-- pre calcs
        sin_omega_T= sin(omega*T);
        cos_omega_T= cos(omega*T);
        a= T*ones(1,L); b= zeros(1,L);
        idx= find( abs(omega) > tol );
        a(idx)= sin_omega_T(idx)./omega(idx);
        b(idx)= (1-cos_omega_T(idx))./omega(idx);
        %--  x/y pos/vel
        X(1,:)= Xd(1,:)+ a.*Xd(2,:)- b.*Xd(4,:);
        X(2,:)= cos_omega_T.*Xd(2,:)- sin_omega_T.*Xd(4,:);
        X(3,:)= b.*Xd(2,:) + Xd(3,:)+ a.*Xd(4,:);
        X(4,:)= sin_omega_T.*Xd(2,:)+ cos_omega_T.*Xd(4,:);
        %-- turn rate
        X(5,:)= Xd(5,:);
        %-- add scaled noise
        X= X+ model.B2*V;
    end
end