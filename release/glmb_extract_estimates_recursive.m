function [est,idxcmp]=glmb_extract_estimates_recursive(glmb,model,meas,est,filter)
%extract estimates via recursive estimator, where
%trajectories are extracted via association history, and
%track continuity is guaranteed with a non-trivial estimator

%extract MAP cardinality and corresponding highest weighted component
[~,mode] = max(glmb.cdn);
M = mode-1;
T= cell(M,1);
J= zeros(2,M);
% mode = NaN(1,M);

tmp = -inf*ones(1,length(glmb.n));
tmp(glmb.n==M) = 0;
[~,idxcmp] = max(glmb.w+tmp);
% [~,idxcmp]= max(exp(glmb.w).*(glmb.n==M));
for m=1:M
    idxptr= glmb.I{idxcmp}(m);
    T{m,1}= glmb.tt{idxptr}.ah;
    J(:,m)= glmb.tt{idxptr}.l;
    %     [~,idxmode] = max(glmb.tt{glmb.I{idxcmp}(m)}.mode);
    %     mode(:,m) = idxmode;
end

H= cell(M,1);
for m=1:M
    H{m}= [num2str(J(1,m)),'.',num2str(J(2,m))];
end

%compute dead & updated & new tracks
[~,io,is]= intersect(est.H,H);
[~,id,in]= setxor(est.H,H);

est.M= M;
est.T= cat(1,est.T(id),T(is),T(in));
est.J= cat(2,est.J(:,id),J(:,is),J(:,in));
est.H= cat(1,est.H(id),H(is),H(in));

%write out estimates in standard format
est.N= zeros(meas.K,1);
est.X= cell(meas.K,1);
est.L= cell(meas.K,1);
est.mode = cell(meas.K,1);
for t=1:length(est.T)
    ks= est.J(1,t);
    bidx= est.J(2,t);
    tah= est.T{t};
    
    w= model.w_birth{bidx}; %log domain
    m= model.m_birth{bidx};
    P= model.P_birth{bidx};
    mode_w = model.mode(bidx,:);%log domain
    for u=1:size(tah,2)
        
        [log_modetemp_predict,log_wtemp_predict,m_predict,P_predict] = ...
            kalman_predict_multiple_IC_MSGLMB(model,mode_w,w,m,P);
        w= log_wtemp_predict;                                                                       %weights of Gaussians for surviving track
        m = m_predict;
        P = P_predict;
        mode_w = log_modetemp_predict; % make sure is in the log domain
        
        %         [m,P] = kalman_predict_multiple(model,m,P);
        %         [m,P] =ekf_predict_multiple(model.dynamics,m,P);
        
        k= ks+u-1;
        
        %         emm= tah(u);
        mc= tah(:,u);
        if strcmp(model.dataset,'CMC1') || strcmp(model.dataset,'CMC2') || strcmp(model.dataset,'CMC3')
            for curren_mode = 1 :  size(mode_w,2)
                log_qz_temp= zeros(model.no_sensors,size(w{:,curren_mode},2));
                %             log_qz_temp= zeros(model.no_sensors,size(w,2));
                for q = 1 : size(mc,1)
                    Z{q} = eval(['meas.bbs',num2str(q),'{',num2str(k),'}']);
                    model.cam_matrix(:,:,q) = eval(['model.cam',num2str(q),'_cam_mat']);
                    if mc(q) == 0 || isnan( mc(q))
                        log_qz_temp(q) = 0;
                        m{:,curren_mode} = m{:,curren_mode};
                        P{:,curren_mode} = P{:,curren_mode};
                        %                     m = m;
                        %                     P = P;
                    else
                        [~,log_qz_temp(q,:),m{:,curren_mode},P{:,curren_mode}] =...
                            ukf_update_per_sensor(Z{q}(:,mc(q)),model,curren_mode,mode_w(:,curren_mode),...
                            m{:,curren_mode},P{:,curren_mode},filter.ukf_alpha,...
                            filter.ukf_kappa,filter.ukf_beta,q); %ukf update for this track and this measuremnent
                        %                     [log_q_mode_temp(q),log_qz_temp(q,:),m,P] =...
                        %                         ukf_update_per_sensor(Z{q}(:,mc(q)),model,curren_mode,mode_w(:,curren_mode),...
                        %                         m,P,filter.ukf_alpha,...
                        %                         filter.ukf_kappa,filter.ukf_beta,q); %ukf update for this track and this measuremnent
                    end
                    
                end
                log_w_temp = sum(log_qz_temp) + w{:,curren_mode} ;                                                                              %unnormalized updated weights
                w{:,curren_mode} = log_w_temp;     % unnormalized
                m{:,curren_mode} = m{:,curren_mode};
                P{:,curren_mode} = P{:,curren_mode};
                %             log_w_temp = sum(log_qz_temp) + w ;                                                                              %unnormalized updated weights
                %             w = sum(log_q_mode_temp) + log_w_temp;     % unnormalized
                %             m = m;
                %             P = P;
                for_cost(curren_mode) = logsumexp(log_w_temp,[],2);
            end
            for  curren_mode = 1:   size(mode_w,2) % normalize here
                w{:,curren_mode} =  w{:,curren_mode} - logsumexp(for_cost,[],2);
                mode_w(:,curren_mode) = log(1);
                w{:,curren_mode} = w{:,curren_mode} - logsumexp(mode_w(:,curren_mode),[],2);
                %             w =  w - logsumexp(for_cost,[],2);
                %             mode_w(:,curren_mode) = logsumexp(w,[],2);
                %             w = w - logsumexp(mode_w(:,curren_mode),[],2);
            end
                   
        else
            for curren_mode = 1 :  size(mode_w,2)
                log_qz_temp= zeros(model.no_sensors,size(w{:,curren_mode},2));
                %             log_qz_temp= zeros(model.no_sensors,size(w,2));
                log_q_mode_temp= zeros(model.no_sensors,size(mode_w(:,curren_mode),2));
                for q = 1 : size(mc,1)
                    Z{q} = eval(['meas.bbs',num2str(q),'{',num2str(k),'}']);
                    model.cam_matrix(:,:,q) = eval(['model.cam',num2str(q),'_cam_mat']);
                    if mc(q) == 0 || isnan( mc(q))
                        log_qz_temp(q) = 0;
                        log_q_mode_temp(q) = 0;
                        m{:,curren_mode} = m{:,curren_mode};
                        P{:,curren_mode} = P{:,curren_mode};
                        %                     m = m;
                        %                     P = P;
                    else
                        [log_q_mode_temp(q),log_qz_temp(q,:),m{:,curren_mode},P{:,curren_mode}] =...
                            ukf_update_per_sensor(Z{q}(:,mc(q)),model,curren_mode,mode_w(:,curren_mode),...
                            m{:,curren_mode},P{:,curren_mode},filter.ukf_alpha,...
                            filter.ukf_kappa,filter.ukf_beta,q); %ukf update for this track and this measuremnent
                        %                     [log_q_mode_temp(q),log_qz_temp(q,:),m,P] =...
                        %                         ukf_update_per_sensor(Z{q}(:,mc(q)),model,curren_mode,mode_w(:,curren_mode),...
                        %                         m,P,filter.ukf_alpha,...
                        %                         filter.ukf_kappa,filter.ukf_beta,q); %ukf update for this track and this measuremnent
                    end
                    
                end
                log_w_temp = sum(log_qz_temp) + w{:,curren_mode} ;                                                                              %unnormalized updated weights
                w{:,curren_mode} = sum(log_q_mode_temp) + log_w_temp;     % unnormalized
                m{:,curren_mode} = m{:,curren_mode};
                P{:,curren_mode} = P{:,curren_mode};
                %             log_w_temp = sum(log_qz_temp) + w ;                                                                              %unnormalized updated weights
                %             w = sum(log_q_mode_temp) + log_w_temp;     % unnormalized
                %             m = m;
                %             P = P;
                for_cost(curren_mode) = sum(log_q_mode_temp) + logsumexp(log_w_temp,[],2);
            end
            for  curren_mode = 1:   size(mode_w,2) % normalize here
                w{:,curren_mode} =  w{:,curren_mode} - logsumexp(for_cost,[],2);
                mode_w(:,curren_mode) = logsumexp(w{:,curren_mode},[],2);
                w{:,curren_mode} = w{:,curren_mode} - logsumexp(mode_w(:,curren_mode),[],2);
                %             w =  w - logsumexp(for_cost,[],2);
                %             mode_w(:,curren_mode) = logsumexp(w,[],2);
                %             w = w - logsumexp(mode_w(:,curren_mode),[],2);
            end
            
        end
        
        [~,idxmode] = max(mode_w);
        [~,idxtrk]= max(w{:,idxmode});
        %         [~,idxtrk]= max(w);
        est.N(k)= est.N(k)+1;
        est.X{k}= cat(2,est.X{k},m{:,idxmode}(:,idxtrk));
        %         est.X{k}= cat(2,est.X{k},m(:,idxtrk));
        est.L{k}= cat(2,est.L{k},est.J(:,t));
        est.mode{k} = cat(2,est.mode{k},idxmode);
        
        if ~(strcmp(model.dataset,'CMC1') || strcmp(model.dataset,'CMC2') || strcmp(model.dataset,'CMC3'))
            for n_mode = 1 : size(mode_w,2)
                [w{:,n_mode},m{:,n_mode},P{:,n_mode}]= glmb_gaus_mode_prune(w{:,n_mode},m{:,n_mode},P{:,n_mode},filter.gauss_elim_threshold);
                [w{:,n_mode},m{:,n_mode},P{:,n_mode}]= glmb_gaus_mode_merge(w{:,n_mode},m{:,n_mode},P{:,n_mode},filter.gauss_merge_threshold);
                [w{:,n_mode},m{:,n_mode},P{:,n_mode}]= glmb_gaus_cap(w{:,n_mode},m{:,n_mode},P{:,n_mode},filter.gauss_max_num_per_track);
            end
        end
    end
end
end
