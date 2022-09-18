function glmb_update = glmb_multisensor_joint_predict_update_approx(glmb,model,filter,meas,k)

num_components = length(glmb.w);
num_sensors = model.no_sensors;
index_sensors = model.index_sensors;
max_num_sensors = num_sensors;
num_meas = zeros(num_sensors,1);
for s = 1:num_sensors
    meas.Z{k,s} = eval(['meas.bbs',num2str(index_sensors(s)),'{',num2str(k),'}']);
    model.cam_matrix(:,:,s) = eval(['model.cam',num2str(index_sensors(s)),'_cam_mat']);
    num_meas(s) = size(meas.Z{k,s},2);
    
end

% ------------------------------------
% STEP 1: Create predicted track table
% ------------------------------------
%   Birth tracks
num_birth = length(model.r_birth);
tt_birth = cell(num_birth,1);
for i = 1:num_birth
    tt_birth{i}.m = model.m_birth{i};
    tt_birth{i}.P = model.P_birth{i};
    tt_birth{i}.w = model.w_birth{i};
    tt_birth{i}.mode = model.mode(i,:); 
    tt_birth{i}.l = [k;i];
    tt_birth{i}.ah = [];
end

% Surviving tracks
num_survive = length(glmb.tt);
tt_survive = cell(num_survive,1);
for i = 1:num_survive
    %     [m_predict,P_predict] = ...
    %       ekf_predict_multiple(model.dynamics,glmb.tt{i}.m,glmb.tt{i}.P);
%         [~,~,m_predict,P_predict] = ...
%             kalman_predict_multiple_IC_MSGLMB(model,1,...
%             glmb.tt{i}.w,glmb.tt{i}.m,glmb.tt{i}.P);
    [log_modetemp_predict,log_wtemp_predict,m_predict,P_predict] = ...
        kalman_predict_multiple_IC_MSGLMB(model,glmb.tt{i}.mode,...
        glmb.tt{i}.w,glmb.tt{i}.m,glmb.tt{i}.P);
    tt_survive{i}.m = m_predict;
    tt_survive{i}.P = P_predict;
    tt_survive{i}.w = log_wtemp_predict; % make sure is in the log domain
    tt_survive{i}.mode = log_modetemp_predict; % make sure is in the log domain
    tt_survive{i}.l = glmb.tt{i}.l;
    tt_survive{i}.ah = glmb.tt{i}.ah;
end

% Union of births and survivals
glmb_predict.tt = cat(1,tt_birth,tt_survive);
num_tracks_predict = num_birth + num_survive;

% --------------------------
% STEP 2: Measurement gating
% --------------------------

if filter.gate_flag
    % Find measurements inside each track's gate for each sensor
    for t = 1:num_tracks_predict
        for s = 1:num_sensors
            glmb_predict.tt{t}.gatemeas{s} = ...
                gate_msmeas_ekf_idx(meas.Z{k,s},...
                filter.gamma,...
                model.sensor,s,...
                glmb_predict.tt{t}.m,...
                glmb_predict.tt{t}.P);
        end
    end
else
    % No gating, include all measurements for all tracks
    for t = 1:num_tracks_predict
        for s = 1:num_sensors
            glmb_predict.tt{t}.gatemeas{s} = 1:size(meas.Z{k,s},2);
        end
    end
end

% --------------------------------------------------
% STEP 3: Compute average survival, death, detection
%         and misdetection probabilities
% -------------------------------------------------

avps = [(model.r_birth) ; zeros(num_survive,1)];
avqs = [log(1-exp(model.r_birth)) ; zeros(num_survive,1)];
for t = 1:num_survive
    temp_ps = compute_pS(model,glmb.tt{t}.mode,glmb.tt{t}.w,glmb.tt{t}.m,glmb.tt{t}.l,k);
    avps(num_birth+t) = (temp_ps);
    avqs(num_birth+t) = log(1-exp(temp_ps));
end

%   avpd = zeros(num_tracks_predict,num_sensors);
%   avqd = zeros(num_tracks_predict,num_sensors);
%   for t = 1:num_tracks_predict
%     for s = 1:num_sensors
%       avpd(t,s) = log(model.sensor.P_D(s));
%       avqd(t,s) = log(1-model.sensor.P_D(s));
%     end
%   end

% -----------------------------
% STEP 4: Create updated tracks
% -----------------------------

allcostc = cell(num_sensors,1);
jointcostc = cell(num_sensors,1);
avpp = zeros(num_tracks_predict,num_sensors);
for s = 1:num_sensors
    allcostc{s} = -inf*ones(length(glmb_predict.tt),1+num_meas(s));
    %     jointcostc{s} = -inf*ones(length(glmb_predict.tt),1+num_meas(s));
    for t = 1:length(glmb_predict.tt)
        for m = glmb_predict.tt{t}.gatemeas{s}
            for curren_mode = 1:  size(glmb_predict.tt{t}.mode,2)
                %       [qz_temp,~,~] = ekf_update_multiple(meas.Z{k,s}(:,m),model.sensor,s,glmb_predict.tt{t}.m,glmb_predict.tt{t}.P);
                %                 [~,qz_temp,~,~] = ukf_update_per_sensor(meas.Z{k,s}(:,m),model,1,1,...
                %                     glmb_predict.tt{t}.m,glmb_predict.tt{t}.P,filter.ukf_alpha,...
                %                     filter.ukf_kappa,filter.ukf_beta,s); %ukf update for this track and this measuremnent
                
                [log_q_mode,log_qz_temp,~,~] = ukf_update_per_sensor(meas.Z{k,s}(:,m),model,curren_mode,glmb_predict.tt{t}.mode(:,curren_mode),...
                    glmb_predict.tt{t}.m{:,curren_mode},glmb_predict.tt{t}.P{:,curren_mode},filter.ukf_alpha,...
                    filter.ukf_kappa,filter.ukf_beta,s); %ukf update for this track and this measuremnent
                w_temp = log_qz_temp + glmb_predict.tt{t}.w{:,curren_mode};
                for_cost(curren_mode) = log_q_mode + logsumexp(w_temp,[],2); 
            end
            
            %             w_temp = qz_temp + glmb_predict.tt{t}.w;
            %             allcostc{s}(t,1+m) = logsumexp(w_temp);
            
            allcostc{s}(t,1+m) = logsumexp(for_cost,[],2); %predictive likelihood
            
        end
    end
    % %     jointcostc{s} = [avps+avqd(:,s) repmat(avps+avpd(:,s),[1 num_meas(s)])] + ...
    % %                     allcostc{s} - log(model.sensor.lambda_c(s)*model.sensor.pdf_c(s));
    %     jointcostc{s} = [ avps+avqd(:,s)  repmat(avps+avpd(:,s),[1 num_meas(s)]) + ...
    %                     allcostc{s}(:,2:end) - log(model.sensor.lambda_c(s)*model.sensor.pdf_c(s))];
    %
    %     %jointcostc{s}(:,1) = avqs;
    %     if (num_meas(s) > 0)
    %       avpp(:,s) = logsumexp(jointcostc{s}(:,2:end),[],2);
    %     else
    %       avpp(:,s) = -inf;
    %     end
end

% Gated measurement index matrix
gatemeasidxs = cell(num_sensors,1);
gatemeasidxc = cell(num_sensors,1);
for s = 1:num_sensors
    gatemeasidxs{s} = zeros(num_tracks_predict,num_meas(s));
    for t = 1:num_tracks_predict
        gatemeasidxs{s}(t,1:length(glmb_predict.tt{t}.gatemeas{s})) = glmb_predict.tt{t}.gatemeas{s};
    end
    gatemeasidxc{s} = (gatemeasidxs{s} > 0);
end

% -------------------------
% STEP 5: Component updates
% -------------------------

%textprogressbar('Processing components: ');

total_comp_request = filter.comp_max_num_request;
lambda_c = model.lambda_c(index_sensors);
pdf_c = model.pdf_c(index_sensors);
weights = glmb.w;
log_sum_sqrt_weights = logsumexp(0.5*weights,[],2);
num_tracks = zeros(1,num_components);
tindices = cell(1,num_components);
mindices = cell(1,num_components);
avps_comp = cell(1,num_components);
avqs_comp = cell(1,num_components);
avpp_comp = cell(1,num_components);
avpd_comp = cell(1,num_components);
avqd_comp = cell(1,num_components);
neglogcostc_tempd = cell(1,num_components);

for pidx = 1:num_components
    num_exists = length(glmb.I{pidx});
    num_tracks(pidx) = num_birth + num_exists;
    tindices{pidx} = [1:num_birth num_birth+glmb.I{pidx}'];
    avps_comp{pidx} = avps(tindices{pidx});
    avqs_comp{pidx} = avqs(tindices{pidx});
    %     avpd_comp{pidx} = avpd(tindices{pidx},:);
    %     avqd_comp{pidx} = avqd(tindices{pidx},:);
    %     avpp_comp{pidx} = avpp(tindices{pidx},:);
    mindices{pidx} = cell(1,num_sensors);
    neglogcostc_tempd{pidx} = cell(1,num_sensors);
    extract =[];
    for a = 1 : size(tindices{pidx},2)
        extract{a}.m_pd = glmb_predict.tt{tindices{pidx}(a)}.m;
        extract{a}.P_pd =  glmb_predict.tt{tindices{pidx}(a)}.P;
        extract{a}.w_pd =  glmb_predict.tt{tindices{pidx}(a)}.w;
        extract{a}.mode_pd = glmb_predict.tt{tindices{pidx}(a)}.mode;
    end
    for s = 1:num_sensors
        lselmask = false(num_tracks_predict,num_meas(s));
        lselmask(tindices{pidx},:) = gatemeasidxc{s}(tindices{pidx},:);
        mindices{pidx}{s} = unique_faster(gatemeasidxs{s}(lselmask));
        mindices{pidx}{s} = [1 1+mindices{pidx}{s}];
        
        num_comp_request = round(exp((log(total_comp_request)+0.5*weights(pidx)-log_sum_sqrt_weights)));
        %       Detection_aka_occlusion_Model
        avpd_comp{pidx}{s}  = Detection_aka_occlusion_Model(extract,model,s,pidx,num_comp_request,filter);
        
        avqd_comp{pidx}{s} = log( 1 - exp(avpd_comp{pidx}{s}) );
        avpd_comp{pidx}{s} = (avpd_comp{pidx}{s});
        jointcostc{s} = -inf*ones(length(glmb_predict.tt),1+num_meas(s));
        
        jointcostc{s}(tindices{pidx},mindices{pidx}{s}) = [ avps_comp{pidx}+avqd_comp{pidx}{s}  repmat(avps_comp{pidx}+avpd_comp{pidx}{s} ,[1 num_meas(s)]) + ...
            allcostc{s}(tindices{pidx},mindices{pidx}{s}(:,2:end)) - log(model.lambda_c(s)*model.pdf_c(s))];
        
        %       if (num_meas(s) > 0)
        avpp(:,s) = logsumexp(jointcostc{s}(:,1:end),[],2);
        %       else
        %           avpp(:,s) = -inf;
        %       end
        neglogcostc_tempd{pidx}{s} = -jointcostc{s}(tindices{pidx},mindices{pidx}{s});
    end
    avpp_comp{pidx} = avpp(tindices{pidx},:);
    
end
glmb_update_w = cell(1,num_components);
glmb_update_I = cell(1,num_components);
glmb_update_n = cell(1,num_components);
tt_update_parent_comp = cell(1,num_components);
tt_update_currah_comp = cell(1,num_components);
tt_update_linidx_comp = cell(1,num_components);
num_comp_gen = zeros(1,num_components);

for pidx = 1:num_components
    
    ndimvec = zeros(1+num_sensors,1);
    ndimvec(1) = length(tindices{pidx});
    for s = 1:num_sensors
        ndimvec(s+1) = size(neglogcostc_tempd{pidx}{s},2);
    end
    
    dcost = avqs_comp{pidx};
    %     scost = logsumexp(avpp_comp{pidx},[],2);
    scost = sum(avpp_comp{pidx},2);
    dprob = dcost - logsumexp(dcost,scost); % ratio  a / a+b ??
    neglogdprob_tempd = -dprob;
    neglogdcost = -dcost;
    
    num_comp_request = round(exp((log(total_comp_request)+0.5*weights(pidx)-log_sum_sqrt_weights)));
    
    seed = pidx;
    if num_comp_request > 0
        %       uasses= gibbs_multisensor_approx(neglogdprob_tempd',...
        %                                        neglogcostc_tempd{pidx},...
        %                                        num_comp_request,...
        %                                        filter.gibbs_burn_in,...
        %                                        filter.gibbs_thinning_factor,...
        %                                        seed);
        
        %         uasses= RANDOMSCAN_gibbs_multisensor_approx_matlab(neglogdprob_tempd',...
        %             neglogdcost',...
        %             neglogcostc_tempd{pidx},...
        %             num_comp_request,...
        %             filter);
        
        uasses= gibbs_multisensor_approx_matlab(neglogdprob_tempd',...
            neglogdcost',...
            neglogcostc_tempd{pidx},...
            num_comp_request,...
            filter);

        %
        
        
        uasses = double(uasses+1);
        uasses(uasses<1) = -inf;
    else
        uasses = zeros(0,num_tracks(pidx),num_sensors);
    end
    
    % Restore original indices of gated measurements,
    % including surivived+missed=0, survived+detected to 1:|Z|
    for s = 1:num_sensors
        tempslice = uasses(:,:,s);
        tempslice(tempslice>0) = mindices{pidx}{s}(tempslice(tempslice>0))-1;
        uasses(:,:,s) = tempslice;
    end
    
    ncg = size(uasses,1);
    glmb_update_w{pidx} = zeros(1,ncg);
    glmb_update_I{pidx} = cell(1,ncg);
    glmb_update_n{pidx} = zeros(1,ncg);
    num_insert = zeros(1,ncg);
    next_index = 1;
    
    % Generate corrresponding jointly predicted/updated hypotheses/components
    update_hypcmp_tmp = cell(1,ncg);
    off_idx = false(num_tracks(pidx),ncg);
    for hidx = 1:ncg
        update_hypcmp_tmp{hidx} = reshape(uasses(hidx,:,:),[num_tracks(pidx) num_sensors]);
        off_idx(:,hidx) = update_hypcmp_tmp{hidx}(:,1)<0;
        num_insert(hidx) = sum(~off_idx(:,hidx));
    end
    
    total_insert = sum(num_insert);
    tt_update_parent_comp{pidx} = zeros(total_insert,1);
    tt_update_currah_comp{pidx} = zeros(total_insert,num_sensors);
    tt_update_linidx_comp{pidx} = zeros(total_insert,1);
    
    for hidx = 1:ncg
        aug_idx = [tindices{pidx}(:) 1+update_hypcmp_tmp{hidx}];
        mis_idx = (update_hypcmp_tmp{hidx}==0);
        det_idx = (update_hypcmp_tmp{hidx}>0);
        %         local_avpdm = reshape([avpd_comp{pidx}{:}],[size(avpd_comp{pidx},1) size(avpd_comp{pidx},2)]);
        %         local_avqdm = reshape([avqd_comp{pidx}{:}],[size(avqd_comp{pidx},1) size(avqd_comp{pidx},2)]);
        extract =[];
        for_occlu_idx = find(~off_idx(:,hidx));
        for a = 1 : size(tindices{pidx}(~off_idx(:,hidx)),2)
            extract{a}.m_pd = glmb_predict.tt{tindices{pidx}(for_occlu_idx(a))}.m;
            extract{a}.P_pd =  glmb_predict.tt{tindices{pidx}(for_occlu_idx(a))}.P;
            extract{a}.w_pd =  glmb_predict.tt{tindices{pidx}(for_occlu_idx(a))}.w;
            extract{a}.mode_pd = glmb_predict.tt{tindices{pidx}(for_occlu_idx(a))}.mode;
        end
        local_avpdm = ([avpd_comp{pidx}{:}]);
        local_avqdm = ([avqd_comp{pidx}{:}]);
        if any(~off_idx(:,hidx))
            for s = 1 : num_sensors
                avpd_comp_temp  = Detection_aka_occlusion_Model(extract,model,s,pidx,num_comp_request,filter);
                avqd_comp_temp = log(1 - exp(avpd_comp_temp));
                local_avpdm(for_occlu_idx,s) = (avpd_comp_temp);
                local_avqdm(for_occlu_idx,s) = (avqd_comp_temp);
            end
        end
        repeated_lambda_c = repmat(lambda_c',[length(tindices{pidx}) 1]);
        repeated_pdf_c = repmat(pdf_c',[length(tindices{pidx}) 1]);
        update_hypcmp_idx = [];
        update_hypcmp_idx(off_idx(:,hidx)) = -inf;
        update_hypcmp_idx(~off_idx(:,hidx)) = ndsub2ind([num_tracks_predict; 1+num_meas],aug_idx(~off_idx(:,hidx),:));
        num_trk = sum(update_hypcmp_idx>0);
        
        glmb_update_w{pidx}(hidx) = ... %-sum(lambda_c) + ...
            num_meas'*log(lambda_c.*pdf_c) + ...
            weights(pidx) + ...
            sum(avps_comp{pidx}(~off_idx(:,hidx))) + ...
            sum(avqs_comp{pidx}(off_idx(:,hidx))) + ...
            sum(sum(local_avpdm(det_idx))) + ...
            sum(sum(local_avqdm(mis_idx))) - ...
            sum(sum(log(repeated_lambda_c(det_idx).*repeated_pdf_c(det_idx))));
        
        if num_trk > 0
            glmb_update_I{pidx}{hidx} = next_index:next_index+num_trk-1;
        end
        glmb_update_n{pidx}(hidx)= num_trk;
        
        if num_insert(hidx) > 0
            tt_update_parent_comp{pidx}(next_index:next_index+num_insert(hidx)-1) = tindices{pidx}(~off_idx(:,hidx))';
            tt_update_currah_comp{pidx}(next_index:next_index+num_insert(hidx)-1,:) = update_hypcmp_tmp{hidx}(~off_idx(:,hidx),:);
            tt_update_linidx_comp{pidx}(next_index:next_index+num_insert(hidx)-1) = update_hypcmp_idx(update_hypcmp_idx>0)';
            next_index = next_index + num_insert(hidx);
        end
        
    end
    
    num_comp_gen(pidx) = ncg;
    
    %percentage_complete = 100*pidx/length(glmb_update.w);
    %textprogressbar(percentage_complete);
    
end

%textprogressbar(' Done');

total_comp_gen = sum(num_comp_gen);
glmb_update.I = cell(1,total_comp_gen);
glmb_update.w = -inf*ones(1,total_comp_gen);
glmb_update.n = zeros(1,total_comp_gen);

comp_idx = 1;
trk_idx = 1;
for pidx = 1:num_components
    
    nc = length(glmb_update_w{pidx});
    for i = 0:nc-1
        glmb_update.I{comp_idx+i} = glmb_update_I{pidx}{i+1} + trk_idx - 1;
    end
    glmb_update.w(comp_idx:comp_idx+nc-1) = glmb_update_w{pidx};
    glmb_update.n(comp_idx:comp_idx+nc-1) = glmb_update_n{pidx};
    comp_idx = comp_idx + nc;
    nt = length(tt_update_parent_comp{pidx});
    tt_update_parent(trk_idx:trk_idx+nt-1) = tt_update_parent_comp{pidx};
    tt_update_currah(trk_idx:trk_idx+nt-1,:) = tt_update_currah_comp{pidx};
    tt_update_linidx(trk_idx:trk_idx+nt-1) = tt_update_linidx_comp{pidx};
    trk_idx = trk_idx + nt;
    
end

% Component updates via posterior weight correction (including generation of track table)
[tt_update_allkey,tt_update_oldidx,tt_update_newidx] = unique(tt_update_linidx);
tt_update_msqz = zeros(length(tt_update_allkey),1);
tt_update = cell(length(tt_update_allkey),1);

tt_predict = glmb_predict.tt(tt_update_parent(tt_update_oldidx));
meas_comb = tt_update_currah(tt_update_oldidx,:);
Z = meas.Z(k,:);

for i = 1:length(tt_update_allkey)
    mc = meas_comb(i,:)';
    %     [qz_temp,m_temp,P_temp] = ekf_msjointupdate_multiple(Z,mc,model.sensor,tt_predict{i}.m,tt_predict{i}.P);
    %     [~,qz_temp,m_temp,P_temp] = ukf_update_joint_sensor(Z,mc,model,1,1,...
    %         tt_predict{i}.m,tt_predict{i}.P,filter.ukf_alpha,...
    %         filter.ukf_kappa,filter.ukf_beta); %ukf update for this track and this measuremnent
    m_temp = tt_predict{i}.m ;
    P_temp  = tt_predict{i}.P ;
    mode_temp = tt_predict{i}.mode;
    
    %     qz_temp(1:model.no_sensors) = 0;
    %     for q = 1 : size(mc,1)
    %         if mc(q) == 0
    %             qz_temp(q) = 0;
    %             m_temp = m_temp;
    %             P_temp = P_temp;
    %         else
    %             [~,qz_temp(q),m_temp,P_temp] = ukf_update_per_sensor(Z{:,q}(:,mc(q)),model,1,1,...
    %                 m_temp,P_temp,filter.ukf_alpha,...
    %                 filter.ukf_kappa,filter.ukf_beta,q); %ukf update for this track and this measuremnent
    %
    %         end
    %
    %     end
    %
    %     w_temp = sum(qz_temp) + tt_predict{i}.w;
    %
    %     tt_update{i}.m = m_temp;
    %     tt_update{i}.P = P_temp;
    %     tt_update{i}.w = w_temp - logsumexp(w_temp);
    %     tt_update_msqz(i) = logsumexp(w_temp);
    for curren_mode = 1 :  size(mode_temp,2)
        log_qz_temp= zeros(model.no_sensors,size(tt_predict{i}.w{:,curren_mode},2));
        log_q_mode_temp= zeros(model.no_sensors,size(tt_predict{i}.mode(:,curren_mode),2));
        for q = 1 : size(mc,1)
            if mc(q) == 0
                log_qz_temp(q) = 0;
                log_q_mode_temp(q) = 0;
                m_temp{:,curren_mode} = m_temp{:,curren_mode};
                P_temp{:,curren_mode} = P_temp{:,curren_mode};
            else
                [log_q_mode_temp(q),log_qz_temp(q,:),m_temp{:,curren_mode},P_temp{:,curren_mode}] =...
                    ukf_update_per_sensor(Z{:,q}(:,mc(q)),model,curren_mode,mode_temp(:,curren_mode),...
                    m_temp{:,curren_mode},P_temp{:,curren_mode},filter.ukf_alpha,...
                    filter.ukf_kappa,filter.ukf_beta,q); %ukf update for this track and this measuremnent
                
            end
            
        end
        log_w_temp = sum(log_qz_temp) + tt_predict{i}.w{:,curren_mode} ;                                                                              %unnormalized updated weights
        tt_update{i}.w{:,curren_mode} = sum(log_q_mode_temp) + log_w_temp;     % unnormalized
        tt_update{i}.m{:,curren_mode} = m_temp{:,curren_mode};
        tt_update{i}.P{:,curren_mode} = P_temp{:,curren_mode};
        for_cost(curren_mode) = sum(log_q_mode_temp) + logsumexp(log_w_temp,[],2);
    end
    for  curren_mode = 1:   size(mode_temp,2) % normalize here
        tt_update{i}.w{:,curren_mode} =  tt_update{i}.w{:,curren_mode} - logsumexp(for_cost,[],2);
        tt_update{i}.mode(:,curren_mode) = logsumexp(tt_update{i}.w{:,curren_mode},[],2);
        tt_update{i}.w{:,curren_mode} = tt_update{i}.w{:,curren_mode} - logsumexp(tt_update{i}.mode(:,curren_mode),[],2);
    end
    
    tt_update{i}.l = tt_predict{i}.l;
    temp = NaN(max_num_sensors,1);
    temp(index_sensors)  = mc;
    tt_update_msqz(i) = logsumexp(for_cost,[],2);

    tt_update{i}.ah = [tt_predict{i}.ah temp];
end

glmb_update.tt = tt_update;

for pidx = 1:length(glmb_update.w)
    glmb_update.I{pidx} = tt_update_newidx(glmb_update.I{pidx});
    glmb_update.w(pidx) = glmb_update.w(pidx) + sum(tt_update_msqz(glmb_update.I{pidx}));
end
glmb_update.w = glmb_update.w - logsumexp(glmb_update.w,[],2);

% Compute cardinality distribution
max_card = max(glmb_update.n);
glmb_update.cdn = zeros(1,max_card+1);
for card = 0:max_card
    idx = (glmb_update.n==card);
    if any(idx)
        glmb_update.cdn(card+1) = exp(logsumexp(glmb_update.w(idx),[],2));
    end
end

% Remove duplicate entries and clean track table
glmb_update = glmb_clean_update(glmb_clean_predict(glmb_update));

end
