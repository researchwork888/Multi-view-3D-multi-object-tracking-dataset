function [est,prior_estimator]  = run_filter(model,meas)

% This is the MATLAB code for the implementation of the Generalized Labeled Multi-Bernoulli filter proposed in
% B.-T. Vo, and B.-N. Vo, "Labeled Random Finite Sets and Multi-Object Conjugate Priors," IEEE Trans. Signal Processing, Vol. 61, No. 13, pp. 3460-3475, 2013.
% http://ba-ngu.vo-au.com/vo/VV_XXXXX_TSP16.pdf
% originally proposed in
% B.-T. Vo, and B.-N. Vo, "Labeled Random Finite Sets and Multi-Object Conjugate Priors," IEEE Trans. Signal Processing, Vol. 61, No. 13, pp. 3460-3475, 2013.
% http://ba-ngu.vo-au.com/vo/VV_Conjugate_TSP13.pdf
% and via Murty's algorithm as per
% B.-N. Vo, B.-T. Vo, and D. Phung, "Labeled Random Finite Sets and the Bayes Multi-Target Tracking Filter," IEEE Trans. Signal Processing, Vol. 62, No. 24, pp. 6554-6567, 2014
% http://ba-ngu.vo-au.com/vo/VVP_GLMB_TSP14.pdf
%
% Note 1: no lookahead PHD/CPHD allocation is implemented in this code, a simple proportional weighting scheme is used for readability
% Note 2: the simple example used here is the same as in the CB-MeMBer filter code for a quick demonstration and comparison purposes
% Note 3: more difficult scenarios require more components/hypotheses (thus exec time) and/or a better lookahead
% ---BibTeX entry
% @ARTICLE{GLMB1,
% author={B.-T. Vo and B.-N. Vo
% journal={IEEE Transactions on Signal Processing},
% title={Labeled Random Finite Sets and Multi-Object Conjugate Priors},
% year={2013},
% month={Jul}
% volume={61},
% number={13},
% pages={3460-3475}}
%
% @ARTICLE{GLMB2,
% author={B.-T. Vo and B.-N. Vo and D. Phung},
% journal={IEEE Transactions on Signal Processing},
% title={Labeled Random Finite Sets and the Bayes Multi-Target Tracking Filter},
% year={2014},
% month={Dec}
% volume={62},
% number={24},
% pages={6554-6567}}
%
% @ARTICLE{GLMB3,
% author={B.-N. Vo and B.-T. Vo and H. Hung},
% journal={IEEE Transactions on Signal Processing},
% title={An Efficient Implementation of the Generalized Labeled Multi-Bernoulli Filter},
% year={XXXX},
% month={XXX}
% volume={XX},
% number={XX},
% pages={XXXX-XXXX}}
%---

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);
est.ah= cell(meas.K,1);

%filter parameters
filter.H_upd= 3000;                 %requested number of updated components/hypotheses
filter.H_max= 100;                 %cap on number of posterior components/hypotheses
filter.H_randomGibss = 50;
filter.temper_min_beta = 0.3;
filter.comp_min_retain_per_card = 0;%0

filter.comp_max_num_request= filter.H_upd;                 %requested number of updated components/hypotheses
filter.comp_max_num_retain= filter.H_max;                 %cap on number of posterior components/hypotheses

filter.comp_min_weight_retain = log(1e-50);


% GM pruning and merging parameters
filter.gauss_max_num_per_track = 10;
filter.gauss_elim_threshold = 1e-5;
filter.gauss_merge_threshold = 4;

filter.P_G= 0.9999999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 0;                             %gating on or off 1/0

% UKF parameters
filter.ukf_alpha= 1;                %scale parameter for UKF - choose alpha=1 ensuring lambda=beta and offset of first cov weight is beta for numerical stability
filter.ukf_beta= 2;                 %scale parameter for UKF
filter.ukf_kappa= 2;                %scale parameter for UKF (alpha=1 preferred for stability, giving lambda=1, offset of beta for first cov weight)

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output


est.filter= filter;

%=== Filtering

%initial prior
glmb_update.tt= cell(0,1);      %track table for GLMB (cell array of structs for individual tracks)
glmb_update.w= log(1);               %vector of GLMB component/hypothesis weights
glmb_update.I= {[]};            %cell of GLMB component/hypothesis labels (labels are indices/entries in track table)
glmb_update.n= 0;               %vector of GLMB component/hypothesis cardinalities
glmb_update.cdn= 1;             %cardinality distribution of GLMB (vector of cardinality distribution probabilities)

prior_estimator.labels = [] ;
prior_estimator.thetas = [] ;
prior_estimator.inits =[] ;

% Output variables
est.X = cell(meas.K,1);
est.N = zeros(meas.K,1);
est.L = cell(meas.K,1);
est.mode = cell(meas.K,1);
est.T= {};
est.M= 0;
est.J= []; est.H= {};




if model.random_sensor 
    
    %recursive filtering
    for k=1:meas.K
        if k <= 100 % all cams
            model.no_sensors = 4;
            model.index_sensors = [1:model.no_sensors];
        elseif (k>100 && k <= 150 )
            model.no_sensors =3;
            model.index_sensors = [1 2 3];
        elseif (k>150 && k <= 200 )  % randomize
            model.no_sensors = 3;
            model.index_sensors=randsample(model.no_sensors,3);
            
        elseif (k>200 && k <= 225 )  %2 cams
            model.no_sensors =2;
            model.index_sensors=[1 3];
            
        elseif (k>225 )  % 2 cams
            model.no_sensors =2;
            model.index_sensors=[2 4];
        end
           
        %joint prediction and update
        glmb_update= msjointpredictupdate(glmb_update,model,filter,meas,k);
        H_posterior= length(glmb_update.w);
        
        if strcmp(model.dataset,'CMC4') || strcmp(model.dataset,'CMC5')
            %pruning, merging, capping Gaussian mixture components
            for ntts = 1 : length(glmb_update.tt)
                for n_mode = 1 : size(glmb_update.tt{ntts}.mode,2)
                    [glmb_update.tt{ntts}.w{:,n_mode},glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode}]= glmb_gaus_mode_prune(glmb_update.tt{ntts}.w{:,n_mode},...
                        glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode},filter.gauss_elim_threshold);
                    [glmb_update.tt{ntts}.w{:,n_mode},glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode}]= glmb_gaus_mode_merge(glmb_update.tt{ntts}.w{:,n_mode},...
                        glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode},filter.gauss_merge_threshold);
                    [glmb_update.tt{ntts}.w{:,n_mode},glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode}]= glmb_gaus_cap(glmb_update.tt{ntts}.w{:,n_mode},...
                        glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode},filter.gauss_max_num_per_track);
                end
            end
        end
        
        %pruning and truncation
        glmb_update = glmb_prune(glmb_update,filter);
        H_prune= length(glmb_update.w);
        glmb_update = glmb_truncate(glmb_update,filter);
        H_cap= length(glmb_update.w);
        glmb_update = glmb_clean_update(glmb_update);
        
        % State estimation
        [est.X{k},est.N(k),est.mode{k},est.L{k},idxcmp] = glmb_extract_estimates(glmb_update,model);
        %     [est,idxcmp]= glmb_extract_estimates_recursive(glmb_update,model,meas,est,filter); %[est.X{k},est.N(k),est.L{k}]= extract_estimates(glmb_update,model);
        
        est.save_sensor_index{k} = model.index_sensors ;
        prior_estimator = collect_trajectories(model,glmb_update,idxcmp,prior_estimator) ;
        
        display_diaginfo(glmb_update,k,est,filter,H_posterior,H_posterior,H_prune,H_cap);
        
        % plot estimates on image
        %     plot_estimates_on_image_(model,meas,est)
        [labelcount,hold_label]= countestlabels(meas,est);
        colorarray= makecolorarray(labelcount);
        est.total_tracks= labelcount;
        est.track_list= cell(k,1);
        for k_temp =meas.start :k
            for eidx=1:size(est.X{k_temp},2)
                [idxx,colorarray]  = assigncolor(colorarray,est.L{k_temp}(:,eidx));
                est.track_list{k_temp} = [est.track_list{k_temp} idxx];
            end
        end
        % [Y_track,l_birth,l_death]= extract_tracks(est.X,est.mode,est.track_list,est.total_tracks,meas);
        [est_temp.Y_track,est_temp.l_birth,est_temp.l_death]= extract_tracks(model,est.X,est.mode,est.track_list,est.total_tracks,meas,k);
        est_temp.colorarray= colorarray;
        est_temp.hold_label = hold_label;
        videofig(k, @redraw,model,meas,est_temp,10);% '4' is the fps
        if (k ~= meas.end)||1 % set a breakpoint here to view result
            close all
        end
    end
    
    
    
else
    %recursive filtering
    for k=1:meas.K
        
        model.index_sensors = [1:model.no_sensors];
        %joint prediction and update
        glmb_update= msjointpredictupdate(glmb_update,model,filter,meas,k);
        H_posterior= length(glmb_update.w);
        
        if strcmp(model.dataset,'CMC4') || strcmp(model.dataset,'CMC5')
            %pruning, merging, capping Gaussian mixture components
            for ntts = 1 : length(glmb_update.tt)
                for n_mode = 1 : size(glmb_update.tt{ntts}.mode,2)
                    [glmb_update.tt{ntts}.w{:,n_mode},glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode}]= glmb_gaus_mode_prune(glmb_update.tt{ntts}.w{:,n_mode},...
                        glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode},filter.gauss_elim_threshold);
                    [glmb_update.tt{ntts}.w{:,n_mode},glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode}]= glmb_gaus_mode_merge(glmb_update.tt{ntts}.w{:,n_mode},...
                        glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode},filter.gauss_merge_threshold);
                    [glmb_update.tt{ntts}.w{:,n_mode},glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode}]= glmb_gaus_cap(glmb_update.tt{ntts}.w{:,n_mode},...
                        glmb_update.tt{ntts}.m{:,n_mode},glmb_update.tt{ntts}.P{:,n_mode},filter.gauss_max_num_per_track);
                end
            end
        end
        
        %pruning and truncation
        glmb_update = glmb_prune(glmb_update,filter);
        H_prune= length(glmb_update.w);
        glmb_update = glmb_truncate(glmb_update,filter);
        H_cap= length(glmb_update.w);
        glmb_update = glmb_clean_update(glmb_update);
        
        % State estimation
        [est.X{k},est.N(k),est.mode{k},est.L{k},idxcmp] = glmb_extract_estimates(glmb_update,model);
        %     [est,idxcmp]= glmb_extract_estimates_recursive(glmb_update,model,meas,est,filter); %[est.X{k},est.N(k),est.L{k}]= extract_estimates(glmb_update,model);
        
        est.save_sensor_index{k} = model.index_sensors ;
        prior_estimator = collect_trajectories(model,glmb_update,idxcmp,prior_estimator) ;
        
        display_diaginfo(glmb_update,k,est,filter,H_posterior,H_posterior,H_prune,H_cap);
        
        % plot estimates on image
        %     plot_estimates_on_image_(model,meas,est)
        [labelcount,hold_label]= countestlabels(meas,est);
        colorarray= makecolorarray(labelcount);
        est.total_tracks= labelcount;
        est.track_list= cell(k,1);
        for k_temp =meas.start :k
            for eidx=1:size(est.X{k_temp},2)
                [idxx,colorarray]  = assigncolor(colorarray,est.L{k_temp}(:,eidx));
                est.track_list{k_temp} = [est.track_list{k_temp} idxx];
            end
        end
        % [Y_track,l_birth,l_death]= extract_tracks(est.X,est.mode,est.track_list,est.total_tracks,meas);
        [est_temp.Y_track,est_temp.l_birth,est_temp.l_death]= extract_tracks(model,est.X,est.mode,est.track_list,est.total_tracks,meas,k);
        est_temp.colorarray= colorarray;
        est_temp.hold_label = hold_label;
        videofig(k, @redraw,model,meas,est_temp,10);% '4' is the fps
        if (k ~= meas.end)||1 % set a breakpoint here to view result
            close all
        end
    end

end
%    vsdv
end

function glmb_nextupdate= msjointpredictupdate(glmb_update,model,filter,meas,k)

num_components = length(glmb_update.w);
num_sensors = model.no_sensors;
index_sensors = model.index_sensors;
max_num_sensors = model.max_num_sensors;
num_meas = zeros(num_sensors,1);
for s = 1:num_sensors
    meas.Z{k,s} = eval(['meas.bbs',num2str(index_sensors(s)),'{',num2str(k),'}']);
    model.cam_matrix(:,:,s) = eval(['model.cam',num2str(index_sensors(s)),'_cam_mat']);
    num_meas(s) = size(meas.Z{k,s},2);
    
end
%create birth tracks
num_birth = length(model.r_birth);
tt_birth = cell(num_birth,1);
for tabbidx=1:num_birth
    tt_birth{tabbidx}.m= model.m_birth{tabbidx};                                   %means of Gaussians for birth track
    tt_birth{tabbidx}.P= model.P_birth{tabbidx};                                   %covs of Gaussians for birth track
    tt_birth{tabbidx}.w= model.w_birth{tabbidx};                                %weights of Gaussians for birth track
    tt_birth{tabbidx}.mode = model.mode(tabbidx,:); 
    tt_birth{tabbidx}.l= [k;tabbidx];                                              %track label
    tt_birth{tabbidx}.ah= [];                                                      %track association history (empty at birth)
end

%create surviving tracks - via time prediction (single target CK)
num_survive = length(glmb_update.tt);
tt_survive = cell(num_survive,1);

for tabsidx=1:num_survive
    [log_modetemp_predict,log_wtemp_predict,m_predict,P_predict] = ...
        kalman_predict_multiple_IC_MSGLMB(model,glmb_update.tt{tabsidx}.mode,...
        glmb_update.tt{tabsidx}.w,glmb_update.tt{tabsidx}.m,glmb_update.tt{tabsidx}.P);
    tt_survive{tabsidx}.m= m_predict;                                                                                   %means of Gaussians for surviving track
    tt_survive{tabsidx}.P= P_predict;                                                                                   %covs of Gaussians for surviving track
    tt_survive{tabsidx}.w= log_wtemp_predict; % make sure is in the log domain                                                                      %weights of Gaussians for surviving track
    tt_survive{tabsidx}.mode = log_modetemp_predict; % make sure is in the log domain
    tt_survive{tabsidx}.l= glmb_update.tt{tabsidx}.l;                                                                       %track label
    tt_survive{tabsidx}.ah= glmb_update.tt{tabsidx}.ah;                                                                     %track association history (no change at prediction)
end

%create predicted tracks - concatenation of birth and survival
glmb_predict.tt= cat(1,tt_birth,tt_survive);                                                                                %copy track table back to GLMB struct
num_tracks_predict = num_birth + num_survive;

%gating by tracks
if filter.gate_flag
    for tabidx=1:length(glmb_predict.tt)
        for s=1:model.N_sensors
%             glmb_predict.tt{tabidx}.gatemeas{s}= gate_msmeas_gms_idx(meas.Z{k,s},filter.gamma(s),model,s,glmb_predict.tt{tabidx}.m,glmb_predict.tt{tabidx}.P);
            glmb_predict.tt{tabidx}.gatemeas{s}= gate_msmeas_ukf_idx(meas.Z{k,s},filter.gamma(s),model,s,glmb_predict.tt{tabidx}.m,glmb_predict.tt{tabidx}.P,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);

        end
    end
else
    for tabidx=1:length(glmb_predict.tt)
        for s=1:num_sensors
            glmb_predict.tt{tabidx}.gatemeas{s}= 1:size(meas.Z{k,s},2);
        end
    end
end

%precalculation loop for average survival/death probabilities
avps = [(model.r_birth) ; zeros(num_survive,1)];
avqs = [log(1-exp(model.r_birth)) ; zeros(num_survive,1)];
for t = 1:num_survive
    temp_ps = compute_pS(model,glmb_update.tt{t}.mode,glmb_update.tt{t}.w,glmb_update.tt{t}.m,glmb_update.tt{t}.l,k);
    avps(num_birth+t) = (temp_ps);
    avqs(num_birth+t) = log(1-exp(temp_ps));
end

%create updated tracks (single target Bayes update)
m= zeros(num_sensors,1);
for s= 1:num_sensors
    m(s)= size(meas.Z{k,s},2);     %number of measurements
end

if strcmp(model.dataset,'CMC1') || strcmp(model.dataset,'CMC2') || strcmp(model.dataset,'CMC3')
    %nested for loop over all predicted tracks and sensors - slow way - Kalman updates on the same prior recalculate all quantities
    allcostc = cell(num_sensors,1);
    jointcostc = cell(num_sensors,1);
    avpp = zeros(num_tracks_predict,num_sensors);
    for s = 1:num_sensors
        allcostc{s} = -inf*ones(length(glmb_predict.tt),1+num_meas(s));
        for t = 1:length(glmb_predict.tt)
            for i = glmb_predict.tt{t}.gatemeas{s}
                for curren_mode = 1:  size(glmb_predict.tt{t}.mode,2)
                    [~,log_qz_temp,~,~] = ukf_update_per_sensor(meas.Z{k,s}(:,i),model,curren_mode,glmb_predict.tt{t}.mode(:,curren_mode),...
                        glmb_predict.tt{t}.m{:,curren_mode},glmb_predict.tt{t}.P{:,curren_mode},filter.ukf_alpha,...
                        filter.ukf_kappa,filter.ukf_beta,s); %ukf update for this track and this measuremnent
                    w_temp = log_qz_temp + glmb_predict.tt{t}.w{:,curren_mode};
                    for_cost(curren_mode) =  logsumexp(w_temp,[],2);
                end
                
                allcostc{s}(t,1+i) = logsumexp(for_cost,[],2); %predictive likelihood
            end
        end
    end
else
    %nested for loop over all predicted tracks and sensors - slow way - Kalman updates on the same prior recalculate all quantities
    allcostc = cell(num_sensors,1);
    jointcostc = cell(num_sensors,1);
    avpp = zeros(num_tracks_predict,num_sensors);
    for s = 1:num_sensors
        allcostc{s} = -inf*ones(length(glmb_predict.tt),1+num_meas(s));
        for t = 1:length(glmb_predict.tt)
            for i = glmb_predict.tt{t}.gatemeas{s}
                for curren_mode = 1:  size(glmb_predict.tt{t}.mode,2)
                    [log_q_mode,log_qz_temp,~,~] = ukf_update_per_sensor(meas.Z{k,s}(:,i),model,curren_mode,glmb_predict.tt{t}.mode(:,curren_mode),...
                        glmb_predict.tt{t}.m{:,curren_mode},glmb_predict.tt{t}.P{:,curren_mode},filter.ukf_alpha,...
                        filter.ukf_kappa,filter.ukf_beta,s); %ukf update for this track and this measuremnent
                    w_temp = log_qz_temp + glmb_predict.tt{t}.w{:,curren_mode};
                    for_cost(curren_mode) = log_q_mode + logsumexp(w_temp,[],2);
                end
                
                allcostc{s}(t,1+i) = logsumexp(for_cost,[],2); %predictive likelihood
            end
        end
    end
end

%gated measurement index matrix
gatemeasidxs = cell(num_sensors,1);
gatemeasidxc = cell(num_sensors,1);
for s = 1:num_sensors
    gatemeasidxs{s} = zeros(num_tracks_predict,num_meas(s));
    for t = 1:num_tracks_predict
        gatemeasidxs{s}(t,1:length(glmb_predict.tt{t}.gatemeas{s})) = glmb_predict.tt{t}.gatemeas{s};
    end
    gatemeasidxc{s} = (gatemeasidxs{s} > 0);
end

total_comp_request = filter.H_upd;
lambda_c = model.lambda_c(index_sensors);
pdf_c = model.pdf_c(index_sensors);
weights = glmb_update.w;
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
    num_exists = length(glmb_update.I{pidx});
    num_tracks(pidx) = num_birth + num_exists;
    tindices{pidx} = [1:num_birth num_birth+glmb_update.I{pidx}'];
    avps_comp{pidx} = avps(tindices{pidx});
    avqs_comp{pidx} = avqs(tindices{pidx});
    mindices{pidx} = cell(1,num_sensors);
    neglogcostc_tempd{pidx} = cell(1,num_sensors);
    % applying pd/occlusion model
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
        
        num_comp_request(pidx) = round(exp((log(total_comp_request)+0.5*weights(pidx)-log_sum_sqrt_weights)));
        %       Detection_aka_occlusion_Model
        avpd_comp{pidx}{s}  = Detection_aka_occlusion_Model(extract,model,s,pidx,num_comp_request(pidx),filter);
        avqd_comp{pidx}{s} = log( 1 - exp(avpd_comp{pidx}{s}) );
        avpd_comp{pidx}{s} = (avpd_comp{pidx}{s});
        jointcostc{s} = -inf*ones(length(glmb_predict.tt),1+num_meas(s));
        
        jointcostc{s}(tindices{pidx},mindices{pidx}{s}) = [ avps_comp{pidx}+avqd_comp{pidx}{s}  repmat(avps_comp{pidx}+avpd_comp{pidx}{s} ,[1 num_meas(s)]) + ...
            allcostc{s}(tindices{pidx},mindices{pidx}{s}(:,2:end)) - log(model.lambda_c(s)*model.pdf_c(s))];
        
        avpp(:,s) = logsumexp(jointcostc{s}(:,1:end),[],2);
        neglogcostc_tempd{pidx}{s} = -jointcostc{s}(tindices{pidx},mindices{pidx}{s});
    end
    avpp_comp{pidx} = avpp(tindices{pidx},:);
    
end

%component updates via fast proposal generation
runidx= 1;
tt_update_parent= [];
tt_update_currah= [];
tt_update_linidx= [];

for pidx=1:num_components
    cpreds= length(glmb_predict.tt);
    nbirths= model.T_birth;
    nexists= length(glmb_update.I{pidx});
    ntracks= nbirths + nexists;
    
    dcost = avqs_comp{pidx};
    scost = sum(avpp_comp{pidx},2);
    dprob = dcost - logsumexp(dcost,scost); % ratio  a / a+b
    neglogdprob = -dprob;
    neglogdcost = -dcost;
    if 0
        [uasses]= gibbs_multisensor_approx_cheap(neglogdcost,neglogdprob,neglogcostc_tempd{pidx}',num_comp_request(pidx)); %#ok<UNRCH>
        uasses(uasses<0)= -inf;         
    else
        if num_comp_request(pidx) > 0               
            uasses= RANDOMSCAN_gibbs_multisensor_approx_matlab(neglogdprob',...
                neglogdcost',...
                neglogcostc_tempd{pidx},...
                num_comp_request(pidx),...
                filter);
            uasses = double(uasses+1);
            uasses(uasses<1) = -inf;
        else
            uasses = zeros(0,num_tracks(pidx),num_sensors);
        end
        
    end
    %set not born/track deaths to -inf assignment
    for s = 1:num_sensors
        tempslice = uasses(:,:,s);
        tempslice(tempslice>0) = mindices{pidx}{s}(tempslice(tempslice>0))-1;
        uasses(:,:,s) = tempslice;
    end
    
    
    %generate corrresponding jointly predicted/updated hypotheses/components
    for hidx=1:size(uasses,1)
        update_hypcmp_tmp= reshape(uasses(hidx,:,:),[ntracks num_sensors]); 
        off_idx= update_hypcmp_tmp(:,1)<0; 
        aug_idx= [tindices{pidx}(:) 1+update_hypcmp_tmp]; 
        mis_idx= (update_hypcmp_tmp==0); 
        det_idx= (update_hypcmp_tmp>0); 
        % applying pd/occlusion model

        extract =[];
        for_occlu_idx = find(~off_idx);
        for a = 1 : size(tindices{pidx}(~off_idx),2)
            extract{a}.m_pd = glmb_predict.tt{tindices{pidx}(for_occlu_idx(a))}.m;
            extract{a}.P_pd =  glmb_predict.tt{tindices{pidx}(for_occlu_idx(a))}.P;
            extract{a}.w_pd =  glmb_predict.tt{tindices{pidx}(for_occlu_idx(a))}.w;
            extract{a}.mode_pd = glmb_predict.tt{tindices{pidx}(for_occlu_idx(a))}.mode;
        end
        local_avpdm = ([avpd_comp{pidx}{:}]);
        local_avqdm = ([avqd_comp{pidx}{:}]);
        
        repeated_lambda_c= repmat(lambda_c',[length(tindices{pidx}) 1]); 
        repeated_pdf_c= repmat(pdf_c',[length(tindices{pidx}) 1]); 
        update_hypcmp_idx=zeros(1,length(off_idx)+length(~off_idx));
        update_hypcmp_idx(off_idx)= -inf;
        update_hypcmp_idx(~off_idx)= ndsub2ind([cpreds; 1+m],aug_idx(~off_idx,:));
        num_trk= sum(update_hypcmp_idx>0);
        glmb_nextupdate.w(runidx)=  ... %-sum(lambda_c) + ...
            num_meas'*log(lambda_c.*pdf_c) + ...
            weights(pidx) + ...
            sum(avps_comp{pidx}(~off_idx)) + ...
            sum(avqs_comp{pidx}(off_idx)) + ...
            sum(sum(local_avpdm(det_idx))) + ...
            sum(sum(local_avqdm(mis_idx))) - ...
            sum(sum(log(repeated_lambda_c(det_idx).*repeated_pdf_c(det_idx))));
        
        if num_trk>0
            glmb_nextupdate.I{runidx}= length(tt_update_parent)+1:length(tt_update_parent)+num_trk;
        end                                                                                              %hypothesis/component tracks (via indices to track table)
        glmb_nextupdate.n(runidx)= num_trk;                                                                                             %hypothesis/component cardinality
        runidx= runidx+1;
        
        tt_update_parent= cat(1,tt_update_parent,tindices{pidx}(~off_idx)');
        tt_update_currah= cat(1,tt_update_currah,update_hypcmp_tmp(~off_idx,:));
        tt_update_linidx= cat(1,tt_update_linidx,update_hypcmp_idx(update_hypcmp_idx>0)');   
    end
end

if strcmp(model.dataset,'CMC1') || strcmp(model.dataset,'CMC2') || strcmp(model.dataset,'CMC3')
 
    %component updates via posterior weight correction (including generation of track table)
    [tt_update_allkey,tt_update_oldidx,tt_update_newidx]= unique(tt_update_linidx);
    tt_update_msqz= zeros(length(tt_update_allkey),1);
    tt_update= cell(length(tt_update_allkey),1);
    Z = meas.Z(k,:);
    clear i
    for tabidx=1:length(tt_update_allkey)
        oldidx= tt_update_oldidx(tabidx);
        preidx= tt_update_parent(oldidx);
        meascomb= tt_update_currah(oldidx,:);
        mc = meascomb' ;
        m_temp = glmb_predict.tt{preidx}.m ;
        P_temp  = glmb_predict.tt{preidx}.P ;
        mode_temp = glmb_predict.tt{preidx}.mode;
        
        for curren_mode = 1 :  size(mode_temp,2)
            log_qz_temp= zeros(model.no_sensors,size(glmb_predict.tt{preidx}.w{:,curren_mode},2));
            for q = 1 : size(mc,1)
                if mc(q) == 0
                    log_qz_temp(q) = 0;
                    m_temp{:,curren_mode} = m_temp{:,curren_mode};
                    P_temp{:,curren_mode} = P_temp{:,curren_mode};
                else
                    [~,log_qz_temp(q,:),m_temp{:,curren_mode},P_temp{:,curren_mode}] =...
                        ukf_update_per_sensor(Z{:,q}(:,mc(q)),model,curren_mode,mode_temp(:,curren_mode),...
                        m_temp{:,curren_mode},P_temp{:,curren_mode},filter.ukf_alpha,...
                        filter.ukf_kappa,filter.ukf_beta,q); %ukf update for this track and this measuremnent
                    
                end
                
            end
            log_w_temp = sum(log_qz_temp) + glmb_predict.tt{preidx}.w{:,curren_mode} ;                                                                              %unnormalized updated weights
            tt_update{tabidx}.w{:,curren_mode} = log_w_temp - logsumexp(log_w_temp,[],2);     % unnormalized
            tt_update{tabidx}.m{:,curren_mode} = m_temp{:,curren_mode};
            tt_update{tabidx}.P{:,curren_mode} = P_temp{:,curren_mode};
            for_cost(curren_mode) = logsumexp(log_w_temp,[],2);
            tt_update{tabidx}.mode(:,curren_mode) = log(1);
        end
        tt_update{tabidx}.l = glmb_predict.tt{preidx}.l;
        temp = NaN(max_num_sensors,1);
        temp(index_sensors)  = mc;
        tt_update{tabidx}.ah = [glmb_predict.tt{preidx}.ah temp];
        tt_update_msqz(tabidx) = logsumexp(for_cost,[],2);
    end

else
    %component updates via posterior weight correction (including generation of track table)
    [tt_update_allkey,tt_update_oldidx,tt_update_newidx]= unique(tt_update_linidx);
    tt_update_msqz= zeros(length(tt_update_allkey),1);
    tt_update= cell(length(tt_update_allkey),1);
    Z = meas.Z(k,:);
    clear i
    for tabidx=1:length(tt_update_allkey)
        oldidx= tt_update_oldidx(tabidx);
        preidx= tt_update_parent(oldidx);
        meascomb= tt_update_currah(oldidx,:);
        mc = meascomb' ;
        m_temp = glmb_predict.tt{preidx}.m ;
        P_temp  = glmb_predict.tt{preidx}.P ;
        mode_temp = glmb_predict.tt{preidx}.mode;
        
        for curren_mode = 1 :  size(mode_temp,2)
            log_qz_temp= zeros(model.no_sensors,size(glmb_predict.tt{preidx}.w{:,curren_mode},2));
            log_q_mode_temp= zeros(model.no_sensors,size(glmb_predict.tt{preidx}.mode(:,curren_mode),2));
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
            log_w_temp = sum(log_qz_temp) + glmb_predict.tt{preidx}.w{:,curren_mode} ;                                                                              %unnormalized updated weights
            tt_update{tabidx}.w{:,curren_mode} = sum(log_q_mode_temp) + log_w_temp;     % unnormalized
            tt_update{tabidx}.m{:,curren_mode} = m_temp{:,curren_mode};
            tt_update{tabidx}.P{:,curren_mode} = P_temp{:,curren_mode};
            for_cost(curren_mode) = sum(log_q_mode_temp) + logsumexp(log_w_temp,[],2);
        end
        for  curren_mode = 1:   size(mode_temp,2) % normalize here
            tt_update{tabidx}.w{:,curren_mode} =  tt_update{tabidx}.w{:,curren_mode} - logsumexp(for_cost,[],2);
            tt_update{tabidx}.mode(:,curren_mode) = logsumexp(tt_update{tabidx}.w{:,curren_mode},[],2);
            tt_update{tabidx}.w{:,curren_mode} = tt_update{tabidx}.w{:,curren_mode} - logsumexp(tt_update{tabidx}.mode(:,curren_mode),[],2);
        end
        
        tt_update{tabidx}.l = glmb_predict.tt{preidx}.l;
        temp = NaN(max_num_sensors,1);
        temp(index_sensors)  = mc;
        tt_update{tabidx}.ah = [glmb_predict.tt{preidx}.ah temp];
        tt_update_msqz(tabidx) = logsumexp(for_cost,[],2);
    end
end
glmb_nextupdate.tt = tt_update;

for pidx = 1:length(glmb_nextupdate.w)
    glmb_nextupdate.I{pidx} = tt_update_newidx(glmb_nextupdate.I{pidx});
    glmb_nextupdate.w(pidx) = glmb_nextupdate.w(pidx) + sum(tt_update_msqz(glmb_nextupdate.I{pidx}));
end
glmb_nextupdate.w = glmb_nextupdate.w - logsumexp(glmb_nextupdate.w,[],2);

% Compute cardinality distribution
max_card = max(glmb_nextupdate.n);
glmb_nextupdate.cdn = zeros(1,max_card+1);
for card = 0:max_card
    idx = (glmb_nextupdate.n==card);
    if any(idx)
        glmb_nextupdate.cdn(card+1) = exp(logsumexp(glmb_nextupdate.w(idx),[],2));
    end
end

% Remove duplicate entries and clean track table
glmb_nextupdate = glmb_clean_update(glmb_clean_predict(glmb_nextupdate));


end

function [X_track,k_birth,k_death]= extract_tracks(model,X,mode,track_list,total_tracks,meas,k_end)

K= size(X,1); 
x_dim= model.x_dim; 
x_dim = x_dim + 1; % for mode 
% kkk=K-1; 
% while x_dim==0,
%     x_dim= size(X{kkk},1); 
%     kkk= kkk-1; 
% end
X_track= NaN(x_dim,K,total_tracks);
k_birth= zeros(total_tracks,1);
k_death= zeros(total_tracks,1);

max_idx= 0;
for kk=meas.start:k_end
    if ~isempty(X{kk})
        X_track(:,kk,track_list{kk})= [X{kk}; mode{kk}];
    end
    if max(track_list{kk})> max_idx %new target born?
        idx= find(track_list{kk}> max_idx);
        k_birth(track_list{kk}(idx))= kk;
    end
    if ~isempty(track_list{kk}), max_idx= max(track_list{kk}); end
    k_death(track_list{kk})= kk;
end

end


function display_diaginfo(glmb,k,est,filter,H_predict,H_posterior,H_prune,H_cap)
if ~strcmp(filter.run_flag,'silence')
    disp([' time= ',num2str(k),...
        ' #eap cdn=' num2str((0:(length(glmb.cdn)-1))*glmb.cdn(:)),...
        ' #var cdn=' num2str((0:(length(glmb.cdn)-1)).^2*glmb.cdn(:)-((0:(length(glmb.cdn)-1))*glmb.cdn(:))^2,4),...
        ' #est card=' num2str(est.N(k),4),...
        ' #comp pred=' num2str(H_predict,4),...
        ' #comp post=' num2str(H_posterior,4),...
        ' #comp updt=',num2str(H_cap),...
        ' #trax updt=',num2str(length(glmb.tt),4)   ]);
end
end

