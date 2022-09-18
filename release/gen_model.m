function model= gen_model(dataset,random_sensor,occ_model_on,detector,performance_eval)
model.dataset = dataset; 
model.random_sensor = random_sensor;
model.detector = detector;
model.performance_eval = performance_eval;
% basic parameters
model.no_sensors = 4; % number of sensors
model.max_num_sensors = 4;
model.x_dim= 9;   %dimension of state vector
model.xv_dim= 5;   %dimension of process noise

model.z_dim=  repmat(4,1,model.no_sensors)';   %dimension of observation vector
model.zv_dim= 4;   %dimension of observation noise

% modes
model.mode_type{1} = ['Upright'];
model.mode_type{2} = ['Fallen'];

% active proposed pd model
model.occ_model_on = occ_model_on;

% boundaries
model.XMAX = [2.03 6.3];
model.YMAX = [0.00 3.41];
model.ZMAX = [0 3];

% ellipsoid plotting params (for plotting at the end)
model.ellipsoid_n = 10;
model.ellipsoid_gen_n = 6;

% camera positions, image size and room dimensions
model.sensor_pos{1} = [0.21; 3.11; 2.24];
model.sensor_pos{2} = [7.17; 3.34; 2.16];
model.sensor_pos{3} = [7.55; 0.47; 2.16];
model.sensor_pos{4} = [0.21; 1.26; 2.20];
model.imagesize = [1920; 1024];
model.room_dim = [7.67; 3.41; 2.7];

% load camera parameters
load('cam1_cam_mat','cam1_cam_mat');
model.cam1_cam_mat = cam1_cam_mat;
model.cam1_homo = [cam1_cam_mat(:,1:2) cam1_cam_mat(:,4)];

load('cam2_cam_mat','cam2_cam_mat');
model.cam2_cam_mat = cam2_cam_mat;
model.cam2_homo = [cam2_cam_mat(:,1:2) cam2_cam_mat(:,4)];

load('cam3_cam_mat','cam3_cam_mat');
model.cam3_cam_mat =cam3_cam_mat;
model.cam3_homo = [cam3_cam_mat(:,1:2) cam3_cam_mat(:,4)];

load('cam4_cam_mat','cam4_cam_mat');
model.cam4_cam_mat = cam4_cam_mat;
model.cam4_homo = [cam4_cam_mat(:,1:2) cam4_cam_mat(:,4)];

% dynamical model parameters (CV model)
model.T= 1; %sampling period (1 time step / unit time step)
model.A0= [ 1 model.T; 0 1 ]; %transition matrix
model.F = blkdiag(kron(eye(3),model.A0),eye(3));

if strcmp(dataset,'CMC1') || strcmp(dataset,'CMC2') || strcmp(dataset,'CMC3')
    model.sigma_v = 0.035; %
    [model.n_mu(1,1),model.n_std_dev(1,1)]=lognormal_with_mean_one(0.06); %0.06 input is std dev of multiplicative lognormal noise.
    [model.n_mu(2,1),model.n_std_dev(2,1)]=lognormal_with_mean_one(0.02); % input is std dev of multiplicative lognormal noise.
    model.sigma_radius = model.n_std_dev(1,1);
    model.sigma_heig =model.n_std_dev(2,1); 
    
    model.B0(:,:,1) = model.sigma_v*[ (model.T^2)/2; model.T ];
    model.B1(:,:,1) = diag([model.sigma_radius,model.sigma_radius,model.sigma_heig]);
    model.B(:,:,1) = blkdiag(kron(eye(3),model.B0(:,:,1)),model.B1(:,:,1) );
    model.Q(:,:,1) = model.B(:,:,1)*model.B(:,:,1)';%     %process noise covariance
    
    % birth parameters (LMB birth model, single component only)
    model.T_birth= 1;         %no. of LMB birth terms
    model.L_birth= zeros(model.T_birth,1);         %no of Gaussians in each LMB birth term
    model.r_birth= zeros(model.T_birth,1);         %prob of birth for each LMB birth term
    model.w_birth= cell(model.T_birth,1);          %weights of GM for each LMB birth term
    model.m_birth= cell(model.T_birth,1);          %means of GM for each LMB birth term
    model.B_birth= cell(model.T_birth,1);          %std of GM for each LMB birth term
    model.P_birth= cell(model.T_birth,1);
    
    model.L_birth(1)=1;     %no of Gaussians in birth term 1
    if strcmp(dataset,'CMC1'), model.r_birth(1)=log(0.012); else model.r_birth(1)=log(0.22); end %prob of birth 0.005
    model.mode(1,1) = log(1); %  in log domain
    model.w_birth{1}{1}(1,1) = log(1); %weight of Gaussians (in log domain) - must be column_vector
    model.scale = 1; % Change this at your own risk!
    [n_mu_hold,n_std_dev_hold]=lognormal_with_mean_one(0.1); % input is std dev of multiplicative lognormal noise.
    model.m_birth{1}{1}(:,1) =  [  2.3;0.01;1.2  ;0.01;0.825;0; (log(0.3))+n_mu_hold;(log(0.3))+n_mu_hold; (log(0.84))+n_mu_hold] ;
    model.B_birth{1}(:,:,1) = diag([ 0.2;0.1;0.2; 0.1;0.15;0.1; n_std_dev_hold;n_std_dev_hold;n_std_dev_hold]);  %diag([ 0.25;0.1;0.25; 0.1;0.15;0.1; n_std_dev_hold;n_std_dev_hold;n_std_dev_hold]);
    model.P_birth{1}{1}(:,:,1) = model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians
    
    % Markov transition matrix for mode 0 is standing 1 is fall
    model.mode_trans_matrix = [log(0.99) log(0.01);...
        log(0.99) log(0.01)];% probability of survival
elseif strcmp(dataset,'CMC4') || strcmp(dataset,'CMC5')
    % transition for standing to standing
    model.sigma_vz(1) = 0.035; 
    model.sigma_vxy(1) = 0.035;
    model.scale = 1; % Change this at your own risk!
    [model.n_mu(1,1),model.n_std_dev(1,1)]=lognormal_with_mean_one(0.06); %input is std dev of multiplicative lognormal noise. 0.06
    [model.n_mu(2,1),model.n_std_dev(2,1)]=lognormal_with_mean_one(0.02); %input is std dev of multiplicative lognormal noise. 0.02
    model.sigma_radius(1) = model.n_std_dev(1,1);
    model.sigma_heig(1) =model.n_std_dev(2,1);  
    model.B0(:,:,1) = [ (model.T^2)/2; model.T ];
    model.B1(:,:,1) = diag([model.sigma_radius(1),model.sigma_radius(1),model.sigma_heig(1)]);
    model.B(:,:,1) = blkdiag(kron(diag([model.sigma_vxy(1),model.sigma_vxy(1),model.sigma_vz(1)]),model.B0(:,:,1)),model.B1(:,:,1) );
    model.Q(:,:,1) = model.B(:,:,1)*model.B(:,:,1)';%     %process noise covariance
    
    % transition for falling to falling
    model.sigma_vz(2) = 0.035; 
    model.sigma_vxy(2) = 0.035; 
    [model.n_mu(1,2),model.n_std_dev(1,2)]=lognormal_with_mean_one(0.4); % input is std dev of multiplicative lognormal noise.  1
    [model.n_mu(2,2),model.n_std_dev(2,2)]=lognormal_with_mean_one(0.2); % input is std dev of multiplicative lognormal noise.  1
    model.sigma_radius(2) = model.n_std_dev(1,2);
    model.sigma_heig(2) =model.n_std_dev(2,2); 
    model.B0(:,:,2) = [ (model.T^2)/2; model.T ];
    model.B1(:,:,2) = diag([model.sigma_radius(2),model.sigma_radius(2),model.sigma_heig(2)]);
    model.B(:,:,2) = blkdiag(kron(diag([model.sigma_vxy(2),model.sigma_vxy(2),model.sigma_vz(2)]),model.B0(:,:,2)),model.B1(:,:,2));
    model.Q(:,:,2) = model.B(:,:,2)*model.B(:,:,2)';%
    
    % transition from standing to fallen (vice versa)
    model.sigma_vz(3) = 0.07; % 0.07
    model.sigma_vxy(3) = 0.07; % 0.07
    
    [model.n_mu(1,3),model.n_std_dev(1,3)]=lognormal_with_mean_one(0.1); % input is std dev of multiplicative lognormal noise.
    model.sigma_radius(3) = model.n_std_dev(1,3);  
    model.sigma_heig(3) = model.n_std_dev(1,3); 
    model.B0(:,:,3) = [ (model.T^2)/2; model.T ];
    model.B1(:,:,3) = diag([model.sigma_radius(3),model.sigma_radius(3),model.sigma_heig(3)]);
    model.B(:,:,3) = blkdiag(kron(diag([model.sigma_vxy(3),model.sigma_vxy(3),model.sigma_vz(3)]),model.B0(:,:,3)),model.B1(:,:,3) );
    model.Q(:,:,3) = model.B(:,:,3)*model.B(:,:,3)';%     %process noise covariance
    
    % birth parameters (LMB birth model, single component only)
    model.T_birth= 1;         %no. of LMB birth terms
    model.L_birth= zeros(model.T_birth,1);                                          %no of Gaussians in each LMB birth term
    model.r_birth= zeros(model.T_birth,1);                                          %prob of birth for each LMB birth term
    model.w_birth= cell(model.T_birth,1);                                           %weights of GM for each LMB birth term
    model.m_birth= cell(model.T_birth,1);                                           %means of GM for each LMB birth term
    model.B_birth= cell(model.T_birth,1);                                           %std of GM for each LMB birth term
    model.P_birth= cell(model.T_birth,1);                                           %cov of GM for each LMB birth term
    
    model.L_birth(1)=1;                                                             %no of Gaussians in birth term 1
    model.r_birth(1)=log(0.07); % 0.005                                        %prob of birth
    model.mode(1,[1 2]) = [log(0.6) log(0.4)];
    model.w_birth{1}{1}(1,:) = log(1);                                                       %weight of Gaussians - must be column_vector
    model.w_birth{1}{2}(1,:) = log(1);
    [n_mu_hold,n_std_dev_hold]=lognormal_with_mean_one(0.2); % 2.52input is std dev of multiplicative lognormal noise.
    model.m_birth{1}{1}(:,1) =  [  2.3;0;1.2 ;0;0.825;0; (log(0.3))+n_mu_hold;(log(0.3))+n_mu_hold; (log(0.84))+n_mu_hold] ;          %mean of Gaussians  [  2.52;0;0.71 ;0;1.65; 0.25; 0.25]    [  2.52;0;0.71 ;0;0.825;0; 0.25; 0.25;0.825] ;
    model.m_birth{1}{2}(:,1) =  [  2.3;0;1.2 ;0;0.825/2;0; (log(0.84))+n_mu_hold;(log(0.84))+n_mu_hold; (log(0.3))+n_mu_hold] ;
    model.B_birth{1}(:,:,1) =  diag([ 0.2;0.1;0.2; 0.1;0.15;0.1; n_std_dev_hold;n_std_dev_hold;n_std_dev_hold]);              %std of Gaussians diag([ 0.1;0.1;0.1;0.1;0.1;0.1;0.1]); diag([ 0.1;0.1;0.1; 0.1;0.1;0.1; 0.1;0.1;0.1]);
    model.P_birth{1}{1}(:,:,1) = model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians
    model.P_birth{1}{2}(:,:,1) = model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';
    % markov transition matrix for mode 0 is standing 1 is fall
    model.mode_trans_matrix = [log(0.6) log(0.4);...
                                log(0.4) log(0.6)];
    
    if strcmp(dataset,'CMC5') % different calib params due to different set of recording
        load('cam4_cam_mat__','cam4_cam_mat__');
        model.cam4_cam_mat = cam4_cam_mat__;
        model.cam4_homo = [cam4_cam_mat__(:,1:2) cam4_cam_mat__(:,4)];
        load('cam3_cam_mat__','cam3_cam_mat__');
        model.cam3_cam_mat = cam3_cam_mat__;
        model.cam3_homo = [cam3_cam_mat__(:,1:2) cam3_cam_mat__(:,4)];
        
    end
end

model.P_S= log(0.9999999); %9999
model.Q_S= log(1-exp(model.P_S));
% measurement parameters
% mode 0 (standing)
[model.meas_n_mu(1,1),model.meas_n_std_dev(1,1)]=lognormal_with_mean_one(0.1);
[model.meas_n_mu(2,1),model.meas_n_std_dev(2,1)]=lognormal_with_mean_one(0.05);
model.D(:,:,1)= diag([20,20,model.meas_n_std_dev(1,1),model.meas_n_std_dev(2,1)]);
% mode 1 (fallen)
[model.meas_n_mu(1,2),model.meas_n_std_dev(1,2)]=lognormal_with_mean_one(0.05);
[model.meas_n_mu(2,2),model.meas_n_std_dev(2,2)]=lognormal_with_mean_one(0.1);
model.D(:,:,2)= diag([20,20,model.meas_n_std_dev(1,2),model.meas_n_std_dev(2,2)]);

model.R(:,:,1)= model.D(:,:,1)*model.D(:,:,1)';
model.R(:,:,2)= model.D(:,:,2)*model.D(:,:,2)';  %observation noise covariance

% detection probabilities
model.P_D= log(0.99);   %probability of detection in measurements
model.Q_D= log(1-exp(model.P_D)); %probability of missed detection in measurements

lambda_c = 5;
model.lambda_c =  repmat(lambda_c,1,model.no_sensors)'; %poisson average rate of uniform clutter (per scan)
model.range_c= [ 1 1920;1 1024 ;(1) (1920);(1) (1024)]; %uniform clutter region
% model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
range_temp = model.range_c(:,2)-model.range_c(:,1);
range_temp(3:4) = log(range_temp(3:4));
model.pdf_c = repmat( 1/prod(range_temp),1,model.no_sensors)';


end