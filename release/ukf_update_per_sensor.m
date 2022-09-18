function [log_q_mode, log_qz_update,m_update,P_update] = ukf_update_per_sensor(z,model,mode,log_w_mode,m,P,alpha,kappa,beta,s)

plength= size(m,2);
zlength= size(z,2);
% theta = m_1(model.x_dim+1,:);
% qz_update= zeros(plength,zlength);
% m_update = zeros(model.x_dim,plength,zlength);
% P_update = zeros(model.x_dim,model.x_dim,plength);

for idxp=1:plength
    m_hold = m([1 3 5],idxp);
    ch1 = m_hold(1) > model.XMAX(1) && m_hold(1) < model.XMAX(2);
    ch2 = m_hold(2) > model.YMAX(1) && m_hold(2) < model.YMAX(2);
    ch3 = m_hold(3) > model.ZMAX(1) && m_hold(3) < model.ZMAX(2);
    if ~(ch1 && ch2 && ch3)
        m_temp = m(:,idxp);
        P_temp = P(:,:,idxp);
        log_qz_temp = log(eps);
    else
        [log_qz_temp,m_temp,P_temp] = ukf_update_single(z,model,m(:,idxp),P(:,:,idxp),alpha,kappa,beta,s,mode);
    end
    
    log_qz_update(:,idxp)   = log_qz_temp;
    m_update(:,idxp,:) = m_temp;
    P_update(:,:,idxp) = P_temp;
end

if mode == 1
    ratio = z(4)/z(3);
    prob = exp(1*(ratio-1));%1
    log_q_mode = log(prob) + log_w_mode; 
    
end

if mode == 2
    ratio = z(4)/z(3);
    prob = exp(-1*(ratio-1));%1
    log_q_mode = log(prob) + log_w_mode; 
    
end



function [log_qz_temp,m_temp,P_temp] = ukf_update_single(z,model,m,P,alpha,kappa,beta,s,mode)

z([3 4],:)= log(z([3 4],:));
[X_ukf,u]= ut( [m; zeros(model.z_dim(s),1) ], blkdiag(P,model.R(:,:,mode)), alpha, kappa );
% [X_ukf,u]= ut( m, P, alpha, kappa );

Z_pred= gen_observation_fn_v2( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.z_dim(s)),:) ,s,mode);
% Z_pred= gen_observation_fn( model, X_ukf(1:model.x_dim,:), X_ukf((model.x_dim+1):(model.x_dim+model.z_dim(s)),:),s,1);
eta= Z_pred*u(:);
S_temp= Z_pred- repmat(eta,[1 length(u)]);
u(1)= u(1)+(1-alpha^2+beta);
S= S_temp*diag(u)*S_temp';
Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
G_temp= X_ukf(1:model.x_dim,:)- repmat(m,[1 length(u)]);
G= G_temp*diag(u)*S_temp';
K  = G*iS;
% qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - ...
% 0.5*dot(z-repmat(eta,[1 size(z,2)]),iS*(z-repmat(eta,[1 size(z,2)]))))'; % real domain
log_qz_temp = -0.5*(size(z,1)*log(2*pi) + log(det_S) + dot((z-repmat(eta,[1 size(z,2)])),iS*(z-repmat(eta,[1 size(z,2)])),1)); % log domain
% iou_likelihood_vector = iou(convert_bbs_to_bbox(z),convert_bbs_to_bbox(eta));
%              if ~isreal(iou_likelihood_vector)
%                  iou_likelihood_vector = 1;
%              end
% %         end
%         qz_temp= iou_likelihood_vector+eps;
m_temp = repmat(m,[1 size(z,2)]) + K*(z-repmat(eta,[1 size(z,2)]));
P_temp = P- G*iS*G';



