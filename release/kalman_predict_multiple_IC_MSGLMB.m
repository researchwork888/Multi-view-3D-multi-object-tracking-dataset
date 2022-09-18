function [mode_predict,w_predict, m_predict, P_predict] = kalman_predict_multiple_IC_MSGLMB(model,log_mode,log_w,m,P)

tlength = 1;% only one gaussian component for transition
mode_trans_matrix = model.mode_trans_matrix;
mode_len = length(log_mode);
% if 1
plength = cell(1,mode_len);

m_predict = cell(1,mode_len);
P_predict = cell(1,mode_len);
w_predict = cell(1,mode_len);

for ct = 1 : tlength  %only one gaussian component for transition
    for next_mode = 1 : mode_len % transition loop  (0 is standing / 1 is fallen)
        for curren_mode = 1 : mode_len
            if curren_mode == next_mode % transition to same mode
                if curren_mode == 1
                    index_for_transition = 1;
                    offset = zeros(size(m{curren_mode})); % this is to offset the normal mean because of lognormal multiplicative noise.
                    offset([7 8 9],:) = offset([7 8 9],:)  + [model.n_mu(1,1);model.n_mu(1,1);model.n_mu(2,1)];
                elseif curren_mode == 2
                    index_for_transition = 2;
                    offset = zeros(size(m{curren_mode})); % this is to offset the normal mean because of lognormal multiplicative noise.
                    offset([7 8 9],:) = offset([7 8 9],:)  + [model.n_mu(1,2);model.n_mu(1,2);model.n_mu(2,2)];
                    
                end
                plength{curren_mode} = size(m{curren_mode},2);
                
                m_per_mode = m{curren_mode} + offset;
                P_per_mode = P{curren_mode};
                w_per_mode = log_w{curren_mode};
                mode_predict_temp(curren_mode) = (mode_trans_matrix(curren_mode,next_mode)) + (log_mode(curren_mode));
                for idxp=1:plength{curren_mode}
                    [m_temp,P_temp] = kalman_predict_single(model.F(:,:,1),model.Q(:,:,index_for_transition),m_per_mode(:,idxp),P_per_mode(:,:,idxp));
                    m_predict{next_mode} = cat(2,m_predict{next_mode},m_temp);
                    P_predict{next_mode} = cat(3,P_predict{next_mode},P_temp);
                    w_predict{next_mode} = cat(2,w_predict{next_mode},w_per_mode(idxp));
                end
            else % transition to diffferent mode
                index_for_transition = 3;
                plength{curren_mode} = size(m{curren_mode},2);
                offset = zeros(size(m{curren_mode})); % this is to offset the normal mean because of lognormal multiplicative noise.
                offset([7 8 9],:) = offset([7 8 9],:)  + [model.n_mu(1,3);model.n_mu(1,3);model.n_mu(1,3)];
                if all([curren_mode next_mode] == [1 2]) && 0
                    m_per_mode = m{curren_mode};
                    
                    temp = max(exp(m{curren_mode}([7,8])));
                    m_per_mode([9]) = log(temp);
                    m_per_mode([7,8]) = m{curren_mode}(9);
                    m_per_mode = m_per_mode + offset;
                else
                    
                    
                    
                    m_per_mode = m{curren_mode} + offset;
                end
                
                P_per_mode = P{curren_mode};
                w_per_mode = log_w{curren_mode};
                
                mode_predict_temp(curren_mode) = (mode_trans_matrix(curren_mode,next_mode)) + (log_mode(curren_mode));
                for idxp=1:plength{curren_mode}
                    [m_temp,P_temp] = kalman_predict_single(model.F(:,:,1),model.Q(:,:,index_for_transition),m_per_mode(:,idxp),P_per_mode(:,:,idxp));
                    m_predict{next_mode} = cat(2,m_predict{next_mode},m_temp);
                    P_predict{next_mode} = cat(3,P_predict{next_mode},P_temp);
                    w_predict{next_mode} = cat(2,w_predict{next_mode},w_per_mode(idxp));
                end
            end
            
        end
        mode_predict(:,next_mode) = logsumexp(mode_predict_temp,[],2);
        %             w_predict{next_mode} = w_predict{next_mode}./sum(w_predict{next_mode}); % real domain
        w_predict{next_mode} = w_predict{next_mode} -  logsumexp(w_predict{next_mode},[],2); % log domain
        
    end
end

% else
%     mode_predict = log(1);% weight in log domain
%     plength= size(m,2);
%     w_predict = log_w;
%     m_predict = zeros(size(m{1}));
%     P_predict = zeros(size(P));
%     offset = zeros(size(m)); % this is to offset the normal mean because of lognormal multiplicative noise.
%
%     offset([7 8 9],:) = offset([7 8 9],:)  + [model.n_mu(1,1);model.n_mu(1,1);model.n_mu(2,1)];
%
%     m_per_mode = m + offset;
%     P_per_mode = P;
%     w_per_mode = w;
%     for idxp=1:plength
%
%         [m_temp,P_temp] = kalman_predict_single(model.F,model.Q(:,:,1),m_per_mode(:,idxp),P_per_mode(:,:,idxp));
%         m_predict(:,idxp) = m_temp;
%         P_predict(:,:,idxp) = P_temp;
%     end
%
% end

if size(w_predict,2) == 1
    w_predict = w_predict' ; % make sure it is along dimension 2.
end

end

function [m_predict,P_predict] = kalman_predict_single(F,Q,m,P)

m_predict = F*m ;
P_predict = Q+F*P*F';


end