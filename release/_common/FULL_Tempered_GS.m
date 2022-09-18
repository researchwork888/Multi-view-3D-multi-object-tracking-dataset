function [assignments,costs]= FULL_Tempered_GS(P0,m,filter)

% Declaration
n1 = size(P0,1);
n2 = size(P0,2);
% Initialization
if m == 0
    currsoln= n1+1:2*n1; %use all missed detections as initial solution
    assignments= currsoln;
    costs= sum(P0([1:n1]+(currsoln-1)*n1));
    return
end

num_iter = filter.H_randomGibss;


assignments= zeros(num_iter,n1);
costs= zeros(num_iter,1)';

currsoln = n1+1:2*n1;
assignments(1,:)= currsoln;
costs(1)=sum(P0([1:n1]+(currsoln-1)*n1));
actual_conditional = zeros(n1,n2);
N_actual_conditional =  zeros(n1,n2);
N_proposal_conditional = zeros(n1,n2);

for var = 1: n1 
    % %
    actual_conditional(var,:) = exp(-P0(var,:));
    
    N_actual_conditional(var,:) = actual_conditional(var,:) ./ sum(actual_conditional(var,:));
    temp_hold =  N_actual_conditional(var,:);
    temp_hold((temp_hold <= eps )) = 0; % zero out what's supposed to be zero
    temp_hold = temp_hold ./ sum(temp_hold); % renormalize, should be the same as N_actual_conditional
    maxx =  max( temp_hold );
    temp_hold(temp_hold == 0) = inf; % to ignore the zero minimum!
    minn = min( temp_hold);
    beta =max( 1-(maxx - minn) , filter.temper_min_beta);
%     N_proposal_conditional(var,:) = N_actual_conditional(var,:).^(beta); % 1-(pd(n))^2
    N_proposal_conditional(var,(N_actual_conditional(var,:)>0)) = N_actual_conditional(var,(N_actual_conditional(var,:)>0)).^(beta);
    %         proposal_conditional(i,:) = actual_conditional(i,:).^(1-pd(i));
    
    N_proposal_conditional(var,:) = N_proposal_conditional(var,:) ./ sum(N_proposal_conditional(var,:));
    % equation 5
    N_proposal_conditional(var,:)=   0.5*N_proposal_conditional(var,:) +   0.5*N_actual_conditional(var,:) ;
    %  %
    unnorm_sel_prob(var) = (N_proposal_conditional(var,currsoln(var))/N_actual_conditional(var,currsoln(var)));

end

flag = 1; % 1 is random scan with liu's modi , 0 is normal Gibbs
for sol= 2:num_iter
    
    sel_prob = unnorm_sel_prob ./ sum(unnorm_sel_prob);
    cdf = cumsum(sel_prob); %
    var = sum(cdf < rand*cdf(end)) + 1;

    tempsamp = N_proposal_conditional(var,:);
    idxold= find(tempsamp>0);
    tempsamp= tempsamp(idxold);
       
    tempsamp_prob = tempsamp ./ sum(tempsamp);
    unnormalized_proposal = tempsamp_prob;
    unnormalized_proposal(idxold == currsoln(var)) =  0 ;
    temp_idxold = find(unnormalized_proposal);
    unnormalized_proposal= unnormalized_proposal(temp_idxold);
    
    denom = 1 - tempsamp_prob(idxold == currsoln(var));
    
    tempsamp_prob2 = unnormalized_proposal ./ denom;
        if flag  
% % % % % %        This is Liu's modification
                    cdf = cumsum(tempsamp_prob2);
%                     if cdf(end) ~= 1 , error('error'); end
                    currsoln_hold  = sum(cdf < rand*cdf(end)) + 1;
                    accep_ratio = (denom)/(1-tempsamp_prob( temp_idxold(currsoln_hold) ) );
                    if rand <= min(1,accep_ratio)
                        currsoln(var) = idxold(temp_idxold(currsoln_hold));
                    else
                        currsoln(var) = currsoln(var);
                    end
            
        else
%        % Uncomment this section for normal Gibbs
                cdf = cumsum(tempsamp); %#ok<UNRCH>
                currsoln(var)  = sum(cdf < rand*cdf(end)) + 1;
                currsoln(var) = idxold(currsoln(var));
        end
%   %%      %%
%     STORE_assignments(sol,:)  = currsoln; 
        unnorm_sel_prob(var) = (N_proposal_conditional(var,currsoln(var))/N_actual_conditional(var,currsoln(var)));
        assignments(sol,:)= currsoln;
        costs(sol)= sum(P0([1:n1]+(currsoln-1)*n1)); %#ok<NBRAK>

        for i = [1:var-1,var+1:n1]
            % %
            
%                 actual_conditional(i,:) = exp(-P0(i,:));

            if (actual_conditional(i,assignments(sol-1,var)) ==  exp(-P0(i,assignments(sol-1,var))) && ...
                    actual_conditional(i,currsoln(var)) == 0) ||  (assignments(sol-1,var) == currsoln(var)) %|| all(currsoln == assignments(sol-1,:))
                
            else
                actual_conditional(i,assignments(sol-1,var)) =  exp(-P0(i,assignments(sol-1,var)));
                actual_conditional(i,currsoln(var))= 0;
                N_proposal_conditional(i,currsoln(var))= 0;
                
                
                N_actual_conditional(i,:) = actual_conditional(i,:) ./ sum(actual_conditional(i,:));
                temp_hold =  N_actual_conditional(i,:);
                temp_hold((temp_hold <= eps )) = 0; % zero out what's supposed to be zero
                temp_hold = temp_hold ./ sum(temp_hold); % renormalize, should be the same as N_actual_conditional
                maxx =  max( temp_hold );
                temp_hold(temp_hold == 0) = inf; % to ignore the zero minimum!
                minn = min( temp_hold);
                beta =max( 1-(maxx - minn) , filter.temper_min_beta);
%                 N_proposal_conditional(i,:) = N_actual_conditional(i,:).^(beta); % 1-(pd(n))^2
                N_proposal_conditional(i,(N_actual_conditional(i,:)>0)) = N_actual_conditional(i,(N_actual_conditional(i,:)>0)).^(beta);

                %         proposal_conditional(i,:) = actual_conditional(i,:).^(1-pd(i));
                
                N_proposal_conditional(i,:) = N_proposal_conditional(i,:) ./ sum(N_proposal_conditional(i,:));
                % equation 5
                    N_proposal_conditional(var,:)=   0.5*N_proposal_conditional(var,:) +   0.5*N_actual_conditional(var,:) ;
                %  %
                unnorm_sel_prob(i) = (N_proposal_conditional(i,currsoln(i))/N_actual_conditional(i,currsoln(i)));
            end
        end
    
    
    
    
end
[C,I,~]= unique(assignments,'rows');
assignments= C;
costs= costs(I);



