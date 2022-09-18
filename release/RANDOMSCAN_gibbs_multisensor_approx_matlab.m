function [uasses] = RANDOMSCAN_gibbs_multisensor_approx_matlab(neglogdprob_tempd,neglogdcost, neglogcostc_tempd, num_comp_request, filter,seed)

assign_costs = neglogcostc_tempd; 

num_tracks = size(neglogdprob_tempd,2);
num_sensors = size(assign_costs,2);
num_samples = filter.H_randomGibss; 

dprob_temp= exp(-neglogdprob_tempd); % changing
dcost =  (-neglogdcost); % fixed 

currsoln = zeros(num_tracks,num_sensors);
assignments = zeros(num_samples,num_tracks,num_sensors);
actual_conditional= cell(num_sensors,1);
N_actual_conditional=cell(num_sensors,1);
N_proposal_conditional=cell(num_sensors,1);
for_death_cost_log_P0 = cell(num_sensors,1);
% sol = 1;
for i = 1 : num_tracks
    for s = 1: num_sensors
        P0 = assign_costs{s};
        % for dprob
        for_death_cost_log_P0{s} = -P0;
        store_for_deathcost(i,s) = log(sum(exp(for_death_cost_log_P0{s}(i,:))));


    end
    scost = sum(store_for_deathcost(i,:));
    dprob_temp(i) = dcost(i) - log(exp(dcost(i)) + exp(scost));
%     dprob_temp(i) =  dcost(i);
%     dprob_temp(i) = 
    for s = 1 : num_sensors
        P0 = assign_costs{s};
        if currsoln(i,s) == -1
            if s == 1
                QQ = exp(dprob_temp(i));
                PP = 1-QQ;
                Prob_DEATH = [PP QQ];
                maxx =  max(  Prob_DEATH );minn = min( Prob_DEATH);
                beta = max( 1-(maxx - minn) ,   filter.temper_min_beta);
                g_Prob_DEATH = Prob_DEATH.^(beta);
                g_Prob_DEATH = g_Prob_DEATH ./ sum(g_Prob_DEATH);
                g_dprob(i,:) = g_Prob_DEATH;
                TEMP_hold(i,s) = (g_Prob_DEATH(2)/Prob_DEATH(2));
            else
                TEMP_hold(i,s) = 1;
            end
        else
            actual_conditional{s}(i,:) = exp(-P0(i,:));
%             if assignments(sol,i,s)>0 
%                 actual_conditional{s}(i,assignments(sol,i,s)) = exp(-P0(i,assignments(sol,i,s)));
%                 
%             end
%             if currsoln(i,s) > 0 
%                 actual_conditional{s}(i,currsoln(i,s)) = 0;
%             end
%             actual_conditional{s}(i,(actual_conditional{s}(i,:) <= eps )) = 0;
%             actual_conditional{s}(i,indx_tempp{s})= 0;
            
            N_actual_conditional{s}(i,:) =  actual_conditional{s}(i,:) ./ sum(actual_conditional{s}(i,:));
            temp_hold =  N_actual_conditional{s}(i,:);
            temp_hold(temp_hold == 0) = inf;
            maxx =  max(  temp_hold );minn = min( temp_hold);
            beta = max( 1-(maxx - minn) , filter.temper_min_beta);
            N_proposal_conditional{s}(i,:) = N_actual_conditional{s}(i,:).^(beta); % 1-(pd(n))^2
            
            N_proposal_conditional{s}(i,:) = N_proposal_conditional{s}(i,:) ./ sum(N_proposal_conditional{s}(i,:));
            % equation 5
            %                 N_proposal_conditional{s}(var,:)=   0.5*N_proposal_conditional{s}(var,:) +   0.5*N_actual_conditional{s}(var,:) ;
            
            if s == 1
                QQ = exp(dprob_temp(i));
                PP = 1-QQ;
                Prob_DEATH = [PP QQ];
                maxx =  max(  Prob_DEATH );minn = min( Prob_DEATH);
                beta = max( 1-(maxx - minn) ,   filter.temper_min_beta);
                g_Prob_DEATH = Prob_DEATH.^(beta);
                g_Prob_DEATH = g_Prob_DEATH ./ sum(g_Prob_DEATH);
                g_dprob(i,:) = g_Prob_DEATH;
                TEMP_hold(i,s) =  (g_Prob_DEATH(1)*N_proposal_conditional{s}(i,currsoln(i,s)+1))/(Prob_DEATH(1)*N_actual_conditional{s}(i,currsoln(i,s)+1));
            else
                TEMP_hold(i,s) =  (N_proposal_conditional{s}(i,currsoln(i,s)+1))/(N_actual_conditional{s}(i,currsoln(i,s)+1));
            end
        end
    end
    unnorm_sel_prob(i) = prod(TEMP_hold(i,:));
    
end

for sol = 2: num_samples 

    flag =1; % 1 is liu's modi % 0 is normal
%     dprob = exp(dprob_temp);% i think is the correct one
%     dprob = exp(dcost); 
%     dprob = exp(-neglogdprob_tempd); % fixed it (original when BT sent me)
    
        
    sel_prob = unnorm_sel_prob ./ sum(unnorm_sel_prob);
    cdf = cumsum(sel_prob); % equiprobable selection (uniform distribution) (ones(1,n1)
    var = sum(cdf < rand*cdf(end)) + 1;
    dprob = g_dprob(var,2);
%     for var = 1: num_tracks
%         cdf= cumsum([dprob(var),1-dprob(var)]); %Bina == 1 is death 
        cdf= cumsum([dprob,1-dprob]); %Bina == 1 is death 
        Bina = sum(cdf < rand*cdf(end)) + 1; %Bina == 2 is alive 
        if Bina == 2
            for s = 1 : num_sensors
               
%                 P0 = assign_costs{s};
%                 expnegP0 = exp(-P0);
%                 tempsamp = expnegP0(var,:);
%                 indx_temp = currsoln([1:var-1,var+1:end],s); % rule out solution in question
%                 indx_temp = indx_temp(indx_temp>0) + 1; % rule out missed d term and correct the index for 'tempsamp'
%                 tempsamp(indx_temp) = 0; 
                tempsamp =  N_proposal_conditional{s}(var,:)    ;
                if currsoln(var,s) == -1 
                    flag = 0;
                end
                if flag 
                    tempsamp_prob = tempsamp ./ sum(tempsamp);
                    unnormalized_proposal = tempsamp_prob;
                    unnormalized_proposal(currsoln(var,s)+1) =  0 ;
                    
                    denom = 1 - tempsamp_prob(currsoln(var,s)+1);
                    
                    tempsamp_prob2 = unnormalized_proposal ./ denom;
                    % % % % % %        This is Liu's modification
                    cdf = cumsum(tempsamp_prob2);
                    %                     if cdf(end) ~= 1 , error('error'); end
                    currsoln_hold  = sum(cdf < rand*cdf(end)) + 1;
                    accep_ratio = (denom)/(1-tempsamp_prob( currsoln_hold ) );
                    if rand <= min(1,accep_ratio)
                        currsoln(var,s) = currsoln_hold-1;
%                         indi_check_same_sol(s) = false; 
                    else
                        currsoln(var,s) = currsoln(var,s);
%                         indi_check_same_sol(s) = true; 
                    end
                else 
                    cdf= cumsum(tempsamp);
                    currsoln(var,s) = sum(cdf < rand*cdf(end)) + 1;
                    currsoln(var,s) = currsoln(var,s) - 1;
                end
                 % for selection prob only dprob no change
                if s == 1
                    QQ = exp(dprob_temp(var));
                    PP = 1-QQ;
                    Prob_DEATH = [PP QQ];
                    maxx =  max(  Prob_DEATH );minn = min( Prob_DEATH);
                    beta = max( 1-(maxx - minn) ,   filter.temper_min_beta);
                    g_Prob_DEATH = Prob_DEATH.^(beta);
                    g_Prob_DEATH = g_Prob_DEATH ./ sum(g_Prob_DEATH);
                    g_dprob(var,:) = g_Prob_DEATH;
                    TEMP_hold(var,s) =  (g_Prob_DEATH(1)*N_proposal_conditional{s}(var,currsoln(var,s)+1))/(Prob_DEATH(1)*N_actual_conditional{s}(var,currsoln(var,s)+1));
                else
                    TEMP_hold(var,s) =  (N_proposal_conditional{s}(var,currsoln(var,s)+1))/(N_actual_conditional{s}(var,currsoln(var,s)+1));
                end
                
            end
        else 
            currsoln(var,:) = -1*ones(1,num_sensors);
%             dprob_temp(var) = -neglogdprob_tempd(var);
%             if s == 1
    % for selection prob only dprob no change
                QQ = exp(dprob_temp(var));
                PP = 1-QQ;
                Prob_DEATH = [PP QQ];
                maxx =  max(  Prob_DEATH );minn = min( Prob_DEATH);
                beta = max( 1-(maxx - minn) ,   filter.temper_min_beta);
                g_Prob_DEATH = Prob_DEATH.^(beta);
                g_Prob_DEATH = g_Prob_DEATH ./ sum(g_Prob_DEATH);
                g_dprob(var,:) = g_Prob_DEATH;
                TEMP_hold(var,1) = (g_Prob_DEATH(2)/Prob_DEATH(2));
%             else
                TEMP_hold(var,2:end) = ones(1,length(TEMP_hold(var,:))-1 );
%             end
        end
%     end
%     dprob = exp(dprob_temp );
    unnorm_sel_prob(var) = prod(TEMP_hold(var,:));
    assignments(sol,:,:) = currsoln;
    
    
    
    for i = [1:var-1,var+1:num_tracks]
%         indx_tempp =[];
        for s = 1: num_sensors
            P0 = assign_costs{s};
            % to update N_proposal_conditional , N_actual_conditional and actual_conditional
% %             %                 actual_conditional{s}(var,:) = exp(-P0(var,:));
% %             %                 actual_conditional{s}(i,:) = exp(-P0(i,:));
            if assignments(sol-1,var,s)>0
                chk1 = actual_conditional{s}(i,assignments(sol-1,var,s)+1) == exp(-P0(i,assignments(sol-1,var,s)+1));
                actual_conditional{s}(i,assignments(sol-1,var,s)+1) = exp(-P0(i,assignments(sol-1,var,s)+1));
                 
            end
            if currsoln(var,s) > 0
                chk2 = actual_conditional{s}(i,currsoln(var,s)+1) == 0; 
                actual_conditional{s}(i,currsoln(var,s)+1) = 0;
                
            end
            if ~(assignments(sol-1,var,s)>0 ) && ~(currsoln(var,s) > 0)
                continue
            else
                if ~(assignments(sol-1,var,s)>0 )
                   chk1 = 1;
                end
                if ~(currsoln(var,s) > 0)
                   chk2 = 1;
                end
                if (chk1 && chk2)   %|| all(currsoln(:,s)' == assignments(sol-1,:,s))
                    continue;
                end
            end
       
    
            N_actual_conditional{s}(i,:) =  actual_conditional{s}(i,:) ./ sum(actual_conditional{s}(i,:));
            temp_hold =  N_actual_conditional{s}(i,:);
            temp_hold(temp_hold == 0) = inf;
            maxx =  max(  temp_hold );minn = min( temp_hold);
            beta = max( 1-(maxx - minn) , filter.temper_min_beta);
            N_proposal_conditional{s}(i,:) = N_actual_conditional{s}(i,:).^(beta); % 1-(pd(n))^2
            
            N_proposal_conditional{s}(i,:) = N_proposal_conditional{s}(i,:) ./ sum(N_proposal_conditional{s}(i,:));
            % equation 5
            %                 N_proposal_conditional{s}(var,:)=   0.5*N_proposal_conditional{s}(var,:) +   0.5*N_actual_conditional{s}(var,:) ;
                
             % for dprob
%             for_death_cost_log_P0{s} = -P0;
%            tempsamp_log_p0 = for_death_cost_log_P0{s}(i,:);
            if  assignments(sol-1,var,s)>0
                for_death_cost_log_P0{s}(i,assignments(sol-1,var,s)+1) = -P0(i,assignments(sol-1,var,s)+1);
            end
            if currsoln(var,s) > 0 
                for_death_cost_log_P0{s}(i,currsoln(var,s)+1) = log(0);
            end
            
            store_for_deathcost(i,s) = log(sum(exp(for_death_cost_log_P0{s}(i,:))));
%             logsumexp(for_death_cost_log_P0{s}(i,:),[],2);

%         indx_tempp{s} = currsoln([1:i-1,i+1:end],s); % rule out solution in question
        
%         indx_tempp{s} = indx_tempp{s}(indx_tempp{s}>0) + 1; % rule out missed d term and correct the index for 'tempsamp'
%         tempsamp_log_p0(indx_tempp{s}) = log(0);
%         store_for_deathcost(s) = logsumexp(tempsamp_log_p0,[],2);


        end
        scost = sum(store_for_deathcost(i,:));
        dprob_temp(i) = dcost(i) - log(exp(dcost(i)) + exp(scost));

%         logsumexp(dcost(i),scost);
    
%         for s = 1: num_sensors
% %             P0 = assign_costs{s};
%             % for dprob
% %             log_P0 = -P0;
%             tempsamp_log_p0 = log_P0(var,:);
%             %             indx_tempp{s} = currsoln([1:var-1,var+1:end],s); % rule out solution in question
%             indx_tempp{s} = indx_tempp{s}(indx_tempp{s}>0) + 1; % rule out missed d term and correct the index for 'tempsamp'
%             tempsamp_log_p0(indx_tempp{s}) = log(0);
%             store_for_deathcost(s) = logsumexp(tempsamp_log_p0,[],2);
%         end
        %             % for selection prob
        %             actual_conditional{s}(var,:) = [exp(dcost(var))  exp(-P0(var,:))];
        %             actual_conditional{s}(var,(actual_conditional{s}(var,:) <= eps )) = 0;
        %             actual_conditional{s}(var,indx_temp+1)= 0;
        %
        %             N_actual_conditional{s}(var,:) =  actual_conditional{s}(var,:) ./ sum(actual_conditional{s}(var,:));
        % %         end
        %             temp_hold =  N_actual_conditional{s}(var,:);
        %             temp_hold(temp_hold == 0) = inf;
        %             maxx =  max(  N_actual_conditional{s}(var,:) );minn = min( temp_hold);
        %             beta = max( 1-(maxx - minn) , 0.3);
        %             N_proposal_conditional{s}(var,:) = N_actual_conditional{s}(var,:).^(beta); % 1-(pd(n))^2
        %             %         proposal_conditional(i,:) = actual_conditional(i,:).^(1-pd(i));
        %
        %             N_proposal_conditional{s}(var,:) = N_proposal_conditional{s}(var,:) ./ sum(N_proposal_conditional{s}(var,:));
        %             % equation 5
        %             N_proposal_conditional{s}(var,:)=   0.5*N_proposal_conditional{s}(var,:) +   0.5*N_actual_conditional{s}(var,:) ;
        %             %         sel_prob(i) = proposal_conditional(i,currsoln(i))/actual_conditional(i,currsoln(i));
        %             TEMP_hold(var,s) =  (N_proposal_conditional{s}(var,currsoln(var,s)+2)/N_actual_conditional{s}(var,currsoln(var,s)+2));
        %         end
%         scost = sum(store_for_deathcost);
%         dprob_temp(var) = dcost(var) - logsumexp(dcost(var),scost);
        % for selection prob
        for s = 1 : num_sensors
            if currsoln(i,s) == -1
                if s == 1
                    QQ = exp(dprob_temp(i));
                    PP = 1-QQ;
                    Prob_DEATH = [PP QQ];
                    maxx =  max(  Prob_DEATH );minn = min( Prob_DEATH);
                    beta = max( 1-(maxx - minn) ,   filter.temper_min_beta);
                    g_Prob_DEATH = Prob_DEATH.^(beta);
                    g_Prob_DEATH = g_Prob_DEATH ./ sum(g_Prob_DEATH);
                    g_dprob(i,:) = g_Prob_DEATH;
                    TEMP_hold(i,s) = (g_Prob_DEATH(2)/Prob_DEATH(2));
                else
                    TEMP_hold(i,s) = 1;
                end
            else
                if s == 1
                    QQ = exp(dprob_temp(i));
                    PP = 1-QQ;
                    Prob_DEATH = [PP QQ];
                    maxx =  max(  Prob_DEATH );minn = min( Prob_DEATH);
                    beta = max( 1-(maxx - minn) ,   filter.temper_min_beta);
                    g_Prob_DEATH = Prob_DEATH.^(beta);
                    g_Prob_DEATH = g_Prob_DEATH ./ sum(g_Prob_DEATH);
                    g_dprob(i,:) = g_Prob_DEATH;
                    TEMP_hold(i,s) =  (g_Prob_DEATH(1)*N_proposal_conditional{s}(i,currsoln(i,s)+1))/(Prob_DEATH(1)*N_actual_conditional{s}(i,currsoln(i,s)+1));
                else
                    TEMP_hold(i,s) =  (N_proposal_conditional{s}(i,currsoln(i,s)+1))/(N_actual_conditional{s}(i,currsoln(i,s)+1));
                end
            end
        end
        unnorm_sel_prob(i) = prod(TEMP_hold(i,:));
    end
    
end
% get_assignment = permute(assignments,[2,3,1]);
make_it_a_vector = reshape(assignments,[num_samples,num_tracks*num_sensors]);
[~,I,~]= unique(make_it_a_vector,'rows');
uasses = assignments(I,:,:); 
% uasses = double(uasses+1);
%  permute(uasses,[2,3,1]);



