function [assignments] = gibbs_multisensor_approx_cheap(neglogdcost,neglogdprob, neglogcostc,num_comp_request)

assign_costs = neglogcostc; 

num_tracks = size(neglogdprob,1);
num_sensors = size(assign_costs,1);
num_samples = num_comp_request; 

dcost =  (-neglogdcost); % fixed 

currsoln = zeros(num_tracks,num_sensors);
assignments = zeros(num_samples,num_tracks,num_sensors);
store_for_deathcost = zeros(num_sensors,1);
dprob_temp = zeros(num_tracks,1);

for sol = 2: num_samples 
    for var = 1: num_tracks 
        for s = 1: num_sensors
            P0 = assign_costs{s};
            log_P0 = -P0;
            tempsamp_log_p0 = log_P0(var,:);
            indx_temp = currsoln([1:var-1,var+1:end],s);                    % rule out solution in question
            indx_temp = indx_temp(indx_temp>0) + 1;                         % rule out missed d term and correct the index for 'tempsamp'
            tempsamp_log_p0(indx_temp) = log(0);
            store_for_deathcost(s) = log(sum(exp(tempsamp_log_p0)));
        end
        scost = sum(store_for_deathcost);
        dprob_temp(var) = dcost(var) - log(exp(dcost(var)) + exp(scost));

    end
  dprob = exp(dprob_temp); 
%     dprob = exp(dcost); 
%     dprob = exp(-neglogdprob);
    
    for var = 1: num_tracks
        cdf= cumsum([dprob(var),1-dprob(var)]);                             %Bina == 1 is death, Bina == 2 is alive
        Bina = sum(cdf < rand*cdf(end)) + 1; 
        if Bina == 2
            for s = 1 : num_sensors
                P0 = assign_costs{s};
                expnegP0 = exp(-P0);
                tempsamp = expnegP0(var,:);
                indx_temp = currsoln([1:var-1,var+1:end],s);                % rule out solution in question
                indx_temp = indx_temp(indx_temp>1) ;                        % rule out missed d term and correct the index for 'tempsamp'
                tempsamp(indx_temp) = 0; 

                cdf= cumsum(tempsamp);
                currsoln(var,s) = sum(cdf < rand*cdf(end)) + 1; 
            end
        else 
            currsoln(var,:) = -1*ones(1,num_sensors);
        end
    end
    assignments(sol,:,:) = currsoln;
    
end
make_it_a_vector = reshape(assignments,[num_samples,num_tracks*num_sensors]);
[~,I,~]= unique(make_it_a_vector,'rows');
assignments = assignments(I,:,:);


