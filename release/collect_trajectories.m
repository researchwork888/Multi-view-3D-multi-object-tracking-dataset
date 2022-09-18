% usage: this function for fixed birth, no spawn model
function prior_estimator = collect_trajectories(model,glmb_update,idxcmp,prior_estimator)
    % Check the current hypothesis
    num_traj =  length(glmb_update.I{idxcmp}) ;
    for tjidx = 1 : num_traj
        tabidx = glmb_update.I{idxcmp}(tjidx) ; 
        curr_label = glmb_update.tt{tabidx}.l ; 
        curr_theta = glmb_update.tt{tabidx}.ah ; 
        
        % check if label exists in the current trajectories list
        [repeat_label_flag,new_traj_index] = check_repeat_label(prior_estimator.labels , curr_label) ; 
        if repeat_label_flag
            % if this is a repeated label, update the trajectories
            prior_estimator.thetas{new_traj_index} = curr_theta ; 
        else
            birth_comp = curr_label(end) ; 
            curr_init.m = model.m_birth{birth_comp} ; 
            curr_init.P = model.P_birth{birth_comp} ; 
            curr_init.w = model.w_birth{birth_comp} ; 
            prior_estimator.labels = [prior_estimator.labels , {curr_label}] ; 
            prior_estimator.thetas = [prior_estimator.thetas , {curr_theta}] ; % put ah to a cell
            prior_estimator.inits = [prior_estimator.inits , {curr_init}] ; 
        end
    end
end



    function [repeat_label_flag,index_prior] = check_repeat_label(prior_labels , current_label) % usage vector of cell vs cell
        repeat_label_flag = false ; 
%         current_label = {current_label} ; % uncomment this line if
%         current_label is a cell 
        plidx = 1 ;
        while plidx <= length(prior_labels) && ~repeat_label_flag
            prior_label = prior_labels{plidx} ; 
            if length(prior_label) == length(current_label)
                if all ((prior_label - current_label) == 0)
                    repeat_label_flag = true ;
                end
            end
            plidx = plidx + 1 ; 
        end
        if repeat_label_flag
            index_prior = plidx - 1 ; % minus 1 as I plused 1 in the previous step 
        else
            index_prior = plidx ; % this is index of the new traejctories in the prior object
        end
    end
            
        
    