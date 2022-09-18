function glmb_out = glmb_truncate(glmb_in,filter)
  
  if length(glmb_in.w) > filter.H_upd
    
    max_card = max(glmb_in.n);
    num_components = length(glmb_in.w);
    keep = false(1,num_components);
    
    % Sort components in descending weight order
    [~,idxsort] = sort(glmb_in.w,'descend');
    
    % Enforce minimum components per cardinality
    card_num_kept = zeros(1,max_card+1);
    for i = 1:num_components
      idx = idxsort(i);
      card = glmb_in.n(idx);
      if (card_num_kept(card+1) < filter.comp_min_retain_per_card)
        keep(idx) = true;
        card_num_kept(card+1) = card_num_kept(card+1) + 1;
      end
    end
    
    % Fill the remaining components with highest weights
    num_remaining = filter.comp_max_num_retain - sum(keep);
    for i = 1:num_components
      idx = idxsort(i);
      if (num_remaining > 0) && ~keep(idx)
        keep(idx) = true;
        num_remaining = num_remaining - 1;
      end
    end
    
    % Create capped GLMB
    glmb_out.tt = glmb_in.tt;
    glmb_out.w = glmb_in.w(keep);
    glmb_out.I = glmb_in.I(keep);
    glmb_out.n = glmb_in.n(keep);
    
    % Renormalise component weights
    glmb_out.w = glmb_out.w - logsumexp(glmb_out.w,[],2);
    
    % Recompute cardinality distribution
    max_card = max(glmb_out.n);
    glmb_out.cdn = zeros(1,max_card+1);
    for card = 0:max_card
      idx = (glmb_out.n==card);
      if any(idx)
        glmb_out.cdn(card+1) = exp(logsumexp(glmb_out.w(idx),[],2));
      end
    end
    
  else
    
    % Number of components is already less than the maximum
    glmb_out = glmb_in;
    
  end

end
