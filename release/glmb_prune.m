function glmb_out = glmb_prune(glmb_in,filter)
  
  % Keep components with weights above the threshold
  idxkeep = find(glmb_in.w > filter.comp_min_weight_retain);
  glmb_out.tt = glmb_in.tt;
  glmb_out.w = glmb_in.w(idxkeep);
  glmb_out.I = glmb_in.I(idxkeep);
  glmb_out.n = glmb_in.n(idxkeep);

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
  
end
