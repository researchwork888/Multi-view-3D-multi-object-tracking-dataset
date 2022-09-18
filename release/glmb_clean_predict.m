function glmb_temp = glmb_clean_predict(glmb_raw)

%   % Hash label sets, find unique ones, merge all duplicates
%   hash_values = int32(zeros(1,length(glmb_raw.w)));
%   for hidx = 1:length(glmb_raw.w)
%     hash_values(hidx) = hash_matrix(int32(sort(glmb_raw.I{hidx})));
%   end
  
  %hash label sets, find unique ones, merge all duplicates
  for hidx= 1:length(glmb_raw.w)
      glmb_raw.hash{hidx}= sprintf('%i*',sort(glmb_raw.I{hidx}(:)'));
  end
  
  [cu,~,ic] = unique(glmb_raw.hash);
  
  glmb_temp.tt = glmb_raw.tt;
  glmb_temp.w = -inf*ones(1,length(cu));
  glmb_temp.I = cell(1,length(cu));
  glmb_temp.n = zeros(1,length(cu));
  for hidx = 1:length(ic)
    glmb_temp.w(ic(hidx)) = logsumexp(glmb_temp.w(ic(hidx)),glmb_raw.w(hidx));
    glmb_temp.I{ic(hidx)} = glmb_raw.I{hidx};
    glmb_temp.n(ic(hidx)) = glmb_raw.n(hidx);
  end
  glmb_temp.cdn = glmb_raw.cdn;
  
end
