function [X,N,mode,L,idxcmp] = glmb_extract_estimates(glmb,model)
  % Extract estimates via best cardinality, then
  % best component/hypothesis given best cardinality, then
  % best means of tracks given best component/hypothesis and cardinality

  [~,mode] = max(glmb.cdn);
  N = mode - 1;
  X = zeros(model.x_dim,N);
  L = zeros(2,N);
  mode = NaN(1,N);

  tmp = -inf*ones(1,length(glmb.n));
  tmp(glmb.n==N) = 0;
  [~,idxcmp] = max(glmb.w+tmp);
%   for n = 1:N
%     [~,idxtrk] = max(glmb.tt{glmb.I{idxcmp}(n)}.w);
%     X(:,n) = glmb.tt{glmb.I{idxcmp}(n)}.m(:,idxtrk);
%     L(:,n) = glmb.tt{glmb.I{idxcmp}(n)}.l;
%   end
  for n = 1:N
      [~,idxmode] = max(glmb.tt{glmb.I{idxcmp}(n)}.mode);
      mode(:,n) = idxmode;
      [~,idxtrk] = max(glmb.tt{glmb.I{idxcmp}(n)}.w{:,idxmode});
      X(:,n) = glmb.tt{glmb.I{idxcmp}(n)}.m{:,idxmode}(:,idxtrk);
      L(:,n) = glmb.tt{glmb.I{idxcmp}(n)}.l;
  end
  
end
