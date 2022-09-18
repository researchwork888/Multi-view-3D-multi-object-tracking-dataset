function [ospa2 , wass_ospa, haus_ospa] = giou_ospa2(X,Y,c,p,wl) 
% Compute GIoU based OSPA2 distance between two sets of tracks
% Inputs: 
%        X,Y -  3-D matrices RxTxN (R: states of 3D rectangles (two corners 
%        representation, e.g.  
%        X = (x1,y1,x2,y2,z1,z2) ,x2>x1 & y2>y1 && z2>z1), T: time , N: number of tracks)
%        c  -   cut-off parameter
%        p  -   p-parameter for the metric
% Output: scalar ospa distance between 2 sets of tracks X and Y 
% (average over time dimension)
  large_number = 1 ; 

  if (size(X,1) ~= size(Y,1)) || (size(X,2) ~= size(Y,2))
    error('Dimensions of X and Y are inconsistent');
  end

  if ~isnumeric(c) || ~isscalar(c) || (c < 0)
    error('c must be a positive scalar');
  end

  if ~isnumeric(p) || ~isscalar(p) || (p <= 0)
    error('p must be a positive scalar');
  end

  wl = floor(wl);
  if ~isnumeric(wl) || ~isscalar(wl) || (wl <= 0)
    error('Window length must be a positive integer');
  end
  

  eval_idx = 1:size(X,2);
  truncated = true ;
  win_off = (-wl+1):0;
  
  num_x = size(X,3);
  num_y = size(Y,3);
  num_step = size(X,2);
  num_eval = length(eval_idx);
%   result = zeros(1,num_eval);
  % First, for each time index compute the matrix of inter-point
  % distances and the track existence flags
  
  distances = zeros(num_x,num_y,num_step);
  x_exists = false(num_x,num_step);
  y_exists = false(num_y,num_step);
  
  for i = 1:num_step
      
    % Compute distance between every pair of points
    x = permute(X(:,i,:),[1 3 2]);
    y = permute(Y(:,i,:),[1 3 2]);
%     x =  x(~isnan(x));
%     y = y(~isnan(y));
    % Re-calculate the distance using the IoU as the distance
    n = size(x,2) ; 
    m = size(y,2) ; 
    xx = repmat(x,[1,m]) ;
    yy= reshape(repmat(y,[n 1]),[size(y,1) n*m]);
    ax = prod(xx(3:4,:)-xx(1:2,:)); % The rectangle areas in X
    ay = prod(yy(3:4,:)-yy(1:2,:)); % The rectangle areas in Y
%     score_XX = repmat(score_X,[1 m]) ;
    VX = ax .* (xx(6,:)-xx(5,:) );
    VY = ay .* (yy(6,:)-yy(5,:)) ; % as detection score = 1
    
    xym = min(xx,yy);
    xyM = max(xx,yy);
    Int = zeros(1,size(xx,2));
    V_Int = zeros(1,size(xx,2));

    ind = all(xyM([1:2,5],:)<xym([3:4,6],:));
%     indz = (xyM(5,:)<xym(6,:));
    Int(1,ind) = prod(xym(3:4,ind)-xyM(1:2,ind));
    
    score_height(1,ind) =xym(6,ind)-xyM(5,ind);%xyM(5,ind) - abs( xym(5,ind)-xyM(5,ind));
    
    V_Int(1,ind) = Int(1,ind) .* score_height(1,ind); 

    Unn = ax+ay-Int;
    IoU = Int./Unn;
    Cc = prod(xyM(3:4,:)-xym(1:2,:));
    V_Cc = Cc .*  (xyM(6,:)-xym(5,:)); % fix thix

    GIoU = IoU - ((Cc-Unn)./Cc);
    
    V_Unn = VX + VY - V_Int ;
    V_IoU = V_Int./V_Unn ;
    V_GIoU = V_IoU - ((V_Cc-V_Unn)./V_Cc) ;
    
%     d = reshape(-log10((GIoU+1)/2) ,[n m]);
%     d = reshape(0.5*(1-GIoU) ,[n m]);
    d = reshape(0.5*(1-V_GIoU) ,[n m]);

%     D = min(c,d).^p; 
%     d = nan(length(x) , length(y)) ;     
%     for xidx = 1 : length(x)
%         for yidx = 1 : length(y)
%             if ~isnan(x(1,xidx)) && ~isnan(y(1,yidx))
%                 XYm = min( x(:,xidx) , y(:,yidx) ) ; 
%                 XYM = max( x(:,xidx) , y(:,yidx) ) ;
%                 Int = 0 ; 
%                 if XYM(1:2,:)<XYm(3:4,:) 
%                     Int = prod(XYm(3:4)-XYM(1:2));
%                 end
%                 AX = prod(x([3,4],xidx) - x([1,2],xidx)) ; 
%                 AY = prod(y([3,4],yidx) - y([1,2],yidx)) ;
%                 Unn = AX + AY - Int ; 
%                 IoU = Int/Unn ; 
%                 Cc = prod(XYM(3:4,:)-XYm(1:2,:));
%                 GIoU = IoU - ((Cc-Unn)./Cc);
%                 d(xidx , yidx) = -log10((GIoU+1)/2) ; 
%             end
%         end
%     end
    
    
    % Compute track existence flags
    x_exists(:,i) = ~isnan(x(1,:));
    y_exists(:,i) = ~isnan(y(1,:));
        
    % Distance between an empty and non-empty state
    one_exists = bsxfun(@xor,x_exists(:,i),y_exists(:,i)');
    d(one_exists) = c^p;     
    
    % Find times when neither object exists
    neither_exists = bsxfun(@and,~x_exists(:,i),~y_exists(:,i)');
    if truncated
      % Truncated window, exclude times when neither objects exists
      d(neither_exists) = NaN;
    else
      % Full window, distance between empty states is zero
      d(neither_exists) = 0;
    end
    
    % Store the distance matrix for this step
    distances(:,:,i) = d;
    
  end
  
  % Cap all inter-point distances at c^p
  if truncated
    % Truncated window
    distances = min(c^p,distances,'includenan');
  else
    % Full window
    distances = min(c^p,distances);
  end
  
  % Compute the OSPA(2) at each evaluation point
    for i = 1:num_eval

%         i = num_eval ;
        
        % Window indices
        win_idx = eval_idx(i) + win_off;
        idx_val = (win_idx > 0) & (win_idx <= num_step);
        win_idx = win_idx(idx_val);
        
        % Compute the matrix of weighted time-averaged
        % OSPA distances between tracks
        trk_dist = mean(distances(:,:,win_idx),3,'omitnan');
        trk_dist(isnan(trk_dist)) = 0;
        
        % Get the number of objects in X and Y that exist
        % at any time inside the current window
        valid_rows = any(x_exists(:,win_idx),2);
        valid_cols = any(y_exists(:,win_idx),2);
        m = sum(valid_rows);
        n = sum(valid_cols);
        
        % Solve the optimal assignment problem
        trk_dist = trk_dist(valid_rows,valid_cols);
        if isempty(trk_dist)
            cost = 0;
            wass_cost =0;
            N_wass=0;
             vx = max(min(trk_dist,[],2));
            vy = max(min(trk_dist,[],1));
        else
            if m > n
                trk_dist = trk_dist';
            end
            s_n = size(trk_dist,1) ;
            s_m = size(trk_dist,2) ;
            [~,cost] = lapjv(trk_dist);
            vx = max(min(trk_dist,[],2));
            vy = max(min(trk_dist,[],1));
            if c>0.99999999
                if s_m == s_n
                    wass_cost = cost ;
                    N_wass = s_m ;
                else
                    %%% Solve wass assignment using linear programming
                    N = s_n*s_m ;
                    f = reshape(trk_dist',[N,1]) ;
                    
                    A1 = zeros(s_n,N) ; A2 = zeros(s_m,N) ;
                    for ii = 1:s_n
                        for j = 1:s_m
                            k = j + (ii - 1) * s_m;
                            A1(ii, k) = 1;
                            A2(j, k) = 1;
                        end
                    end
                    A = [A1; A2];
                    W1 = (1/s_n)*ones(s_n,1) ; W2 = (1/s_m)*ones(s_m,1) ;
                    b = [W1 ; W2];
                    Aeq = ones(s_m + s_n, s_m * s_n);
                    beq = ones(s_m + s_n, 1) * min(sum(W1), sum(W2));
                    
                    % lower bound
                    lb = zeros(1, s_m * s_n);
                    ub = ones(1, s_m * s_n);
                    options = optimoptions('linprog','Display','none');
                    [~, wass_cost] = linprog(f, A, b, Aeq, beq, lb, ub, options);
                    N_wass = 1 ;
                    %             N_wass_lin = sum(N_wass_lin) ;
                end
            else
                wass_cost = 0 ;
                N_wass = 1 ;
            end
            
            
        end
        
        % Compute the OSPA(2) distance for the current time index
        if max(m,n) == 0
            ospa2(:,i) = 0;
            wass_ospa = 0 ;
            haus_ospa = 0 ;
        else
            %       result(1,i) = ( ( c^p * abs(m-n) + cost ) / max(m,n) ) ^ (1/p);
            %       result(2,i) = ( cost / max(m,n) ) ^ (1/p);
            %       result(3,i) = ( c^p * abs(m-n) / max(m,n) ) ^ (1/p);
            ospa2(1,i) = ( ( c^p * abs(m-n) + cost ) / max(m,n) ) ^ (1/p);
            ospa2(2,i) = ( cost / max(m,n) ) ^ (1/p);
            ospa2(3,i) = ( c^p * abs(m-n) / max(m,n) ) ^ (1/p);
            wass_ospa = wass_cost/N_wass ;
            haus_ospa = max(vx,vy);
        end
    end
    
  end
  
