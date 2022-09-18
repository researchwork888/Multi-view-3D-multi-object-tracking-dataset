function pS = compute_pS(model,mode,w,XX,L,k)


if isempty(XX)
    pS= [];
else
%     if length(mode) == 1
%         [~,ind]=max(mode); % taking MAP of mode.
%         [~,indx]=max(w(:,ind)) ; % taking MAP of GM.
%         X = XX(:,indx);
%     else
        [~,ind]=max(mode); % taking MAP of mode.
        [~,indx]=max(w{:,ind}) ; % taking MAP of GM.
        X = XX{:,ind}(:,indx);
        
%     end
    xbound_min = model.XMAX(1);
    xbound_max = model.XMAX(2);
    ybound_min = model.YMAX(1);
    ybound_max = model.YMAX(2);
    control_param = 0.5; %0.6
    for i = 1 :size(X,2)
        X_hold = X(:,i);
        
        L_hold = L(:,i);
        
        if X_hold(1) > xbound_min && X_hold(1) < xbound_max && X_hold(3) > ybound_min  && X_hold(3) < ybound_max
            scene_mask = exp(model.P_S);
        else
            scene_mask = 1 - exp(model.P_S);
        end
        
        pS(i) = scene_mask / (1 + exp(-control_param*(k-(L_hold'*[1;0]))));
        pS(i) = log(pS(i));
    end
    %     pS= model.P_S*ones(size(X,2),1);
    %     pS= pS(:);
end