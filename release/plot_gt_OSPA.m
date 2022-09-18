function  plot_gt_OSPA(model,truth,meas,est,estimate_handle,performance_eval)


% name_title = 'lalalala'; 
[labelcount hold_label] = countestlabels(meas,truth);
global colorarray 
colorarray= makecolorarray(labelcount);
truth.total_tracks= labelcount;
truth.track_list= cell(meas.end,1);
for k =meas.start : meas.end 
    for eidx=1:size(truth.X{k},2)
        [indx , colorarray] = assigncolor(colorarray, truth.L{k}(:,eidx));
        truth.track_list{k} = [truth.track_list{k} indx];
    end
end
[X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks,meas);
% % %
% % %plot ground truths
% % limit= [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
% % figure;  gtvsest = gca  ;hold on;
% % for i=1:truth.total_tracks
% %     Pt= X_track(:,k_birth(i):1:k_death(i),i); Pt=Pt([1 2 3],:);
% %     plot3( Pt(1,:),Pt(2,:),Pt(3,:),'k-'); 
% %     plot3( Pt(1,1), Pt(2,1),Pt(3,:), 'ko','MarkerSize',6);
% %     plot3( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)),Pt(3,(k_death(i)-k_birth(i)+1)), 'k^','MarkerSize',6); hold on; 
% % end
% %  box on; grid on;   axis([0 8 0 4 0 2.5]);
% %  title('Ground Truths'); hold on; 
% %  for t=1:size(estimate_handle.Y_track,3)
% %     hline2= line(estimate_handle.l_birth(t):1:estimate_handle.l_death(t),estimate_handle.Y_track(1,estimate_handle.l_birth(t):1:estimate_handle.l_death(t),t),...
% %         'LineStyle','none','Marker','.','Markersize',8,'Color',estimate_handle.colorarray_est.rgb(t,:));
% %     plot3(estimate_handle.Y_track(1,meas.start:meas.end,t),estimate_handle.Y_track(3,meas.start:meas.end,t),estimate_handle.Y_track(5,meas.start:meas.end,t), ...
% %         'Marker','.','Markersize',15,'Color',estimate_handle.colorarray_est.rgb(t,:)); hold on; 
% %     pause 
% %  end
% % ylabel('y-coordinate (m)');
% % xlabel('x-coordinate (m)');
% % zlabel('z-coordinate (m)');
% %  gtvsest.FontWeight ='bold';
% % gtvsest.FontSize = 25;
% %  
% %   convert_bbs_to_world_gp(meas.Z_bbs{k});

%%
%
% %plot tracks and measurments in x/y
% figure; clf;  tracking= gca; hold on;
% 
% subplot(3,1,1) 
% box on; 
% %plot x track
% for i=1:truth.total_tracks
%     Px= X_track(:,k_birth(i):1:k_death(i),i); Px=Px([1 2 3],:);
%     hline1= line(k_birth(i):1:k_death(i),Px(1,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3)); hold on;
%     jonah1 =  plot(  k_birth(i),Px(1,1), 'ko','MarkerSize',10); hold on;
%     jonah2 = plot( k_death(i),Px(1,end), 'k^','MarkerSize',10); hold on; 
% 
% end
% 
% %plot x estimate
% for t=1:size(estimate_handle.Y_track,3)
% %     hline2= line(estimate_handle.l_birth(t):1:estimate_handle.l_death(t),estimate_handle.Y_track(1,estimate_handle.l_birth(t):1:estimate_handle.l_death(t),t),...
% %         'LineStyle','none','Marker','.','Markersize',8,'Color',estimate_handle.colorarray_est.rgb(t,:));
%         hline2= line(meas.start:meas.end,estimate_handle.Y_track(1,meas.start:meas.end,t),'LineStyle','none','Marker','.','Markersize',20,'Color',estimate_handle.colorarray_est.rgb(t,:));
% end
% 
% %plot y measurement
% subplot(3,1,2); box on;
% 
% %plot y track
% for i=1:truth.total_tracks
%         Py= X_track(:,k_birth(i):1:k_death(i),i);  Py=Py([1 2 3],:);
%         yhline1= line(k_birth(i):1:k_death(i),Py(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3)); hold on;
%     plot(  k_birth(i),Py(2,1), 'ko','MarkerSize',10); hold on;
%     plot( k_death(i),Py(2,end), 'k^','MarkerSize',10); hold on; 
% end
% 
% %plot y estimate
% for t=1:size(estimate_handle.Y_track,3)
% %     hline2= line(estimate_handle.l_birth(t):1:estimate_handle.l_death(t),estimate_handle.Y_track(2,estimate_handle.l_birth(t):1:estimate_handle.l_death(t),t),...
% %         'LineStyle','none','Marker','.','Markersize',8,'Color',estimate_handle.colorarray_est.rgb(t,:));
%     hline2= line(meas.start:meas.end,estimate_handle.Y_track(3,meas.start:meas.end,t),'LineStyle','none','Marker','.','Markersize',20,'Color',estimate_handle.colorarray_est.rgb(t,:));
% end
% %plot y measurement
% subplot(3,1,3); box on;
% 
% %plot z track
% for i=1:truth.total_tracks
%         Pz=X_track(:,k_birth(i):1:k_death(i),i);  Pz=Pz([1 2 3],:);
%         zhline1= line(k_birth(i):1:k_death(i),Pz(3,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3)); hold on;
%     plot(  k_birth(i),Pz(3,1), 'ko','MarkerSize',10); hold on;
%     plot( k_death(i),Pz(3,end), 'k^','MarkerSize',10); hold on; 
% end
% 
% % plot z estimate
% for t=1:size(estimate_handle.Y_track,3)
% %     hline2= line(estimate_handle.l_birth(t):1:estimate_handle.l_death(t),estimate_handle.Y_track(2,estimate_handle.l_birth(t):1:estimate_handle.l_death(t),t),...
% %         'LineStyle','none','Marker','.','Markersize',8,'Color',estimate_handle.colorarray_est.rgb(t,:));
%     hline2= line(meas.start:meas.end,estimate_handle.Y_track(5,meas.start:meas.end,t),'LineStyle','none','Marker','.','Markersize',20,'Color',estimate_handle.colorarray_est.rgb(t,:));
% end
% 
% subplot(3,1,1);  ylabel('x-coordinate (m)');
% set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[model.XMAX(1)-0.4 model.XMAX(2)]);
% legend([jonah1 jonah2],'Start','End');tracking= gca; tracking.FontWeight ='bold';
% tracking.FontSize = 25;
% 
% subplot(3,1,2);  ylabel('y-coordinate (m)');tracking= gca; tracking.FontWeight ='bold';
% tracking.FontSize = 25;
% set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',model.YMAX);
% subplot(3,1,3); xlabel('Frame'); ylabel('z-coordinate (m)');tracking= gca; tracking.FontWeight ='bold';
% tracking.FontSize = 25;
% set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 1.5]);
% 


%%
% %plot error - OSPA
% ospa_vals= zeros(meas.end,3);
% ospa_c= 1;
% ospa_p= 1;
% for k=meas.start:meas.end
%     if isempty(truth.X{k}) && isempty(est.X{k})
%                      [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist([],[],ospa_c,ospa_p);
%     elseif isempty(truth.X{k})
%              [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist([],est.X{k}([1 3 5],:),ospa_c,ospa_p);
%     elseif isempty(est.X{k})
%              [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(truth.X{k}([1 2 3],:),[],ospa_c,ospa_p);
%     else
%              [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(truth.X{k}([1 2 3],:),est.X{k}([1 3 5],:),ospa_c,ospa_p);
% 
% 
%     end
% end
% 
% figure; ospa= gcf; hold on;
% subplot(3,1,1); plot(meas.start:meas.end,ospa_vals(meas.start:meas.end,1),'k','linewidth',2); grid on; set(gca, 'XLim',[meas.start meas.end]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist'); ospa = gca ; ospa.FontWeight = 'bold'; ospa.FontSize = 25;
% 
% subplot(3,1,2); plot(meas.start:meas.end,ospa_vals(meas.start:meas.end,2),'k','linewidth',2); grid on; set(gca, 'XLim',[meas.start meas.end]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc'); ospa = gca ; ospa.FontWeight = 'bold'; ospa.FontSize = 25;
% 
% subplot(3,1,3); plot(meas.start:meas.end,ospa_vals(meas.start:meas.end,3),'k','linewidth',2); grid on; set(gca, 'XLim',[meas.start meas.end]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
% xlabel('Frame');
%  ospa = gca ; ospa.FontWeight = 'bold'; ospa.FontSize = 25;

%%
%plot error - OSPA^(2)
% order = 1;
% cutoff = 1;
% win_len= 10;
% 
% pick_max_width = max(exp(estimate_handle.Y_track([7 8],:,:)));
% % pick_max_width = max(exp(estimate_handle.Y_track([7 8],:,:)));
% if all(all(pick_max_width<0)), error(''); end
% 
% est_z_axis = [ ((estimate_handle.Y_track([5],:,:)) - exp(estimate_handle.Y_track([9],:,:)));((estimate_handle.Y_track([5],:,:)) + exp(estimate_handle.Y_track([9],:,:))) ];
% gt_z_axis = [(X_track(3,:,:) -X_track(6,:,:));(X_track(3,:,:) +X_track(6,:,:))];
% 
% % giou_estX = [estimate_handle.Y_track([1 3],:,:) - estimate_handle.Y_track([7 8],:,:); ...
% %     estimate_handle.Y_track([1 3],:,:) + estimate_handle.Y_track([7 8],:,:);...
% %     est_z_axis] ;
% % 
% giou_estX = [estimate_handle.Y_track([1 3],:,:) - pick_max_width; ...
%     estimate_handle.Y_track([1 3],:,:) + pick_max_width;...
%     est_z_axis] ;
% 
% giou_gtX = [X_track([1:2],:,:) -   X_track(4,:,:);X_track([1:2],:,:) + X_track(4,:,:);gt_z_axis];%X_track(4,:,:)
% 
% ospa2_cell = cell(1,length(win_len));
% for i = 1:length(win_len)
%     [  ospa2_cell{i} , ~, ~]  = giou_ospa2(giou_gtX,giou_estX,cutoff,order,win_len);
% 
% %     ospa2_cell{i} = compute_ospa2(X_track([1 2 3],meas.start:meas.end,:),estimate_handle.Y_track([1 3 5],meas.start:meas.end,:),cutoff,order,win_len);
% end
% figure('name',name_title);
order = 1;
cutoff = 1;
win_len= 10;

pick_max_width = max(exp(estimate_handle.Y_track([7 8],:,:)));
% pick_max_width = max(exp(estimate_handle.Y_track([7 8],:,:)));
if all(all(pick_max_width<0)), error(''); end

est_z_axis = [ ((estimate_handle.Y_track([5],:,:)) - exp(estimate_handle.Y_track([9],:,:)));((estimate_handle.Y_track([5],:,:)) + exp(estimate_handle.Y_track([9],:,:))) ];
gt_z_axis = [(X_track(3,:,:) -X_track(6,:,:));(X_track(3,:,:) +X_track(6,:,:))];

% giou_estX = [estimate_handle.Y_track([1 3],:,:) - estimate_handle.Y_track([7 8],:,:); ...
%     estimate_handle.Y_track([1 3],:,:) + estimate_handle.Y_track([7 8],:,:);...
%     est_z_axis] ;
% 
giou_estX = [estimate_handle.Y_track([1 3],:,:) - pick_max_width; ...
    estimate_handle.Y_track([1 3],:,:) + pick_max_width;...
    est_z_axis] ;

giou_gtX = [X_track([1:2],:,:) -   X_track(4,:,:);X_track([1:2],:,:) + X_track(4,:,:);gt_z_axis];%X_track(4,:,:)

ospa2_cell = cell(1,length(win_len));
for i = 1:length(win_len)
    if performance_eval == "3D_pos_with_extent"
        [  ospa2_cell{i} , ~, ~]  = giou_ospa2(giou_gtX,giou_estX,cutoff,order,win_len);
    elseif performance_eval == "3D_pos"
        ospa2_cell{i} = compute_ospa2(X_track([1 2 3],meas.start:meas.end,:),estimate_handle.Y_track([1 3 5],meas.start:meas.end,:),cutoff,order,win_len);
    end
end

 

% figure; ospa2= gcf; hold on;
hold on;
windowlengthlabels = cell(1,length(win_len));
subplot(3,1,1);
for i = 1:length(win_len)
    plot(meas.start:meas.end,ospa2_cell{i}(1,:),'r','linewidth',2); grid on; set(gca, 'XLim',[meas.start meas.end]);  ylabel('OSPA$^{(2)}$-gIoU','interpreter','latex');set(gca, 'YLim',[0 cutoff]);
    windowlengthlabels{i} = ['$L_w = ' int2str(win_len(i)) '$'];
end
% legend(windowlengthlabels,'interpreter','latex');
% legend(estimate_handle.name);
% legend('YOLOv3+MVGLMB-OC*');
% legend('Faster-RCNN(VGG16)+MVGLMB-OC*');
 set(gca, 'FontSize',25);
 set(gca, 'FontWeight','bold');
hold off;
%   ospa2 = gca; 
% ospa2.FontWeight = 'bold'; 
% ospa2.FontSize = 25;
% subplot(3,1,2);
% for i = 1:length(win_len)
%     plot(meas.start:meas.end,ospa2_cell{i}(2,:),'k','linewidth',2); grid on; set(gca, 'XLim',[meas.start meas.end]);  ylabel('OSPA$^{(2)}$ Loc','interpreter','latex');set(gca, 'YLim',[0 cutoff]);
%     windowlengthlabels{i} = ['$L_w = ' int2str(win_len(i)) '$'];
% end
% ospa2 = gca; 
% ospa2.FontWeight = 'bold'; 
% ospa2.FontSize = 25;
% subplot(3,1,3);
% for i = 1:length(win_len)
%     plot(meas.start:meas.end,ospa2_cell{i}(3,:),'k','linewidth',2); grid on; set(gca, 'XLim',[meas.start meas.end]);  ylabel('OSPA$^{(2)}$ Card','interpreter','latex');set(gca, 'YLim',[0 cutoff]);
%     windowlengthlabels{i} = ['$L_w = ' int2str(win_len(i)) '$'];
% end
% xlabel('Frame','interpreter','latex');
%   ospa2 = gca; 
% ospa2.FontWeight = 'bold'; 
% ospa2.FontSize = 25;


%% GioU OSPA 2 




 








%%
  %%
% %plot cardinality
% figure; cardinality= gcf; 
% subplot(2,1,1); box on; hold on;
% stairs(meas.start:meas.end,truth.N(meas.start:meas.end),'k','linewidth',2); 
% plot(meas.start:meas.end,est.N(meas.start:meas.end),'r.','Markersize',10);
% 
% grid on;
% legend(gca,'True','Estimated');
% set(gca, 'XLim',[meas.start meas.end]); set(gca, 'YLim',[0 max(truth.N)+2]);
% xlabel('Frame'); ylabel('Cardinality');card = gca; 
% card.FontWeight = 'bold'; card.FontSize = 25;
end

