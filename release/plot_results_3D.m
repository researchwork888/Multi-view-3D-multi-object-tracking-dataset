function  [handle,estimate_handle]=  plot_results_3D(model,meas,est,do_print_videos)

[labelcount,hold_label]= countestlabels();
colorarray= makecolorarray(labelcount);
est.total_tracks= labelcount;
est.track_list= cell(meas.end,1);
for k =meas.start : meas.end
    for eidx=1:size(est.X{k},2)
        est.track_list{k} = [est.track_list{k} assigncolor(est.L{k}(:,eidx))];
    end
end
% [Y_track,l_birth,l_death]= extract_tracks(est.X,est.mode,est.track_list,est.total_tracks,meas);
[Y_track,l_birth,l_death]= extract_tracks(est.X,est.track_list,est.total_tracks,meas);

estimate_handle.Y_track = Y_track;
estimate_handle.l_birth = l_birth;
estimate_handle.l_death = l_death;
estimate_handle.colorarray_est = colorarray;
if do_print_videos == true
    %% for 3D plot
    rectan = [model.XMAX(1) model.YMAX(1);
        model.XMAX(2) model.YMAX(1);
        model.XMAX(2) model.YMAX(2);
        model.XMAX(1) model.YMAX(2);
        model.XMAX(1) model.YMAX(1)];
    rundidx = 1;
    for i = meas.start :meas.end
        figure('visible','off');  threedplot = gcf; threedplot_gca =gca;
        for t = 1:size(Y_track,3)
            if ~any(isnan(Y_track(:,i,t)))
                
                temp =  Y_track([1 3 5 7 8 9],i,t);
                temp([4 5 6],:) = exp(temp([4 5 6],:));
                %             [x, y, z]  = ellipsoid(temp(1),temp(2),temp(3)/2,temp(4),temp(5),temp(3)/2,(model.ellipsoid_n));
                [x, y, z]  = ellipsoid(temp(1),temp(2),temp(3),temp(4),temp(5),temp(6),(model.ellipsoid_n));
                
                vet = [];
                facs = [];
                for asd = 1 : (model.ellipsoid_n + 1 )
                    tempp = [x(:,asd)';y(:,asd)';z(:,asd)'];
                    %             plot3(tempp(1,:),tempp(2,:),tempp(3,:))
                    if asd == 2
                        facs = [1:(model.ellipsoid_n + 1 ), flip((model.ellipsoid_n + 1 + 1):(asd*(model.ellipsoid_n + 1 ))) ];
                    end
                    if  asd > 2
                        facs = [facs; [(model.ellipsoid_n + 1 )*(asd-1)-((model.ellipsoid_n + 1 )*(asd-1))/(asd-1)+1:(model.ellipsoid_n + 1 )*(asd-1), flip(((model.ellipsoid_n + 1 )*(asd-1)+1):((model.ellipsoid_n + 1 )*asd))]];
                    end
                    vet = [vet tempp];
                end
                patch('Faces',facs,'Vertices',vet','FaceColor',colorarray.rgb(t,:),'FaceAlpha',0.2,'EdgeColor',colorarray.rgb(t,:),'LineWidth',1);
                hold on;
                grid on; box on;
                plot3(temp(1),temp(2),temp(3),'*k','MarkerSize',12); hold on;
                plot(temp(1),temp(2),'*r','MarkerSize',12); hold on;
                text(temp(1), temp(2),temp(3)+temp(6), ['ID:(' num2str(hold_label(t,1)) ',' num2str(hold_label(t,2)) ')'], 'Color','red', 'FontWeight','bold', 'FontSize',10); hold on;
                plot3(rectan(:,1),rectan(:,2),zeros(size(rectan,1),1),'-k','linewidth',2);
                axis([0 8 0 4 0 2]);
                xlabel('x-axis (m)'); ylabel('y-axis (m)');zlabel('z-axis (m)');
                view(-64,72)
                %                       legend(asdasd,'Ground Plane')
                threedplot_gca.FontWeight ='bold';
            else
                
            end
        end
        
        plot3(rectan(:,1),rectan(:,2),zeros(size(rectan,1),1),'-k','linewidth',2);
        axis([0 8 0 4 0 2]);
        grid on; box on;
        view(-64,72)
        xlabel('x-axis (m)'); ylabel('y-axis (m)');zlabel('z-axis (m)');
        
        %                       legend(asdasd,'Ground Plane')
        threedplot_gca.FontWeight ='bold';
        
        handle.threedplot(rundidx) = getframe(threedplot);
        %       pause();
        hold off;
        clf
        rundidx = rundidx + 1;
        
    end
    
    %% For each camera
    
    runindx = 1;
    alph = 0.1;
    for k = meas.start : meas.end
        % Cam 1
        figure('visible','off'); img1 = gcf; set(img1,'visible','off');
        imshow(uint8(meas.img1{k}),'Border','tight');
        text(5, 25, 'Frame', 'Color','r', 'FontWeight','bold', 'FontSize',30);
        text(5, 70, num2str(k), 'Color','r', 'FontWeight','bold', 'FontSize',30);
        hold on;axis on; xlabel('x-axis');ylabel('y-label');
        for t=1:size(Y_track,3)
            if ~any(isnan(Y_track(:,k,t)))
                
                temp =  Y_track([1 3 5 7 8 9],k,t);
                temp([4 5 6],:) = exp(temp([4 5 6],:));
                [x, y, z]  = ellipsoid(temp(1),temp(2),temp(3),temp(4),temp(5),temp(6),(model.ellipsoid_n));
                vet = [];
                facs = [];
                for asd = 1 : (model.ellipsoid_n + 1 )
                    tempp = [x(:,asd)';y(:,asd)';z(:,asd)'];
                    if asd == 2
                        facs = [1:(model.ellipsoid_n + 1 ), flip((model.ellipsoid_n + 1 + 1):(asd*(model.ellipsoid_n + 1 ))) ];
                    end
                    if  asd > 2
                        facs = [facs; [(model.ellipsoid_n + 1 )*(asd-1)-((model.ellipsoid_n + 1 )*(asd-1))/(asd-1)+1:(model.ellipsoid_n + 1 )*(asd-1), flip(((model.ellipsoid_n + 1 )*(asd-1)+1):((model.ellipsoid_n + 1 )*asd))]];
                    end
                    vet = [vet tempp];
                end
                cam_mat = eval(['model.cam',num2str(1),'_cam_mat']);
                temp = cam_mat*[vet; ones(1,size(vet,2))];
                img1_vert =  temp([1 2],:) ./ temp(3,:)   ;
                patch('Faces',facs,'Vertices',img1_vert','FaceColor',colorarray.rgb(t,:),'FaceAlpha',alph,'EdgeColor',colorarray.rgb(t,:),'LineWidth',1); hold on;
                text((min(img1_vert(1,:))+max(img1_vert(1,:)))/2, min(img1_vert(2,:)), ['ID:(' num2str(hold_label(t,1)) ',' num2str(hold_label(t,2)) ') ',model.mode_type{1}], 'Color','r', 'FontWeight','bold', 'FontSize',15); hold on;
                
                
                
            else
                
            end
        end
        handle.img1_click(runindx) = getframe(img1);
        
        
        % Cam 2
        figure('visible','off'); img2  = gcf; set(img2,'visible','off');
        imshow(uint8(meas.img2{k}),'Border','tight');
        text(5, 25, 'Frame', 'Color','r', 'FontWeight','bold', 'FontSize',30);
        text(5, 70, num2str(k), 'Color','r', 'FontWeight','bold', 'FontSize',30);
        hold on;axis on; xlabel('x-axis');ylabel('y-label');
        for t=1:size(Y_track,3)
            if ~any(isnan(Y_track(:,k,t)))
                
                temp =  Y_track([1 3 5 7 8 9],k,t);
                temp([4 5 6],:) = exp(temp([4 5 6],:));  
                [x, y, z]  = ellipsoid(temp(1),temp(2),temp(3),temp(4),temp(5),temp(6),(model.ellipsoid_n));
                vet = [];
                facs = [];
                for asd = 1 : (model.ellipsoid_n + 1 )
                    tempp = [x(:,asd)';y(:,asd)';z(:,asd)'];
                    if asd == 2
                        facs = [1:(model.ellipsoid_n + 1 ), flip((model.ellipsoid_n + 1 + 1):(asd*(model.ellipsoid_n + 1 ))) ];
                    end
                    if  asd > 2
                        facs = [facs; [(model.ellipsoid_n + 1 )*(asd-1)-((model.ellipsoid_n + 1 )*(asd-1))/(asd-1)+1:(model.ellipsoid_n + 1 )*(asd-1), flip(((model.ellipsoid_n + 1 )*(asd-1)+1):((model.ellipsoid_n + 1 )*asd))]];
                    end
                    vet = [vet tempp];
                end
                cam_mat = eval(['model.cam',num2str(2),'_cam_mat']);
                temp = cam_mat*[vet; ones(1,size(vet,2))];
                img1_vert =  temp([1 2],:) ./ temp(3,:)   ;
                
                patch('Faces',facs,'Vertices',img1_vert','FaceColor',colorarray.rgb(t,:),'FaceAlpha',alph,'EdgeColor',colorarray.rgb(t,:),'LineWidth',1); hold on;
                text((min(img1_vert(1,:))+max(img1_vert(1,:)))/2, min(img1_vert(2,:)), ['ID:(' num2str(hold_label(t,1)) ',' num2str(hold_label(t,2)) ')'], 'Color','r', 'FontWeight','bold', 'FontSize',15); hold on;
                
            else
                
            end
        end
        handle.img2_click(runindx)  = getframe(img2);
        
        % Cam 3
        figure('visible','off');img3 = gcf; set(img3,'visible','off');
        imshow(uint8(meas.img3{k}),'Border','tight');
        text(5, 25, 'Frame', 'Color','r', 'FontWeight','bold', 'FontSize',30);
        text(5, 70, num2str(k), 'Color','r', 'FontWeight','bold', 'FontSize',30);
        hold on;axis on; xlabel('x-axis');ylabel('y-label');
        for t=1:size(Y_track,3)
            if ~any(isnan(Y_track(:,k,t)))
                temp =  Y_track([1 3 5 7 8 9],k,t);
                temp([4 5 6],:) = exp(temp([4 5 6],:));
                [x, y, z]  = ellipsoid(temp(1),temp(2),temp(3),temp(4),temp(5),temp(6),(model.ellipsoid_n));
                vet = [];
                facs = [];
                for asd = 1 : (model.ellipsoid_n + 1 )
                    tempp = [x(:,asd)';y(:,asd)';z(:,asd)'];
                    if asd == 2
                        facs = [1:(model.ellipsoid_n + 1 ), flip((model.ellipsoid_n + 1 + 1):(asd*(model.ellipsoid_n + 1 ))) ];
                    end
                    if  asd > 2
                        facs = [facs; [(model.ellipsoid_n + 1 )*(asd-1)-((model.ellipsoid_n + 1 )*(asd-1))/(asd-1)+1:(model.ellipsoid_n + 1 )*(asd-1), flip(((model.ellipsoid_n + 1 )*(asd-1)+1):((model.ellipsoid_n + 1 )*asd))]];
                    end
                    vet = [vet tempp];
                end
                cam_mat = eval(['model.cam',num2str(3),'_cam_mat']);
                temp = cam_mat*[vet; ones(1,size(vet,2))];
                img1_vert =  temp([1 2],:) ./ temp(3,:)   ;
                patch('Faces',facs,'Vertices',img1_vert','FaceColor',colorarray.rgb(t,:),'FaceAlpha',alph,'EdgeColor',colorarray.rgb(t,:),'LineWidth',1); hold on;
                text((min(img1_vert(1,:))+max(img1_vert(1,:)))/2, min(img1_vert(2,:)), ['ID:(' num2str(hold_label(t,1)) ',' num2str(hold_label(t,2)) ')'], 'Color','r', 'FontWeight','bold', 'FontSize',15); hold on;
            else
                
            end
        end
        handle.img3_click(runindx) = getframe(img3);
        
        % Cam 4
        figure('visible','off');img4  = gcf; set(img4,'visible','off');
        imshow(uint8(meas.img4{k}),'Border','tight');
        text(5, 25, 'Frame', 'Color','r', 'FontWeight','bold', 'FontSize',30);
        text(5, 70, num2str(k), 'Color','r', 'FontWeight','bold', 'FontSize',30);
        hold on;axis on; xlabel('x-axis');ylabel('y-label');
        for t=1:size(Y_track,3)
            if ~any(isnan(Y_track(:,k,t)))
                temp =  Y_track([1 3 5 7 8 9],k,t);
                temp([4 5 6],:) = exp(temp([4 5 6],:)); % [1 3 5 7 9]
                [x, y, z]  = ellipsoid(temp(1),temp(2),temp(3),temp(4),temp(5),temp(6),(model.ellipsoid_n));
                vet = [];
                facs = [];
                for asd = 1 : (model.ellipsoid_n + 1 )
                    tempp = [x(:,asd)';y(:,asd)';z(:,asd)'];
                    if asd == 2
                        facs = [1:(model.ellipsoid_n + 1 ), flip((model.ellipsoid_n + 1 + 1):(asd*(model.ellipsoid_n + 1 ))) ];
                    end
                    if  asd > 2
                        facs = [facs; [(model.ellipsoid_n + 1 )*(asd-1)-((model.ellipsoid_n + 1 )*(asd-1))/(asd-1)+1:(model.ellipsoid_n + 1 )*(asd-1), flip(((model.ellipsoid_n + 1 )*(asd-1)+1):((model.ellipsoid_n + 1 )*asd))]];
                    end
                    vet = [vet tempp];
                end
                cam_mat = eval(['model.cam',num2str(4),'_cam_mat']);
                temp = cam_mat*[vet; ones(1,size(vet,2))];
                img1_vert =  temp([1 2],:) ./ temp(3,:)   ;
                patch('Faces',facs,'Vertices',img1_vert','FaceColor',colorarray.rgb(t,:),'FaceAlpha',alph,'EdgeColor',colorarray.rgb(t,:),'LineWidth',1); hold on;
                text((min(img1_vert(1,:))+max(img1_vert(1,:)))/2, min(img1_vert(2,:)), ['ID:(' num2str(hold_label(t,1)) ',' num2str(hold_label(t,2)) ')'], 'Color','r', 'FontWeight','bold', 'FontSize',15); hold on;
            else
                
            end
        end
        handle.img4_click(runindx)  = getframe(img4);
        runindx = runindx + 1 ;
        
    end
else
        handle = []; 
end



function ca= makecolorarray(nlabels)
    lower= 0.1;
    upper= 0.9;
    rrr= rand(1,nlabels)*(upper-lower)+lower;
    ggg= rand(1,nlabels)*(upper-lower)+lower;
    bbb= rand(1,nlabels)*(upper-lower)+lower;
    ca.rgb= [rrr; ggg; bbb]';
    ca.lab= cell(nlabels,1);
    ca.cnt= 0;   
end

function idx= assigncolor(label)
    str= sprintf('%i*',label);
    tmp= strcmp(str,colorarray.lab);
    if any(tmp)
        idx= find(tmp);
    else
        colorarray.cnt= colorarray.cnt + 1;
        colorarray.lab{colorarray.cnt}= str;
        idx= colorarray.cnt;
    end
end

function [count,c]= countestlabels
    labelstack= [];
    for k=meas.start:meas.end
        labelstack= [labelstack est.L{k}];
    end
    [c,~,~]= unique(labelstack','rows');
    count=size(c,1);
end

end

function [X_track,k_birth,k_death]= extract_tracks(X,track_list,total_tracks,meas)

K= size(X,1); 
x_dim= size(X{1},1); 
% x_dim = x_dim + 1 ; % for mode 
k=K-1; while x_dim==0, x_dim= size(X{k},1); k= k-1; end
X_track= NaN(x_dim,K,total_tracks);
k_birth= zeros(total_tracks,1);
k_death= zeros(total_tracks,1);

max_idx= 0;
for k=meas.start:meas.end
    if ~isempty(X{k})
        X_track(:,k,track_list{k})= [X{k}];%; mode{k}];
    end
    if max(track_list{k})> max_idx %new target born?
        idx= find(track_list{k}> max_idx);
        k_birth(track_list{k}(idx))= k;
    end
    if ~isempty(track_list{k}), max_idx= max(track_list{k}); end
    k_death(track_list{k})= k;
end
end

function Xc= get_comps(X,c)

if isempty(X)
    Xc= [];
else
    Xc= X(c,:);
end
end