function redraw(k,model,meas,est)
alph = 0.1;
im = meas.img1{k};
imshow(im,'Border','tight');
hold on;
text(5, 25, 'Frame', 'Color','r', 'FontWeight','bold', 'FontSize',30);
text(5, 70, num2str(k), 'Color','r', 'FontWeight','bold', 'FontSize',30);

% Cam 1
for t=1:size(est.Y_track,3)
    if ~any(isnan(est.Y_track(:,k,t)))
        temp =  est.Y_track([1 3 5 7 8 9],k,t);
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
        patch('Faces',facs,'Vertices',img1_vert','FaceColor',est.colorarray.rgb(t,:),'FaceAlpha',alph,'EdgeColor',est.colorarray.rgb(t,:),'LineWidth',1);
        hold on;
        text((min(img1_vert(1,:))+max(img1_vert(1,:)))/2, min(img1_vert(2,:)), ['ID:(' num2str(est.hold_label(t,1)) ',' num2str(est.hold_label(t,2)) ') ',model.mode_type{est.Y_track(model.x_dim+1,k,t)}], 'Color','r', 'FontWeight','bold', 'FontSize',15);
        hold on;
        
        
        
    else
        
    end
end
hold off;

end