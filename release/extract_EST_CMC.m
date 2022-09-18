
%%
if strcmp(model.dataset,'CMC2')
    fileID = fopen('./Results/EST_CMC2_WORLD_CENTROID.txt','w');
    for i = meas.start : meas.end
        store_L = [est.L{i}' ];
        [c,~,~]= unique(store_L,'rows');
        store_L = c' ;
        store_X = [];
        for dsa = 1 :size(store_L,2)
            store_X = [];
            indx1 = find(all(est.L{i} == store_L(:,dsa) ));
            %          indx2 = find(all(L2{i} == store_L(:,dsa) ));
            %          indx3 = find(all(L3{i} == store_L(:,dsa) ));
            %          indx4 = find(all(L4{i} == store_L(:,dsa) ));
            
            %             if isempty(indx1)
            %
            %             else
            %                 store_X = [store_X X1{i}(:,indx1)];
            %
            %             end
            %             if isempty(indx2)
            %             else
            %                 store_X = [store_X X2{i}(:,indx2)];
            %             end
            %             if isempty(indx3)
            %             else
            %                 store_X = [store_X X3{i}(:,indx3)];
            %             end
            %             if isempty(indx4)
            %             else
            %                 store_X = [store_X X4{i}(:,indx4)];
            %             end
            
            %         truth.X{i,:}(:,dsa) = sum(store_X,2)/size(store_X,2);
            %         truth.L{i,:}(:,dsa) = store_L(:,dsa);
            %         truth.N(i,:) =  size(store_L,2);
            pick_xy_lengths = max(exp(est.X{i,:}(7,indx1)),exp(est.X{i,:}(8,indx1)));
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,store_L(1,dsa),-1,-1,-1,-1,-1,est.X{i,:}(1,indx1),est.X{i,:}(3,indx1),est.X{i,:}(5,indx1)...
                ,pick_xy_lengths,pick_xy_lengths,exp(est.X{i,:}(9,indx1)));

        end
        if  est.N(i,:) == 0
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
            
        end
    end
    fclose(fileID);
    %%
elseif strcmp(model.dataset, 'CMC1')
    fileID = fopen('./Results/EST_CMC1_WORLD_CENTROID.txt','w');
    for i = meas.start : meas.end
        store_L = [est.L{i}' ];
        [c,~,~]= unique(store_L,'rows');
        store_L = c' ;
        store_X = [];
        for dsa = 1 :size(store_L,2)
            store_X = [];
            indx1 = find(all(est.L{i} == store_L(:,dsa) ));
            %          indx2 = find(all(L2{i} == store_L(:,dsa) ));
            %          indx3 = find(all(L3{i} == store_L(:,dsa) ));
            %          indx4 = find(all(L4{i} == store_L(:,dsa) ));
            
            %             if isempty(indx1)
            %
            %             else
            %                 store_X = [store_X X1{i}(:,indx1)];
            %
            %             end
            %             if isempty(indx2)
            %             else
            %                 store_X = [store_X X2{i}(:,indx2)];
            %             end
            %             if isempty(indx3)
            %             else
            %                 store_X = [store_X X3{i}(:,indx3)];
            %             end
            %             if isempty(indx4)
            %             else
            %                 store_X = [store_X X4{i}(:,indx4)];
            %             end
            
            %         truth.X{i,:}(:,dsa) = sum(store_X,2)/size(store_X,2);
            %         truth.L{i,:}(:,dsa) = store_L(:,dsa);
            %         truth.N(i,:) =  size(store_L,2);
            pick_xy_lengths = max(exp(est.X{i,:}(7,indx1)),exp(est.X{i,:}(8,indx1)));
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,store_L(1,dsa),-1,-1,-1,-1,-1,est.X{i,:}(1,indx1),est.X{i,:}(3,indx1),est.X{i,:}(5,indx1)...
                ,pick_xy_lengths,pick_xy_lengths,exp(est.X{i,:}(9,indx1)));
            
        end
        if  est.N(i,:) == 0
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
            
        end
    end
    fclose(fileID);
    
    
    %%
elseif strcmp(model.dataset, 'CMC4')    
    fileID = fopen('./Results/EST_CMC4_WORLD_CENTROID.txt','w');
    for i = meas.start : meas.end
        store_L = [est.L{i}' ];
        [c,~,~]= unique(store_L,'rows');
        store_L = c' ;
        store_X = [];
        for dsa = 1 :size(store_L,2)
            store_X = [];
            indx1 = find(all(est.L{i} == store_L(:,dsa) ));
            %          indx2 = find(all(L2{i} == store_L(:,dsa) ));
            %          indx3 = find(all(L3{i} == store_L(:,dsa) ));
            %          indx4 = find(all(L4{i} == store_L(:,dsa) ));
            
            %             if isempty(indx1)
            %
            %             else
            %                 store_X = [store_X X1{i}(:,indx1)];
            %
            %             end
            %             if isempty(indx2)
            %             else
            %                 store_X = [store_X X2{i}(:,indx2)];
            %             end
            %             if isempty(indx3)
            %             else
            %                 store_X = [store_X X3{i}(:,indx3)];
            %             end
            %             if isempty(indx4)
            %             else
            %                 store_X = [store_X X4{i}(:,indx4)];
            %             end
            
            %         truth.X{i,:}(:,dsa) = sum(store_X,2)/size(store_X,2);
            %         truth.L{i,:}(:,dsa) = store_L(:,dsa);
            %         truth.N(i,:) =  size(store_L,2);
            pick_xy_lengths = max(exp(est.X{i,:}(7,indx1)),exp(est.X{i,:}(8,indx1)));
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,store_L(1,dsa),-1,-1,-1,-1,-1,est.X{i,:}(1,indx1),est.X{i,:}(3,indx1),est.X{i,:}(5,indx1)...
                ,pick_xy_lengths,pick_xy_lengths,exp(est.X{i,:}(9,indx1)));
            
        end
        if  est.N(i,:) == 0
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
            
        end
    end
    fclose(fileID);
    
    %%
elseif strcmp(model.dataset, 'CMC3')
    fileID = fopen('./Results/EST_CMC3_WORLD_CENTROID.txt','w');
    for i = meas.start : meas.end
        store_L = [est.L{i}' ];
        [c,~,~]= unique(store_L,'rows');
        store_L = c' ;
        store_X = [];
        for dsa = 1 :size(store_L,2)
            store_X = [];
            indx1 = find(all(est.L{i} == store_L(:,dsa) ));
            %          indx2 = find(all(L2{i} == store_L(:,dsa) ));
            %          indx3 = find(all(L3{i} == store_L(:,dsa) ));
            %          indx4 = find(all(L4{i} == store_L(:,dsa) ));
            
            %             if isempty(indx1)
            %
            %             else
            %                 store_X = [store_X X1{i}(:,indx1)];
            %
            %             end
            %             if isempty(indx2)
            %             else
            %                 store_X = [store_X X2{i}(:,indx2)];
            %             end
            %             if isempty(indx3)
            %             else
            %                 store_X = [store_X X3{i}(:,indx3)];
            %             end
            %             if isempty(indx4)
            %             else
            %                 store_X = [store_X X4{i}(:,indx4)];
            %             end
            
            %         truth.X{i,:}(:,dsa) = sum(store_X,2)/size(store_X,2);
            %         truth.L{i,:}(:,dsa) = store_L(:,dsa);
            %         truth.N(i,:) =  size(store_L,2);
            pick_xy_lengths = max(exp(est.X{i,:}(7,indx1)),exp(est.X{i,:}(8,indx1)));
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,store_L(1,dsa),-1,-1,-1,-1,-1,est.X{i,:}(1,indx1),est.X{i,:}(3,indx1),est.X{i,:}(5,indx1)...
                ,pick_xy_lengths,pick_xy_lengths,exp(est.X{i,:}(9,indx1)));            
        end
        if  est.N(i,:) == 0
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
            
        end
    end
    fclose(fileID);

elseif strcmp(model.dataset, 'CMC5')
    fileID = fopen('./Results/EST_CMC5_WORLD_CENTROID.txt','w');
    for i = meas.start : meas.end
        store_L = [est.L{i}' ];
        [c,~,~]= unique(store_L,'rows');
        store_L = c' ;
        store_X = [];
        for dsa = 1 :size(store_L,2)
            store_X = [];
            indx1 = find(all(est.L{i} == store_L(:,dsa) ));
            %          indx2 = find(all(L2{i} == store_L(:,dsa) ));
            %          indx3 = find(all(L3{i} == store_L(:,dsa) ));
            %          indx4 = find(all(L4{i} == store_L(:,dsa) ));
            
            %             if isempty(indx1)
            %
            %             else
            %                 store_X = [store_X X1{i}(:,indx1)];
            %
            %             end
            %             if isempty(indx2)
            %             else
            %                 store_X = [store_X X2{i}(:,indx2)];
            %             end
            %             if isempty(indx3)
            %             else
            %                 store_X = [store_X X3{i}(:,indx3)];
            %             end
            %             if isempty(indx4)
            %             else
            %                 store_X = [store_X X4{i}(:,indx4)];
            %             end
            
            %         truth.X{i,:}(:,dsa) = sum(store_X,2)/size(store_X,2);
            %         truth.L{i,:}(:,dsa) = store_L(:,dsa);
            %         truth.N(i,:) =  size(store_L,2);
            pick_xy_lengths = max(exp(est.X{i,:}(7,indx1)),exp(est.X{i,:}(8,indx1)));
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,store_L(1,dsa),-1,-1,-1,-1,-1,est.X{i,:}(1,indx1),est.X{i,:}(3,indx1),est.X{i,:}(5,indx1)...
                ,pick_xy_lengths,pick_xy_lengths,exp(est.X{i,:}(9,indx1)));            
            
        end
        if  est.N(i,:) == 0
            fprintf(fileID,'%d %d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f\n',i,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
            
        end
    end
    fclose(fileID);
end



