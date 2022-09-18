function [assignments_a,costs_a]= mbestwrap_updt_gibbsamp_multisensor(P0,m,p)
    

n1 = size(P0{1},1);
n2 = size(P0,2);
no_sensors = size(P0,2);

assignments= zeros(no_sensors,n1,m);
costs= zeros(m,1)';

currsoln= repmat(n1+1:2*n1,no_sensors,[1 1]); %use all missed detections as initial solution
assignments(:,:,1)= currsoln;
temp = 0; 
for q = 1: no_sensors
    temp = temp + sum(P0{q}(sub2ind(size(P0{q}),1:n1,currsoln(q,:))));   
end
costs(1) = temp; 
clear temp;

for sol= 2:m
    for var= 1:n1 % this is looping through all labels 
        for q = 1 : no_sensors
            tempsamp= exp(-P0{q}(var,:)); %grab row of costs for current association variable (grab label by label)
            tempsamp(currsoln(q,[1:var-1,var+1:end]))= 0; %lock out current and previous iteration step assignments except for the one in question
            idxold= find(tempsamp>0); 
            tempsamp= tempsamp(idxold);
            [~,currsoln(q,var)]= histc(rand(1,1),[0;cumsum(tempsamp(:))/sum(tempsamp)]); 
            currsoln(q,var)= idxold(currsoln(q,var));
        end
        
    end
    assignments(:,:,sol)= currsoln;
    temp = 0; 
    for q = 1: no_sensors
        temp = temp + sum(P0{q}(sub2ind(size(P0{q}),1:n1,currsoln(q,:))));   
    end
    costs(sol) = temp; 
    clear temp;
end


check_flag = ones(1,size(assignments,3));
temp_flag = 1:size(assignments,3);
indx = [];
flag  = true;
for i = 1:size(assignments,3)
%     temp = assignments(:,:,count);
    if check_flag(1,i) == 1
        check_flag(1,i) = 0;
        check_array = temp_flag.*check_flag;
    for count = 1: size(check_array,2)
        if check_flag(1,count) == 0
            
        else 
            if all(all( assignments(:,:,i) ==  assignments(:,:,count)))
                check_flag(1,count) = 0;
                if flag 
                    indx = [indx; i ];
                    flag = false; 
                end
            end
        end

    end
                if flag 
                    indx =[indx; i ];
                    flag = false; 

                end
    flag = true; 
    end
end
% [C,I,~]= unique(assignments,'rows');
assignments_a= assignments(:,:,indx);
costs_a= costs(indx);



