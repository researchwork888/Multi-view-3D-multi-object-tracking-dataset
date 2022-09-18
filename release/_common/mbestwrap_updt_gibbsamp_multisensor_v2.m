function [assignments,costs]= mbestwrap_updt_gibbsamp_multisensor_v2(P0,m,p)
    

n1 = size(P0{1},1);
n2 = size(P0,2);
no_sensors = size(P0,2);

assignments= zeros(m,n1);
costs= zeros(m,1)';

currsoln= n1+1:2*n1; %use all missed detections as initial solution
assignments(1,:)= currsoln;
costs(1)=sum(P0{1}(sub2ind(size(P0{1}),1:n1,currsoln)));


for sol= 2:m 
    for var= 1:n1 % this is looping through all labels 
        tempsamp = [];
        for q = 1 : no_sensors
            temp_asd = exp(-P0{q}(var,:)); %grab row of costs for current association variable (grab label by label)
            tempsamp = cat(2,tempsamp,temp_asd);
        end
        tepsamp(currsoln([1:var-1,var+1:end]))= 0; %lock out current and previous iteration step assignments except for the one in question
        idxold= find(tempsamp>0); 
        tempsamp= tempsamp(idxold);
        [~,currsoln(var)]= histc(rand(1,1),[0;cumsum(tempsamp(:))/sum(tempsamp)]); 
        currsoln(var)= idxold(currsoln(var));
    end
       assignments(sol,:)= currsoln;
    costs(sol)= sum(P0(sub2ind(size(P0),1:n1,currsoln)));
end
[C,I,~]= unique(assignments,'rows');
assignments= C;
costs= costs(I);

end