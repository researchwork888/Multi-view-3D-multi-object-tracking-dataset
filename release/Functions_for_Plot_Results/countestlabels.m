function [count,c]= countestlabels(meas,est)
    labelstack= [];
    for k=meas.start:meas.end
        labelstack= [labelstack est.L{k}];
    end
    [c,~,~]= unique(labelstack','rows');
    count=size(c,1);
end