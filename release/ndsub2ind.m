function linidx= ndsub2ind(siz,idx)

if isempty(idx)
    linidx= [];
else
    linidx=idx(:,1);
    totals=cumprod(siz);
    for i = 2:size(idx,2)
        linidx = linidx + (idx(:,i)-1)*totals(i-1);
    end
end

