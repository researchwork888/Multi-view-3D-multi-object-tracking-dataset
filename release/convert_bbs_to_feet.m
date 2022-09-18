function feet = convert_bbs_to_feet(bbs)
    if size(bbs,2) == 0 
        if ~(size(bbs,1) == 2) 
            feet = bbs(1:2,:);
        else
            feet = bbs;
        end
    else 
        
    feet(1,:) = bbs(1,:) + (bbs(3,:)/2);
    feet(2,:) = bbs(2,:) + bbs(4,:);
    end
end