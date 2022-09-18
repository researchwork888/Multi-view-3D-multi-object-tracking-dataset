function bbs =   convert_bbox_to_bbs(bbox)
bbs(1,:)  = bbox(1,:);
bbs(2,:)  = bbox(2,:);
bbs(3,:) =  bbox(3,:) - bbox(1,:);
bbs(4,:) =  bbox(4,:) - bbox(2,:);

end
