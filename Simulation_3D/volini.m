function vol = volini(sole)
detJ = zeros(size(sole.elements_vol,1),1);
vol = 0;
for j = 1:size(sole.elements_vol,1)
    vertices = sole.coor(sole.elements_vol(j,:),:);
    detJ(j) = det([1,1,1,1;vertices']);
    vol = vol + abs(detJ(j)/6);
end