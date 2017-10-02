function coornew = deformationFE(s_coor,p_mid,spl)
% p = [ones(1,13); ones(11,1) p_mid ones(11,1); ones(1,13)];
% p = [ones(1,23); ones(21,1) p_mid ones(21,1); ones(1,23)];
% With fix Dirichlet node
p = [ones(1,size(p_mid,1)+2); ones(size(p_mid,1),1) p_mid ones(size(p_mid,1),1); ones(1,size(p_mid,1)+2)];
% With moving Dirichlet node
%p = p_mid;
coornew = zeros(s_coor);
coor_tmp = zeros(1,3);
for i=1:s_coor(1)
    s1 = sum(sum(p(spl(i).idx_u, spl(i).idx_v).*spl(i).basisValues));
    [coor_tmp(1,1),coor_tmp(1,2),coor_tmp(1,3)] = sph2cart(spl(i).u,spl(i).v,s1*spl(i).r);
    coornew(i,:) = [coor_tmp(1,3),coor_tmp(1,2),-coor_tmp(1,1)];    
end
% coornew(78,1)=coornew(89,1);
% dPsurf = coornew(sole.nodesFreeSurf,:) - sole.coor(sole.nodesFreeSurf,:);
% coornew = sole.coor + sole.getdP(dPsurf);
% scatter3(coornew(:,1),coornew(:,2),coornew(:,3),'r')
% hold on
end