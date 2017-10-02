function cost = optimization(p_mid,sole,friction,Fdes,Zdes,spl,coorini)
% a = load('polynomer.mat');
% p_mid = a.p_mid;
% if norm(p_mid)<7.9 || norm(p_mid)>8.1 
%     p_mid
% end
coornew = deformation(size(sole.coor),p_mid,spl);

%%% translate because spl is in coornew-sole.trasl
coornew(:,1)= coornew(:,1) + sole.trasl(1);
coornew(:,2) = coornew(:,2) + sole.trasl(2);
coornew(:,3) = coornew(:,3) + sole.trasl(3) + sole.zpoles;
%coornew(sole.nodesDirichlet,3) = sole.coor(sole.nodesDirichlet,3); % reset z coordinater of Dirichlet nodes

dPsurf = coornew(sole.nodesFreeSurf,:) - coorini(sole.nodesFreeSurf,:);
coornew = coorini + sole.getdP(dPsurf);

% figure(4); clf;
% scatter3(coornew(:,1),coornew(:,2),coornew(:,3),'r')
% hold on;
% scatter3(coorini(:,1),coorini(:,2),coorini(:,3),'g')

sole.coor = coornew;

sole.stiffness();
sole.stiffnessSurface();
% 
% dPsurf = coornew(sole.nodesFreeSurf,:) - sole.coor(sole.nodesFreeSurf,:);
% coornew = sole.coor + sole.getdP(dPsurf);

% %sole.coor = coornew;
% % if sole.IteOpt>1
% %     a = load('polynome.mat');
% %     p_mid_Old = a.p_mid;
% %     if norm(p_mid-p_mid_Old)>0.1
%         pname = 'input/semelle1 L=0.23, l=0.13, e=0.03 m new centre/';
%         coornewtmp = remesh(sole,pname,coornew);
%         [D,I] = pdist2(coornew,coornewtmp,'euclidean','Smallest',1);       
%         sole.coor(I,:) = coornewtmp; 
% %     else
% sole.stiffness();
% sole.stiffnessSurface();   
        
%     end
% end

% stressVM0 = zeros(sole.nTot,1);
% plotsole(4,sole.elements_surf,sole.coor,stressVM0,[],[],[],[]);  
    

save 'polynome.mat' p_mid
cost = simulationOpt(sole,friction,Fdes,Zdes);
end