function stressVM = stressVonMises(sole,dP,ABe)
%--------------------------------------------------------------------------
%         Computation of the VonMises' Stress
%--------------------------------------------------------------------------
displ = zeros(sole.nTot,3);
dP(sole.nodesDirichlet3) = [];  
ic = find(sole.dof>0);                              % finds active dofs

displ(ic) = dP(sole.dof(ic));

stress = zeros(sole.nTot,6);                 % initializes "stress"
stressG2 = zeros(4,6);
counter = zeros(sole.nTot,1);                % initializes "counter"
for e = 1:size(sole.elements_vol,1)                      % computes nodal stresses
    nodes = sole.elements_vol(e,:);
    stressG2(1,:) = [displ(nodes(1),:) displ(nodes(2),:) displ(nodes(3),:) displ(nodes(4),:)]*ABe{e}';
    stressG2(2,:) = stressG2(1,:);
    stressG2(3:4,:) = stressG2(1:2,:);
    stress(nodes,:) = stress(nodes,:) + stressG2;   % adds stresses to triangle nodes
    counter(nodes) = counter(nodes) + 1;        
end
for icomp=1:6
    stress(:,icomp)=stress(:,icomp)./counter;  % naive average of stresses
end

%--------------------------------------------------------------------------
%             Calcul de la contrainte equivalente de Von Mises
%--------------------------------------------------------------------------
stressVM=sqrt((3/2)*(((2*stress(:,1)-stress(:,2)-stress(:,3))/3).^2 ...
        +((2*stress(:,2)-stress(:,1)-stress(:,3))/3).^2 ...
        +((2*stress(:,3)-stress(:,2)-stress(:,1))/3).^2 ...
        +(stress(:,4)).^2+(stress(:,5)).^2+(stress(:,6)).^2));
end
