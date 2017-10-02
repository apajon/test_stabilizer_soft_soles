function [dOl_dPree_mc,dY_dPree_mc] = compute_gradient(sole,ind_cont,angleact,PSurf3,Fc_mc3,PabsOld_mc,ind_slip,Zdes,fric,FcOld3,displ_first,psi_first,Pfree_c,Cs,contAngle)
% dOl = ddispl
% dY = dangleact
D=3;
ind_cont3 = sort(([D*ind_cont-2; D*ind_cont-1; D*ind_cont]));
ind_cont3 = reshape(ind_cont3,D*length(ind_cont),1); 
theta = angleact(1);
phi = angleact(2);
psi = angleact(3);
R = Rot(theta,phi,psi);
RBlock_c = [];
for i=1:length(ind_cont)
    RBlock_c = blkdiag(RBlock_c,R);
end

Rini = Rot(theta,phi,psi_first);
% RBlock_old_c = [];
% for i=1:length(ind_cont)
%     RBlock_old_c = blkdiag(RBlock_old_c,ROld);
% end
    
P_mc3 = PSurf3(ind_cont3);
P_mc = reshape(P_mc3,3,length(ind_cont));
% FcOld3_mc= FcOld3(ind_cont3);
% PabsOld_mc = PabsOld(:,ind_cont);
Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
Wc = zeros(3*length(ind_cont),3*length(ind_cont));
for i = 1:length(ind_cont)
    for j = 1:length(ind_cont)
        Wc((3*i)-2:3*i,(3*j)-2:3*j) = R * Cs((3*ind_cont(i))-2:3*ind_cont(i),(3*ind_cont(j))-2:3*ind_cont(j)) * R';
    end
end
der_C_s = cell(3*sole.nFreeSurf,1);
KNodesFreeSurf = reshape(sole.dof(sole.nodesFreeSurf,:)',[1,3*sole.nFreeSurf]);
for i = 1:(3*sole.nFreeSurf)      
    der_C_s{i} = sole.der_C{KNodesFreeSurf(i)}(KNodesFreeSurf,KNodesFreeSurf);
end

der_C_cc = cell(3*length(ind_cont),1);
for i = 1:(3*length(ind_cont))
    der_C_cc{i} = der_C_s{ind_cont3(i)}(ind_cont3,ind_cont3);
end
%%% Gradient Simulation
E = gradient_simu_E_contact(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
Ccc = Cs(ind_cont3,ind_cont3);
B = gradient_simu_B_contact(angleact,P_mc3,Fc_mc3,ind_cont,Ccc,ind_slip,psi_first,contAngle);
% A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric);
A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric,contAngle,displ_first,psi_first,angleact);

Lambda = gradient_simu_Lambda_contact(ind_cont,der_C_cc,ind_slip,Fc_mc3,RBlock_c);
dF_delta_mc = cell(3*length(ind_cont),1);
dF_dPree_mc = cell(3*length(ind_cont),1);
ddeltat_dPree_mc = cell(3*length(ind_cont),1);
% dOl_dPree_mc = cell(3*length(ind_cont),1);

%dY_dPree_mc = cell(3*length(ind_cont),1);
dP_dPfree_mc = cell(3*length(ind_cont),1);
% dPfree_dp_mc = dPfree_dp(ind_cont3,:);
for i=1:(3*length(ind_cont))
    Grad1 = inv(E*inv(A)*B)*(-E*inv(A)*Lambda);
    dOl_dPree_mc = Grad1(1:3,:);
    dY_dPree_mc = Grad1(4:6,:);
    %dOl_dPree_mc = repmat(Grad1(1:3,:),length(ind_cont),1);
    %dY_dPree_mc = Grad1(4:6,:);
%     dF_delta_mc{i} = inv(A) * B * [dOl_dPree_mc{i};dY_dPree_mc{i}] + Lambda{i};
%     dF_dPree_mc{i} = dF_delta_mc{i}(1:(3*length(ind_cont)),:);
%     ddeltat_dPree_mc{i} = dF_delta_mc{i}((3*length(ind_cont)+1):end,:);
%     dP_dPfree_mc{i} = dY_dPree_mc{i} * Ccc * RBlock_c' * Fc_mc3 + RBlock_c * der_C_cc{i} * RBlock_c' * Fc_mc3 + Wc * dF_dPree_mc{i} + dOl_dPree_mc{i};
%     
end
end