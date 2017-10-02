function [Pg,ind_Cont,FcFreeSurf,Ftot,Z,displ,angleact,Kcart,no_conv,dPabstot,J,Mo]=contacts(sole,P0,center,displ,Pg,fric,angleact,contAngle,FcFreeSurf,Fdes,Zdes)
% Contact with friction based on the works of Christian Duriez, Frédéric Dubois,
% Abderrahmane Kheddar and Claude Andriot
% "Realistic Haptic Rendering of Interacting Deformable Obects in Virtual
% Environments"
D = 3; % Dimension of the contact problem
% theta = angle(1);
% phi = angle(2);
m = sole.nFreeSurf;
Pgini = Pg;

[displ,angleact,Fc_mc3,PcSurf3,ind_Cont,ind_Cont3,Ftot,Z,Kcart,no_conv,J,Mo] = PsoleAngle(m,sole.Cs,sole,P0,Pg,fric,contAngle,FcFreeSurf,displ,angleact,Fdes,Zdes);
save 'Kcart.mat' Kcart
% if contAngle ~= 1
%     PcSurf = reshape(PcSurf3,D,m)';
%     figure(8); clf;
%     scatter3(PcSurf(:,1),PcSurf(:,2),PcSurf(:,3),'b')
%     hold on
% end
% PcSurf = reshape(PcSurf3,D,m)';
% theta = angleact(1);
% phi = angleact(2);
% psi = angleact(3);
% R = Rot(theta,phi,psi);
% mc = length(ind_Cont);
% % ind_NonCont3 = setdiff(sole.nodesFreeSurf_noDir3,ind_Cont3);
% 
% if contAngle==1
%     dPsurf = PcSurf - P0(:,sole.nodesFreeSurf)';
% else
%     dPsurf = PcSurf - Pgini(sole.nodesFreeSurf,:);
% end    
% PosRot = R*P0 + center;
% PosRot = PosRot';
% displacement = repmat(displ', sole.nTot, 1);
% Pg = PosRot + sole.getdP(dPsurf) + displacement;
% Pg(sole.nodesFreeSurf,:) = PcSurf;
% Fc_mc = reshape(Fc_mc3,D,mc)';
% FcFreeSurf = zeros(D,m);
% FcFreeSurf(:,ind_Cont) = Fc_mc';
% FcFreeSurf = FcFreeSurf';

theta = angleact(1);
phi = angleact(2);
psi = angleact(3);
R = Rot(theta,phi,psi);
%%% Find the displacement of internal nodes and position of all nodes
FcSurf3 = zeros(D*sole.nFreeSurf,1);
FcSurf3(ind_Cont3) = Fc_mc3;
FcFreeSurf = reshape(FcSurf3,D,m)';
Fcloc = R'*FcFreeSurf';
Fcloc3 = reshape(Fcloc,D*m,1);
dPlocSurf3 = sole.Cs * Fcloc3;
% Displacement of all nodes
dPInt3 = -(sole.m_invKii_Kis) * dPlocSurf3;
dPloc3 = zeros(D*sole.nTot,1);
dPloc3(sole.nodesFreeSurf3) = dPlocSurf3;
dPloc3(sole.nodesInt3) = dPInt3;
dPloc = reshape(dPloc3,D,sole.nTot)';
dPabstot = R*dPloc';
displacement = repmat(displ', sole.nTot, 1);
PosRot = R*P0 + center;
Pg = PosRot' + dPabstot' + displacement;
% figure(5); clf;
% scatter3(Pg(:,1),Pg(:,2),Pg(:,3),'g')
% hold on
% 
% scatter3(Pg(sole.nodesInt,1),Pg(sole.nodesInt,2),Pg(sole.nodesInt,3),'r')
% hold on
% RPfree_ntot = R*P0;
% RPfree_nf = RPfree_ntot(:,sole.nodesFree);
% 
% %delta_mc3 = displ_mc3 + RPfree_mc(:) + W_mc * Fc_mc3;
% % RtFc_mc3 = R'*reshape(Fc_mc3,D,mc);
% % RtFc_mc3 = RtFc_mc3(:);
% % Ccc_Cont3 = Ccc(ind_Cont3,ind_Cont3);
% % CccRtFc_mc3 = Ccc_Cont3 * RtFc_mc3;
% % RCccRtFc_mc3 = R*reshape(CccRtFc_mc3,D,mc);
% % RCccRtFc_mc3 = RCccRtFc_mc3(:);
% % delta_mc3 = displ_mc3 + RPfree_mc(:) + RCccRtFc_mc3;
% RtFcFree = zeros(D,sole.nFree);
% Fc_mc = reshape(Fc_mc3,D,mc)';
% Nc_free = sole.nodesFreeSurf(ind_Cont);
% RtFcFree(:,Nc_free) = R' * Fc_mc';
% RtFcFree3 = RtFcFree(:);
% CccRtFcFree3 = sole.C * RtFcFree3;
% 
% for i=1:sole.nTot
%     
% end

% 
% Knc_nc = sole.Ks(ind_NonCont3,ind_NonCont3);
% Knc_c = sole.Ks(ind_NonCont3,ind_Cont3);
% Kc_nc = sole.Ks(ind_Cont3,ind_NonCont3);
% Kc_c = sole.Ks(ind_Cont3,ind_Cont3);
% invKnc_nc_Knc_c = inv(Knc_nc)*Knc_c;
% Kc = (Kc_c - Kc_nc*invKnc_nc_Knc_c);
% % m_invKConCon_KNconCon = full(KNconNcon\KNconCon);
% %%% Find the displacement of internal nodes and position of all nodes
% % FcSurf3 = zeros(D*sole.nFreeSurf,1);
% % FcSurf3(ind_Cont3) = Fc_mc3;
% % FcSurf = reshape(FcSurf3,D,m)';
% % Fcloc = R'*FcSurf';
% % Fcloc3 = reshape(Fcloc,D*m,1);
% % FcontLoc = ind_Cont
% Fc_mc = reshape(Fc_mc3,D,mc)';
% FcFreeSurf = zeros(D,m);
% FcFreeSurf(:,ind_Cont) = Fc_mc';
% FcFreeSurf = FcFreeSurf';
% FcLoc = R' * Fc_mc';
% FcLoc3 = reshape(FcLoc,D*mc,1);
% % dPcLoc = inv(Kc_c)*FcLoc3;
% dPcLoc = sole.Cs(ind_Cont3,ind_Cont3)*FcLoc3;
% dPNcLoc = -invKnc_nc_Knc_c * dPcLoc;
% %dPNcLoc = sole.Cs(ind_NonCont3,ind_Cont3) * FcLoc3;
% dPlocSurf3(ind_Cont3,1) = dPcLoc;
% dPlocSurf3(ind_NonCont3,1) = dPNcLoc;
% % dPlocSurf3 = sole.Cs * Fcloc3;
% 
% %%%%%%% Maybe is like that. If not delete this line
% % dPlocSurf3(ind_NonCont3) = -(m_invKConCon_KNconCon) * dPlocSurf3(ind_Cont3);
% 
% % Displacement of all nodes
% dPInt3 = -(sole.m_invKii_Kis) * dPlocSurf3;
% dPloc3 = zeros(D*sole.nTot,1);
% dPloc3(sole.nodesFreeSurf3) = dPlocSurf3;
% dPloc3(sole.nodesInt3) = dPInt3;
% dPloc = reshape(dPloc3,D,sole.nTot)';
% dPabstot = R*dPloc';
% displacement = repmat(displ', sole.nTot, 1);
% PosRot = R*P0 + center;
% Pg = PosRot' + dPabstot' + displacement;
% Pg(sole.nodesFreeSurf,:)=reshape(PcSurf,D,m)';