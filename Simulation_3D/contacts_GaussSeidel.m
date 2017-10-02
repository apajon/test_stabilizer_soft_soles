function [Pg,Pc,dPabstot,contact,FcSurf,Ftot,Z,displ,angleact,Kcart]=contacts_GaussSeidel(sole,P0,center,displ,Pg,fric,angleact,contAngle,FcSurf,Fdes,Zdes)
% Contact with friction based on the works of Christian Duriez, Frédéric Dubois,
% Abderrahmane Kheddar and Claude Andriot
% "Realistic Haptic Rendering of Interacting Deformable Obects in Virtual
% Environments"
D = 3; % Dimension of the contact problem
% theta = angle(1);
% phi = angle(2);
m = sole.nFreeSurf;
Ccc = sole.Cs;

[displ,angleact,Fc_mc3,ind_Cont,ind_Cont3,Ftot,Z,Kcart] = PsoleAngle_GaussSeidel(m,Ccc,sole,P0,Pg,fric,contAngle,FcSurf,displ,angleact,Fdes,Zdes);
theta = angleact(1);
phi = angleact(2);
psi = angleact(3);
R = Rot(theta,phi,psi);
%%% Find the displacement of internal nodes and position of all nodes
FcSurf3 = zeros(D*sole.nFreeSurf,1);
FcSurf3(ind_Cont3) = Fc_mc3;
FcSurf = reshape(FcSurf3,D,m)';
Fcloc = R'*FcSurf';
Fcloc3 = reshape(Fcloc,D*m,1);
dPlocSurf3 = Ccc * Fcloc3;
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
% Position of contacts
Pc = Pg(sole.nodesFreeSurf(ind_Cont),:);
contact = ind_Cont;