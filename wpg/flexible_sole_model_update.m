function [sole,Z,Ftot,param_sopt,ABe,Pg,FcFreeSurf,D,Mz]=flexible_sole_model_update(sole,angleact,displ,contAngle,param_sopt,ABe,Pg,FcFreeSurf,D,need_plot)

addpath ./Model_3D/
% addpath ./Model_3D/input
addpath ./Model_3D/FEM
addpath ./Model_3D/gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PsoleAngle.m
m = sole.nFreeSurf;
Cs = sole.Cs;
no_conv = [];
P0 = sole.coor';
if contAngle==1
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld(1,:) = P0(1,sole.nodesFreeSurf);
    PabsOld(2,:) = P0(2,sole.nodesFreeSurf);
    PabsOld(3,:) = P0(3,sole.nodesFreeSurf);    
    FcNew3 = zeros(D*m,1);
    displ_first= [0;0;0];
    psi_first = 0;
else
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld = Pg(sole.nodesFreeSurf,:)';
    FcNew = FcFreeSurf'; 
    FcNew3 = reshape(FcNew,D*m,1);
    displ_first= [0;0;0];
    psi_first = 0;
end
FcOld3 = FcNew3;
FcOld = reshape(FcOld3,3,sole.nFreeSurf);
diplini = displ;
angleactini = angleact;
angleactOld = angleact;

[displ,angleact,FtotZMP,Fc_mc3_out,ind_cont_out,ind_slip_out,Mzmp,PSurf3] = Gauss(contAngle,param_sopt.friction,m,displ,angleact,FcNew3,Pfree,PabsOld,Cs);

if ind_cont_out(1)==0
    a = find(ind_cont_out==0,2);
    Fc_mc3 = Fc_mc3_out(1:D*(a(2)-1));
    ind_cont = (sort(ind_cont_out(1:(a(2)-1))+1))';      
else
    a = find(ind_cont_out==0,1);
    Fc_mc3 = Fc_mc3_out(1:D*(a-1));
    ind_cont = (sort(ind_cont_out(1:(a-1))+1))';              
end
if norm(ind_slip_out)>0
    if ind_slip_out(1)==0
        b = find(ind_slip_out==0,2);
        ind_slip = (sort(ind_slip_out(1:(b(2)-1))+1))';  
    else
        b = find(ind_slip_out==0,1);
        ind_slip = (sort(ind_slip_out(1:(b-1))+1))';      
    end
else
    ind_slip = [];
end

Ftot = FtotZMP(1:3);
Z = FtotZMP(4:5);
Mz=Mzmp;
ind_cont3 = sort(([D*ind_cont-2; D*ind_cont-1; D*ind_cont]));
ind_cont3 = reshape(ind_cont3,D*length(ind_cont),1); 
Pfree_c = Pfree(:,ind_cont);
Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
P_mc3 = PSurf3(ind_cont3);
P_mc = reshape(P_mc3,3,length(ind_cont));
theta = angleact(1);
phi = angleact(2);
psi = angleact(3);
R = Rot(theta,phi,psi);
Wc = zeros(3*length(ind_cont),3*length(ind_cont));
for j = 1:length(ind_cont)
    for k = 1:length(ind_cont)
        Wc((3*j)-2:3*j,(3*k)-2:3*k) = R * Cs((3*ind_cont(j))-2:3*ind_cont(j),(3*ind_cont(k))-2:3*ind_cont(k)) * R';
    end
end
PabsOld_mc = PabsOld(:,ind_cont);
A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,param_sopt.friction,contAngle,displ_first,psi_first,angleact);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_cont3 = sort([D*ind_cont-2 D*ind_cont-1 D*ind_cont]);
ind_cont3 = reshape(ind_cont3,D*length(ind_cont),1); 
theta = angleact(1);
phi = angleact(2);
psi = angleact(3);
R = Rot(theta,phi,psi);
%%% Find the displacement of internal nodes and position of all nodes
FcSurf3 = zeros(D*sole.nFreeSurf,1);
FcSurf3(ind_cont3) = Fc_mc3;
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
% PosRot = R*P0 + center;
PosRot = R*sole.coor';
Pg = PosRot' + dPabstot' + displacement;

Fc = FcFreeSurf(ind_cont,:);
Pc = Pg(sole.nodesFreeSurf(ind_cont),:);

stressVM = stressVonMises(sole,dPabstot,ABe); % VonMises' Stress
if need_plot
    plotsole(2,sole.elements_surf,Pg,stressVM,Pc,Fc,Z,Ftot,-127.5,30);
end

rmpath ./Model_3D/
% rmpath ./Model_3D/input
rmpath ./Model_3D/FEM
rmpath ./Model_3D/gradient