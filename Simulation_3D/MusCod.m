format long
clc
% clear all
close all
addpath ./input
addpath ./t4
addpath ./utility

D = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Sole FEM                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'semelle1.inp';
pname = 'input/semelle1 L=0.25, l=0.19, e=0.03 m new centre/';
 
sole = soleFEM_newStiff(pname,fname);
friction = 0.8;
sole.stiffness(sole.coor);
sole.stiffnessSurface();
m = sole.nFreeSurf;
Ccc = sole.Cs;

contAngle = 1;
displ = [-0.01;0.1;-0.1];
angleact = [degtorad(0);degtorad(10);degtorad(0)];
% [VarI1, VarI2, VarI3] = textread('/media/63DFD35C3FD8A3D0/Giappone/Code/sim3d/Simulation 3D - Desired ZMP - Position and Force - 3 Angles - meters/RightFootPosAndOri.txt', '%f %f %f', 'headerlines', 0);
% displ1 = [VarI1(1);VarI2(1);VarI3(1)];
% R = [VarI1(2) VarI2(2) VarI3(2);VarI1(3) VarI2(3) VarI3(3);VarI1(4) VarI2(4) VarI3(4)];
FcNew3 = zeros(D*m,1);
%acos(VarI1(2))
% angleact1 = [0;acos(VarI1(2));0];
Pfree = sole.coor(sole.nodesFreeSurf,:)';
PabsOld = Pfree;
[displ,angleact,FtotZMP,Fc_mc3_out,ind_cont_out] = Gauss(contAngle,friction,m,displ,angleact,[0;0;0],FcNew3,Pfree,PabsOld,Ccc);
Ftot = FtotZMP(1:3);
Z = FtotZMP(4:6);

fileID = fopen('ForceAndZMP.txt','w');
fprintf(fileID,'%f, %f, %f,\n',Ftot);
fprintf(fileID,'%f, %f, %f,\n',Z);
fclose(fileID);

exit;