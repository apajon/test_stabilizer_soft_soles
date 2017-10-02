format long
clc
% clear all
close all
addpath ./input
addpath ./results
addpath ./FEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Sole FEM                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'semelle1.msh';
pname = 'input/semelle1 L=0.23, l=0.13, e=0.021 m new centre/';
%%% foot size %%%
l=0.13;
L=0.23;
e=0.021;
sole = soleFEM_newStiff(pname,fname,l,L,e);
coorini = sole.coor;
friction = 0.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ZMP                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linux
% [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('/media/63DFD35C3FD8A3D0/Giappone/Code/sim3d/Simulation 3D - Desired ZMP - Position and Force - 3 Angles - meters - Linux/trajectory/exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
% Windows
if exist('trajectory/exemple_trajectoire.txt','file')~=0
    % Linux
    [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('trajectory/exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
elseif exist('trajectory\exemple_trajectoire.txt','file')~=0
    % Windows
    [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('trajectory\exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
end
[Zdes,Fdes] = changeRef(ZMPx, ZMPy, Fx, Fy, Fz);
% Fdes(2,:)=0;
% Zdes(2,:)=0;

%Zdes(:,end-1:end)=[];
% Zdes(:,1)=[];
Zdes = Zdes(:,1:1:length(Zdes));
% Zdes(:,1:5)=[];
% Zdes(:,end-1:end)=[];

%Fdes(:,end-1:end)=[];
% Fdes(:,1)=[];
Fdes = Fdes(:,1:1:length(Fdes));
% Fdes(:,1:5)=[];
% Fdes(:,end-1:end)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           B-spline                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coor = zeros(sole.nTot,3);
coor(:,1)= sole.coor(:,1) - sole.trasl(1);
coor(:,2) = sole.coor(:,2) - sole.trasl(2);
coor(:,3) = sole.coor(:,3) - sole.trasl(3) - sole.zpoles;
azimuth = zeros(1,size(coor,1));
elevation = zeros(1,size(coor,1));
r = zeros(1,size(coor,1));
for i=1:size(sole.coor,1);
    [azimuth(1,i),elevation(1,i),r(1,i)] = cart2sph(-coor(i,3),coor(i,2),coor(i,1));
end
spline_res = -pi/2:pi/5:pi/2;
spl = splineBasis([azimuth;elevation;r], spline_res, spline_res);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Shape Optimization                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p_mid = ones(length(spline_res)+2,length(spline_res)+2);
p_mid = ones(length(spline_res),length(spline_res));
%a = load('polynomesmoothdetailedmaterial2.mat');
% % % % % %  a = load('polynome24-03-2015-1st-part.mat');
% % % % % % % a = load('polynome24-03-2015-2nd-part.mat');
%p_mid = a.p_mid;
% % % % p_mid = p_mid - 0.05;
%p_mid_Old = p_mid;
% % % 
% b = load('polynomestraightdetailed1-25.mat');
% p_mid2 = b.p_mid;
% %p_mid = b.p_mid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Remesh                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coornew = deformation(size(sole.coor),p_mid,spl);
% %coornew = deformation(size(sole.coor),p_mid2-(ones(length(spline_res)+2,length(spline_res)+2)-p_mid),spl);
% coornew(:,1)= coornew(:,1) + sole.trasl(1);
% coornew(:,2) = coornew(:,2) + sole.trasl(2);
% coornew(:,3) = coornew(:,3) + sole.trasl(3);
% coornewtmp = remesh(sole,pname,coornew);
% [D,I] = pdist2(coornew,coornewtmp,'euclidean','Smallest',1);       
% sole.coor(I,:) = coornewtmp; 
% writemshfile(sole,pname)
% stressVM0 = zeros(sole.nTot,1);
% plotsole(1,sole.elements_surf,sole.coor,stressVM0,[],[],[],[]);
% 
% 
% coor = zeros(sole.nTot,3);
% coor(:,1) = sole.coor(:,1) - sole.trasl(1);
% coor(:,2) = sole.coor(:,2) - sole.trasl(2);
% coor(:,3) = sole.coor(:,3) - sole.trasl(3);
% azimuth = zeros(1,size(coor,1));
% elevation = zeros(1,size(coor,1));
% r = zeros(1,size(coor,1));
% for i=1:size(sole.coor,1);
%     [azimuth(1,i),elevation(1,i),r(1,i)] = cart2sph(-coor(i,3),coor(i,2),coor(i,1));
% end
% spl = splineBasis([azimuth;elevation;r], spline_res, spline_res);
% stressVM0 = zeros(sole.nTot,1);
% plotsole(2,sole.elements_surf,sole.coor,stressVM0,[],[],[],[]);

%p_mid = ones(length(spline_res)+2,length(spline_res)+2);

% Poisson ratio caucciu = 0.5
% Young Modulus rubber = 10^7 Pa; rubber = 0.01-0.1 * 10^9 Pa
% Elastomer
% Butyl Rubber 0.001-0.002 GPa
% Sylicon Elastomers 0.005-0.02 GPa
% Neoprene (CR) 0.0007-0.002 GPa
% Neoprene (CR)
%Young = 1500000;
Young = 700000;
Poisson = 0.3;
sole.setMaterial(Young,Poisson);

stressVM0 = zeros(sole.nTot,1);
plotsole(1,sole.elements_surf,sole.coor,stressVM0,[],[],[],[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Optimization                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%u0 = (zeros(length(spline_res)+2,length(spline_res)+2)+1.4)+(ones(length(spline_res)+2,length(spline_res)+2)-p_mid_Old);
%l0 = (zeros(length(spline_res)+2,length(spline_res)+2)+0.4)+(ones(length(spline_res)+2,length(spline_res)+2)-p_mid_Old);
%u0 = (zeros(length(spline_res)+2,length(spline_res)+2)+1.05);
%l0 = (zeros(length(spline_res)+2,length(spline_res)+2)+0.95);
u0 = (zeros(length(spline_res),length(spline_res))+1.5);
l0 = (zeros(length(spline_res),length(spline_res))+0.4);


%optimization(p_mid,sole,soleini,friction,Fdes,Zdes,spl)
%fmincon stopped because it exceeded the function evaluation limit,
%options.MaxFunEvals = 3000 (the default value).
%options = optimset('Display','iter','Algorithm','Interior-Point','MaxFunEvals',3000000000000000000000000000000000,'FinDiffType','central');
options = optimset('Display','iter','Algorithm','Interior-Point','MaxFunEvals',3000000000000000000000000000000000,'TolX',1e-4,'TolFun',1e-4);
%options = optimset('Display','iter','Algorithm','Interior-Point','SubproblemAlgorithm','cg');
% options = optimset('Display','iter','Algorithm','sqp','TolCon',1e-10);
% options = optimset('Display','iter','Algorithm','active-set');
%options = optimset('Display','iter','Algorithm','active-set','DiffMinChange',1e-8,'RelLineSrchBnd',0.03,'RelLineSrchBndDuration',50);
% options = optimset('Display','iter','Algorithm','active-set','DiffMinChange',1e-8,'RelLineSrchBnd',0.03);
%options = optimset('Display','iter','Algorithm','active-set','RelLineSrchBnd',3,'RelLineSrchBndDuration',2);
%options = optimset('Display','iter','Algorithm','interior-point','RelLineSrchBnd',0.03);
% vol = volini(sole);
% X = fmincon(@(p_mid)costFunc(p_mid,sole,friction,Fdes,Zdes,spl,coorini),p_mid,[],[],[],[],l0,u0,@(p_mid)fminconstr_detJ(p_mid,sole,spl,vol),options);
%X = fmincon(@(p_mid)costFunc(p_mid,sole,friction,Fdes,Zdes,spl),p_mid,[],[],[],[],l0,u0,[],options);
costFunc(p_mid,sole,friction,Fdes,Zdes,spl,coorini)

% fseminf 
%[X,fval,exitflag] = fmincon(@(p_mid)costFunc(p_mid,sole,soleini,friction,Fdes,Zdes,spl),p_mid,[],[],[],[],l0,u0,[],options);
%psoptions = psoptimset('Display','iter');
%X = patternsearch(func,p_mid,[],[],[],[],l0,u0,[],psoptions);
% opts = optimset('Display','iter');
% opts.Display = 'iter';
% opts.TolX = 1.e-12;
% opts.TolFun = 1.e-12;
% opts.MaxFunEvals = 100;
% X=fminsearchcon(@(p_mid)costFunc(p_mid,sole,soleini,friction,Fdes,Zdes,spl),p_mid,l0,u0,[],[],@(p_mid)fminconstr_detJ(p_mid,soleini,spl),opts);
% X=fminsearch(@(p_mid)costFunc(p_mid,sole,soleini,friction,Fdes,Zdes,spl),p_mid,opts);