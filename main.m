clear all
close all
clc
format long

% addpath ./Simulation_3D
% addpath ./Simulation_3D/input
% addpath ./Simulation_3D/results
% addpath ./Simulation_3D/FEM
% addpath ./Simulation_3D/iso2mesh
% addpath ./Simulation_3D/gradient
addpath ./wpg
% addpath ./fsqp
% addpath ./input
%%
%%%%%%%%%%%%%%%%%%%%%%%%%% WPG parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% create walking parameters
robot = 2;
type_traj=1;
firstSS = -1;
frequency = 200;%% WARNING : all phase duration (tss, tds) * frequency / nbpoly(ssp ou dsp) = positive integer
lambda = 0.9;
epsilon = 100;
% e=0.021;
e = 0.03;
rightorleft = -1; % Chose foot (left or right),+1 for right foot and -1 for left foot
poly_degree=5;
nbpolyssp=6; %% WARNING : must be an even integer
nbpolydsp=8;
nbpolypi=2;
nbpolypf=2;
wpg_param = wpg_parameters(robot,type_traj,firstSS,frequency,lambda,epsilon,e,rightorleft,poly_degree,nbpolyssp,nbpolydsp,nbpolypi,nbpolypf);
clear robot type_traj firstSS frequency lambda epsilon e rightorleft poly_degree nbpolyssp nbpolydsp nbpolypi nbpolypf
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of wpg parameters using the WPG of ICAR 2015 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate QP problem
zmp = wpg_qp_problem_Adrien(wpg_param);
scalf = 1/(1e+6);
zmp.H = zmp.H*scalf;
zmp.f = zmp.f*scalf;
zmp.C = zmp.C*scalf;
clear scalf
% zmp.H = zmp.H/wpg_param.nbpointdiscret;
% zmp.f = zmp.f/wpg_param.nbpointdiscret;
% zmp.C = zmp.C/wpg_param.nbpointdiscret;
% cond(zmp.H)
%zmp.reduce_constraint(1);
% zmp.extend_constraint_symetry(wpg_param);

opt = optimset('Algorithm', 'interior-point-convex','Display','iter');
% opt = optimset('Algorithm', 'interior-point-convex','Display','iter','MaxIter',500);
psa_abcdDSP = quadprog(zmp.H,zmp.f,sparse(zmp.A),sparse(zmp.b),sparse(zmp.Aeq),sparse(-zmp.beq),[],[],[],opt);
clear opt

wpg_param.psa_abcdDSP = psa_abcdDSP;
clear psa_abcdDSP
wpg_param.nbcontrolpointzmp = wpg_param.nbcontrolpointzmp+4;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trajectories_zmp = wpg_trajectories();
trajectories_zmp.drawing(wpg_param,zmp,1,1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute best ankle movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flexible_sole_optim_ankle(wpg_param,zmp,trajectories_zmp,1,0.65,80000*4,0.3)
% %Compute ZMP and Fzmp from ankle angle and position
% flexible_sole_model([],[],[],[],[]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulate the ZMP/COM error and its compensation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath ./test' stabilizer'/
% [x_simu y_simu]=control_layer_zmp_com_init(wpg_param,zmp,trajectories_zmp);
% rmpath ./test' stabilizer'/