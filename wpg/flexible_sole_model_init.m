function [sole param_sopt ABe Pg FcFreeSurf D]=flexible_sole_model_init(friction,Young,Poisson)

format long
% clc
% clear all
% close all

addpath ./Model_3D/
% addpath ./Model_3D/input
addpath ./Model_3D/FEM
addpath ./Model_3D/gradient

clear sole soleini
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Sole FEM                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fname = 'semelle1.msh';
% pname = 'Model_3D/input/semelle1 L=0.23, l=0.13, e=0.03 m new centre/';
pname = 'Simulation_3D/input/semelle1 L=0.23, l=0.13, e=0.021 m new centre/';
% pname = 'Simulation_3D/input/semelle1 L=0.23, l=0.13, e=0.021 m new centre_low/';
% pname = 'Simulation_3D/input/semelle1 L=0.23, l=0.13, e=0.015 m new centre/';
% pname = 'Simulation_3D/input/semelle1 L=0.23, l=0.13, e=0.006 m new centre/';

%%% foot size %%%
l=0.13;
L=0.23;
% e=0.03;
e=0.021;
% e=0.006;
%%
%%% sole FEM creation %%%
sole = soleFEM_newStiff(pname,fname,l,L,e);
if isempty(Young)
    Young = 1000000;
end
if isempty(Poisson)
    Poisson = 0.3;
end
sole.setMaterial(Young,Poisson);
% Fdes = 1.0e+02 * [-0.052709310000000;0.030650450000000;1.279851840000000];
% Zdes = [0.020625645823457;0.049298000000611];

%%
%%% Create shape parameters struct
move_dirichlet = 1;
param_sopt = struct;
if isempty(friction)
    friction=0.8;
end
param_sopt.friction = friction;
% param_sopt.Fdes = Fdes;
% param_sopt.Zdes = Zdes;
param_sopt.move_dirichlet = move_dirichlet;

sole.stiffness();
sole.stiffnessSurface();
sole.derStiff();

ABe = prep_stressVonMises(sole);
%%% Test Gradient algo with given position and orientation
Pg = [];
FcFreeSurf = [];
D = 3;

rmpath ./Model_3D/
% rmpath ./Model_3D/input
rmpath ./Model_3D/FEM
rmpath ./Model_3D/gradient