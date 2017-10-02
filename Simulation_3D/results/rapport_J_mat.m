clear all
close all
clc
format long

load('Simulation_3D/results/J_f10.mat');
Jmat_f10=Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1;

load('Simulation_3D/results/J_f50.mat');
Jmat_f50=Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1;

load('Simulation_3D/results/J_f100.mat');
Jmat_f100=Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1;

load('Simulation_3D/results/J_f210.mat');
Jmat_f200=Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1;

load('Simulation_3D/results/J_f400.mat');
Jmat_f400=Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1;

load('Simulation_3D/results/J_f600.mat');
Jmat_f600=Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1;

clear Jmat
%%
Jmat_f50(4:6,4:6)./Jmat_f10(4:6,4:6)
Jmat_f100(4:6,4:6)./Jmat_f50(4:6,4:6)
Jmat_f200(4:6,4:6)./Jmat_f100(4:6,4:6)
Jmat_f400(4:6,4:6)./Jmat_f200(4:6,4:6)
Jmat_f600(4:6,4:6)./Jmat_f400(4:6,4:6)


%%
Jmat_f10(4:6,4:6)./Jmat_f200(4:6,4:6)
Jmat_f50(4:6,4:6)./Jmat_f200(4:6,4:6)
Jmat_f100(4:6,4:6)./Jmat_f200(4:6,4:6)
Jmat_f400(4:6,4:6)./Jmat_f200(4:6,4:6)
Jmat_f600(4:6,4:6)./Jmat_f200(4:6,4:6)

Jmat_f10(3,3)./Jmat_f200(3,3)
Jmat_f50(3,3)./Jmat_f200(3,3)
Jmat_f100(3,3)./Jmat_f200(3,3)
Jmat_f400(3,3)./Jmat_f200(3,3)
Jmat_f600(3,3)./Jmat_f200(3,3)

%%
figure()
hold on
plot([10 50 100 200 400 600],[0.058 0.29 0.58 1 2.31 3.45])
hold off

%%
figure()
hold on
plot([20 60 220 360 500 1670],[0.09 0.06 0.04 0.03 0.025 0.007])
hold off
