function [x_simu y_simu]=control_layer_zmp_com_init(wpg_param,zmp,trajectories_zmp,f,w)

% f=wpg_param.frequency;
time=[0:length(trajectories_zmp.xpzmp)-1]'*1/f;

% w=wpg_param.w;

q=[-13,-3,-w];
k1=-q(1)*q(2)-1;
k2=(q(1)+q(2))/-q(3);
k3=w;

perturb=zeros(length(trajectories_zmp.xpzmp),1);
perturb=[time perturb];
%%
p_d_x=trajectories_zmp.xpzmp-trajectories_zmp.xpzmp(1);
p_d=[time p_d_x];

x_d_x=trajectories_zmp.xpcom-trajectories_zmp.xpzmp(1);
x_d_x=[time x_d_x];

sx_d_x=zmp.A_xcom_spd*wpg_param.psa_abcdDSP(1:size(zmp.A_xcom_spd,2))+zmp.B_xcom_spd;
sx_d=[time sx_d_x];

simset('srcworkspace','current')
simOut=sim('test stabilizer/simu_control_robot','StartTime','0','StopTime',sprintf('%d',time(end)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f));
x_simu.p_=simOut.get('p_');
x_simu.p=simOut.get('p');
x_simu.x=simOut.get('x');
x_simu.sx=simOut.get('sx');
x_simu.ax=simOut.get('ax');

%%
p_d_y=trajectories_zmp.ypzmp;
p_d=[time p_d_y];
% p_d=[p_d p_d_y];

x_d_y=trajectories_zmp.ypcom;
x_d=[time x_d_y];
% x_d=[x_d x_d_y];

sx_d_y=zmp.A_ycom_spd*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_ycom_spd,2))+zmp.B_ycom_spd;
sx_d=[time sx_d_y];
% sx_d=[sx_d sx_d_y];

simOut=sim('test stabilizer/simu_control_robot','StartTime','0','StopTime',sprintf('%d',time(end)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f));
y_simu.p_=simOut.get('p_');
y_simu.p=simOut.get('p');
y_simu.x=simOut.get('x');
y_simu.sx=simOut.get('sx');
y_simu.ax=simOut.get('ax');

% %%
% figure()
% hold on
% plot(sx_d(:,2),'r')
% plot((x_d(2:end,2)-x_d(1:end-1,2))*100)
% hold off
% %%
% figure()
% hold on
% plot(zmp.A_yzmp_spd*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_yzmp_spd,2))+zmp.B_yzmp_spd,'r')
% plot((p_d(2:end,2)-p_d(1:end-1,2))*100)
% hold off

end