function [simu]=control_layer_zmp_com_init_2(ZMP,COM,S_COM,f,w)

% f=wpg_param.frequency;
time=[0:length(ZMP)-1]'*1/f;

% w=wpg_param.w;

q=[-13,-3,-w];
k1=-q(1)*q(2)-1;
k2=(q(1)+q(2))/-q(3);
k3=w;

perturb=zeros(length(ZMP),1);
perturb=[time perturb];
%%
p_d_x=ZMP-ZMP(1);
p_d=[time p_d_x];

x_d_x=COM-ZMP(1);
x_d=[time x_d_x];

sx_d_x=S_COM;
sx_d=[time sx_d_x];

% simset('srcworkspace','current');
simOut=sim('test stabilizer/simu_control_robot','StartTime','0','StopTime',sprintf('%d',time(end)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f),'srcworkspace','current');
simu.p_=simOut.get('p_')+ZMP(1);
simu.p=simOut.get('p')+ZMP(1);
simu.x=simOut.get('x')+ZMP(1);
simu.sx=simOut.get('sx');
simu.ax=simOut.get('ax');

end