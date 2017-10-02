f=200;
time=[0:length(trajectories_zmp.xpzmp)-1]'*1/f;
% time=[0:160-1]'*1/f;

w=wpg_param.w;

q=[-13,-3,-w];
k1=-q(1)*q(2)-1;
k2=(q(1)+q(2))/-q(3);
k3=w;

perturb=zeros(length(trajectories_zmp.xpzmp),1);
perturb=[time perturb];
%%
p_d=trajectories_zmp.xpzmp-trajectories_zmp.xpzmp(1);
% p_d=repmat(p_d(1),length(time),1);
p_d=[time p_d];

x_d=trajectories_zmp.xpcom-trajectories_zmp.xpzmp(1);
% x_d=repmat(x_d(1),length(time),1);
x_d=[time x_d];

sx_d=zmp.A_xcom_spd*wpg_param.psa_abcdDSP(1:size(zmp.A_xcom_spd,2))+zmp.B_xcom_spd;
% sx_d=repmat(sx_d(1),length(time),1);
sx_d=[time sx_d];

% options=simset('srcworkspace','current');
% simOut=sim('test stabilizer/simu_control_robot','StartTime','0','StopTime',sprintf('%d',time(end)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f),'srcworkspace','current');
simOut=sim('test stabilizer/simu_control_robot_2','StartTime','0','StopTime',sprintf('%d',time(end)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f),'srcworkspace','current');

x_simu.p_=simOut.get('p_');
x_simu.p=simOut.get('p');
x_simu.x=simOut.get('x');
x_simu.sx=simOut.get('sx');
x_simu.ax=simOut.get('ax');

figure()
hold on
plot(time,[p_d(:,2) x_d(:,2) x_simu.p x_simu.x x_simu.p_])
hold off
% %%
% addpath ./test' stabilizer'/
% figure()
% hold on
% x_simu.p=p_d(1,2);
% x_simu.x=x_d(1,2);
% x_simu.sx=sx_d(1,2);
% for i=2:length(trajectories_zmp.xpzmp)
%     simu=control_layer_zmp_com_init_2(p_d(i-1:i,2),x_d(i-1:i,2),sx_d(i-1:i,2),x_simu.p,x_simu.x,x_simu.sx,wpg_param.frequency,wpg_param.w);
%     plot(i,[p_d(i,2) x_d(i,2) x_simu.p x_simu.x],'*')
% end
% hold off
% rmpath ./test' stabilizer'/
%%
sp=[];
p=[x_simu.p(1)];
for i=1:80
    sp=[sp;(x_simu.p_(i)-x_simu.p(i))/0.05];
    p=[p;p(end)+sp(end)*1/f];
end
figure()
plot(time(1:80),[p_d(1:80,2) x_simu.p_(1:80) x_simu.p(1:80) p(1:80)])
%%
p_d=trajectories_zmp.ypzmp;
p_d=[time p_d];

x_d=trajectories_zmp.ypcom;
x_d=[time x_d];

sx_d=zmp.A_ycom_spd*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_ycom_spd,2))+zmp.B_ycom_spd;
sx_d=[time sx_d];

simOut=sim('test stabilizer/simu_control_robot','StartTime','0','StopTime',sprintf('%d',time(end)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f));
y_simu.p_=simOut.get('p_');
y_simu.p=simOut.get('p');
y_simu.x=simOut.get('x');
y_simu.sx=simOut.get('sx');
y_simu.ax=simOut.get('ax');

figure()
grid on
hold on
plot(time,[p_d(:,2) x_d(:,2) y_simu.p y_simu.x])
hold off


%%
figure()
hold on
plot(sx_d(:,2),'r')
plot((x_d(2:end,2)-x_d(1:end-1,2))*100)
hold off
%%
figure()
hold on
plot(zmp.A_yzmp_spd*psa_abcdDSP(end/2+1:end/2+size(zmp.A_yzmp_spd,2))+zmp.B_yzmp_spd,'r')
plot((p_d(2:end,2)-p_d(1:end-1,2))*100)
hold off

%%
time_r=time(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));

k=diag(zmp.k_diag);
k(sum(wpg_param.discretization(1:wpg_param.nbpolypi))+1:end)=[];
% xp_d=trajectories_zmp.xpzmp(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))-trajectories_zmp.xpzmp(1);
xp_d=trajectories_zmp.xpzmp(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
yp_d=trajectories_zmp.ypzmp(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
xp=x_simu.p(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))+trajectories_zmp.xpzmp(1);
yp=y_simu.p(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
xp_=x_simu.p_(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))+trajectories_zmp.xpzmp(1);
yp_=y_simu.p_(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
% xp1_d=trajectories_zmp.xpzmp1(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))-trajectories_zmp.xpzmp1(1);
% xp2_d=trajectories_zmp.xpzmp2(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))-trajectories_zmp.xpzmp2(1);
xp1_d=trajectories_zmp.xpzmp1(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
xp2_d=trajectories_zmp.xpzmp2(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
yp1_d=trajectories_zmp.ypzmp1(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
yp2_d=trajectories_zmp.ypzmp2(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));

xp1=xp1_d+k.*(xp-xp_d);
xp2=xp2_d+(1-k).*(xp-xp_d);
yp1=yp1_d+k.*(yp-yp_d);
yp2=yp2_d+(1-k).*(xp-xp_d);

xp1_=xp1_d+k.*(xp_-xp_d);
xp2_=xp2_d+(1-k).*(xp_-xp_d);
yp1_=yp1_d+k.*(yp_-yp_d);
yp2_=yp2_d+(1-k).*(xp_-xp_d);

% %%
% figure()
% grid on
% hold on
% plot(time_r,[xp_d xp xp1_d xp1])
% hold off
% 
% figure()
% grid on
% hold on
% plot(time_r,[xp_d xp xp2_d xp2])
% hold off
% 
% figure()
% grid on
% hold on
% plot(time_r,[yp_d yp yp1_d yp1])
% hold off
% 
% figure()
% grid on
% hold on
% plot(time_r,[yp_d yp yp2_d yp2])
% hold off

Lfooty=wpg_param.inttoankle+wpg_param.exttoankle;
traslx=wpg_param.step_number_pankle_fixed(1,2)-wpg_param.backtoankle;
trasly=wpg_param.step_number_pankle_fixed(1,3)-(Lfooty/2-wpg_param.inttoankle);
xp1_d=xp1_d-traslx;
yp1_d=yp1_d-trasly;
xp1=xp1-traslx;
yp1=yp1-trasly;
xp1_=xp1_-traslx;
yp1_=yp1_-trasly;

trasly=wpg_param.step_number_pankle_fixed(2,3)+(Lfooty/2-wpg_param.inttoankle);
xp2_d=xp2_d-traslx;
yp2_d=yp2_d-trasly;
xp2=xp2-traslx;
yp2=yp2-trasly;
xp2_=xp2_-traslx;
yp2_=yp2_-trasly;
%%
Fx_d=-zmp.A_xfcom*wpg_param.psa_abcdDSP(1:size(zmp.A_xfcom,2))+zmp.B_xfcom;
Fy_d=-zmp.A_yfcom*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_yfcom,2))+zmp.B_xfcom;
Fx_d=Fx_d(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
Fy_d=Fy_d(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
Fz_d=wpg_param.mg;

Fx=wpg_param.m*x_simu.ax(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
Fy=wpg_param.m*y_simu.ax(1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)));
Fz=wpg_param.mg;

Fx1_d=k.*Fx_d;
Fy1_d=k.*Fy_d;
Fz1_d=k.*Fz_d;

Fx2_d=(1-k).*Fx_d;
Fy2_d=(1-k).*Fy_d;
Fz2_d=(1-k).*Fz_d;

Fx1=k.*Fx;
Fy1=k.*Fy;
Fz1=k.*Fz;

Fx2=(1-k).*Fx;
Fy2=(1-k).*Fy;
Fz2=(1-k).*Fz;
%%
load('Simulation_3D/results/J.mat');

%%
% angleact1 = [0.;0.;0];
% displ1 = [0;0;0];
angleact1 = Jmat.step_1.angleact_tot(:,1);
displ1 = Jmat.step_1.displ_tot(:,1);
angleact_tot1 = [angleact1];
displ_tot1 = [displ1];

angleact2 = Jmat.step_2.angleact_tot(:,1);
displ2 = Jmat.step_2.displ_tot(:,1);
angleact_tot2 = [angleact2];
displ_tot2 = [displ2];
for i=2:sum(wpg_param.discretization(1:wpg_param.nbpolypi))
    delta=[Fx1_d(i-1);Fy1_d(i-1);Fz1_d(i-1);xp1_d(i-1);yp1_d(i-1);0]-[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0];
    displangleact1=[displ1;angleact1]-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta;
    displ1=displangleact1(1:3,:);
    angleact1=displangleact1(4:6,:);
    
    angleact_tot1 = [angleact_tot1 angleact1];
    displ_tot1 = [displ_tot1 displ1];
    
    delta=[Fx2_d(i-1);Fy2_d(i-1);Fz2_d(i-1);xp2_d(i-1);yp2_d(i-1);0]-[Fx2_d(i);Fy2_d(i);Fz2_d(i);xp2_d(i);yp2_d(i);0];
    displangleact2=[displ2;angleact2]-Jmat.step_2.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta;
    displ2=displangleact2(1:3,:);
    angleact2=displangleact2(4:6,:);
    
    angleact_tot2 = [angleact_tot2 angleact2];
    displ_tot2 = [displ_tot2 displ2];
end
%%
% angleact1 = [0.;0.;0];
% displ1 = [0;0;0];
angleact1 = Jmat.step_1.angleact_tot(:,1);
displ1 = Jmat.step_1.displ_tot(:,1);
angleact_tot1 = [angleact1];
displ_tot1 = [displ1];

angleact2 = Jmat.step_2.angleact_tot(:,1);
displ2 = Jmat.step_2.displ_tot(:,1);
angleact_tot2 = [angleact2];
displ_tot2 = [displ2];
for i=2:sum(wpg_param.discretization(1:wpg_param.nbpolypi))
    delta=[-Fx1(i-1);-Fy1(i-1);Fz1(i-1);xp1(i-1);yp1(i-1);0]-[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0];
    deltadisplangeact=-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta;
    displangleact1=[displ1;angleact1]+deltadisplangeact;
    displ1=displangleact1(1:3,:);
    angleact1=displangleact1(4:6,:);
    
    angleact_tot1 = [angleact_tot1 angleact1];
    displ_tot1 = [displ_tot1 displ1];
    
    delta=[-Fx2(i-1);-Fy2(i-1);Fz2(i-1);xp2(i-1);yp2(i-1);0]-[Fx2_d(i);Fy2_d(i);Fz2_d(i);xp2_d(i);yp2_d(i);0];
    displangleact2=[displ2;angleact2]-Jmat.step_2.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta;
    displ2=displangleact2(1:3,:);
    angleact2=displangleact2(4:6,:);
    
    angleact_tot2 = [angleact_tot2 angleact2];
    displ_tot2 = [displ_tot2 displ2];
end
%%
% angleact1 = [0.;0.;0];
% displ1 = [0;0;0];
angleact1 = Jmat.step_1.angleact_tot(:,1);
displ1 = Jmat.step_1.displ_tot(:,1);
angleact_tot1 = [angleact1];
displ_tot1 = [displ1];

angleact2 = Jmat.step_2.angleact_tot(:,1);
displ2 = Jmat.step_2.displ_tot(:,1);
angleact_tot2 = [angleact2];
displ_tot2 = [displ2];
for i=2:sum(wpg_param.discretization(1:wpg_param.nbpolypi))
    delta=[Fx1_d(i-1);Fy1_d(i-1);Fz1_d(i-1);xp1_d(i-1);yp1_d(i-1);0]-[-Fx1_d(i);-Fy1_d(i);Fz1_d(i);xp1_(i);yp1_(i);0];
    deltadisplangeact=-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta;
    displangleact1=[Jmat.step_1.displ_tot(:,i-1);Jmat.step_1.angleact_tot(:,i-1)]+deltadisplangeact;
    displ1=displangleact1(1:3,:);
    angleact1=displangleact1(4:6,:);
    
    angleact_tot1 = [angleact_tot1 angleact1];
    displ_tot1 = [displ_tot1 displ1];
    
    delta=[Fx2_d(i-1);Fy2_d(i-1);Fz2_d(i-1);xp2_d(i-1);yp2_d(i-1);0]-[-Fx2(i);-Fy2(i);Fz2(i);xp2(i);yp2(i);0];
    displangleact2=[displ2;angleact2]-Jmat.step_2.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta;
    displ2=displangleact2(1:3,:);
    angleact2=displangleact2(4:6,:);
    
    angleact_tot2 = [angleact_tot2 angleact2];
    displ_tot2 = [displ_tot2 displ2];
end
%%
angleact1 = Jmat.step_1.angleact_tot(:,1);
displ1 = Jmat.step_1.displ_tot(:,1);
angleact_tot1 = [angleact1];
displ_tot1 = [displ1];

angleact2 = Jmat.step_2.angleact_tot(:,1);
displ2 = Jmat.step_2.displ_tot(:,1);
angleact_tot2 = [angleact2];
displ_tot2 = [displ2];
for i=2:sum(wpg_param.discretization(1:wpg_param.nbpolypi))
    delta=[Jmat.step_1.Ftot_tot(1,i-1);Jmat.step_1.Ftot_tot(2,i-1);Jmat.step_1.Ftot_tot(3,i-1);Jmat.step_1.Z_tot(1,i-1);Jmat.step_1.Z_tot(2,i-1);Jmat.step_1.Mo_tot(i-1)]-[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0];
    deltadisplangeact=-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta;
%     norm_delta_displ=norm(deltadisplangeact(1:3));
%     norm_delta_angle=norm(deltadisplangeact(4:6));
%     if (norm_delta_displ>0.005)
%         1
%         if (norm_delta_displ*10*3.14159/180>0.005*norm_delta_angle)
%             delta_displangle = delta_displangle*0.005/norm_delta_displ;
% %             //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_displ/max_norm_ddispl);
%         else 
%             delta_displangle = delta_displangle*10*3.14159/180/norm_delta_angle;
% %             //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_angle/max_norm_dangle);
%         end
%     elseif (norm_delta_angle>10*3.14159/180)
%         2
%             delta_displangle = delta_displangle*10*3.14159/180/norm_delta_angle;
% %             //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_angle/max_norm_dangle);
%     end
    displangleact1=[displ1;angleact1]+deltadisplangeact;
    displ1=displangleact1(1:3,:);
    angleact1=displangleact1(4:6,:);
    
    angleact_tot1 = [angleact_tot1 angleact1];
    displ_tot1 = [displ_tot1 displ1];
    
    delta=[Jmat.step_2.Ftot_tot(1,i-1);Jmat.step_2.Ftot_tot(2,i-1);Jmat.step_2.Ftot_tot(3,i-1);Jmat.step_2.Z_tot(1,i-1);Jmat.step_2.Z_tot(2,i-1);Jmat.step_2.Mo_tot(i-1)]-[Fx2_d(i);Fy2_d(i);Fz2_d(i);xp2_d(i);yp2_d(i);0];
    displangleact2=[displ2;angleact2]-Jmat.step_2.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta;
    displ2=displangleact2(1:3,:);
    angleact2=displangleact2(4:6,:);
    
    angleact_tot2 = [angleact_tot2 angleact2];
    displ_tot2 = [displ_tot2 displ2];
end

%%
close all

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.displ_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' displ_tot1(1,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.displ_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' displ_tot1(2,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.displ_tot(3,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' displ_tot1(3,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.angleact_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' angleact_tot1(1,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.angleact_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' angleact_tot1(2,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.angleact_tot(3,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' angleact_tot1(3,:)'],'LineWidth',2)
hold off

%%
close all

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.displ_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' displ_tot2(1,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.displ_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' displ_tot2(2,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.displ_tot(3,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' displ_tot2(3,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.angleact_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' angleact_tot2(1,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.angleact_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' angleact_tot2(2,:)'],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.angleact_tot(3,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' angleact_tot2(3,:)'],'LineWidth',2)
hold off

%%
close all
figure()
grid on
hold on
plot(time_r,[Jmat.step_1.Ftot_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fx1_d],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.Ftot_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fy1_d],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.Ftot_tot(3,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fz1_d],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.Z_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' xp1_d],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_1.Z_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' yp1_d],'LineWidth',2)
hold off

%%
close all
figure()
grid on
hold on
plot(time_r,[Jmat.step_2.Ftot_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fx2_d],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.Ftot_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fy2_d],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.Ftot_tot(3,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fz2_d],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.Z_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' xp2_d],'LineWidth',2)
hold off

figure()
grid on
hold on
plot(time_r,[Jmat.step_2.Z_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' yp2_d],'LineWidth',2)
hold off
%%
close all
flexible_sole_model(angleact_tot1,displ_tot1,0.65,80000*4,0.3,true);
close all
flexible_sole_model(angleact_tot2,displ_tot2,0.65,80000*4,0.3);

flexible_sole_model(Jmat.step_1.angleact_tot,Jmat.step_1.displ_tot,0.65,80000*4,0.3);

flexible_sole_model([],[],0.65,80000*4,0.3);

[sole param_sopt ABe Pg FcFreeSurf D]=flexible_sole_model_init(0.65,80000*4,0.3);
for i=1:size(angleact_tot1,2)
    angleact = angleact_tot1(:,i);
    displ = displ_tot1(:,i);
    contAngle = i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PsoleAngle.m
    [sole Z Ftot param_sopt ABe Pg FcFreeSurf D]=flexible_sole_model_update(sole,angleact,displ,contAngle,param_sopt,ABe,Pg,FcFreeSurf,D,true);
end
% %%
% F_tot1 = [];
% X_tot1 = [];
% for i=2:sum(wpg_param.discretization(1:wpg_param.nbpolypi))
%     delta=[Jmat.step_2.displ_tot(1,i-1);Jmat.step_2.displ_tot(2,i-1);Jmat.step_2.displ_tot(3,i-1);Jmat.step_2.angleact_tot(1,i-1);Jmat.step_2.angleact_tot(2,i-1);Jmat.step_2.angleact_tot(3,i-1)]-[Jmat.step_2.displ_tot(1,i);Jmat.step_2.displ_tot(2,i);Jmat.step_2.displ_tot(3,i);Jmat.step_2.angleact_tot(1,i);Jmat.step_2.angleact_tot(2,i);Jmat.step_2.angleact_tot(3,i)];
%     FX=[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0]+Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)*delta;
%     F=FX(1:3,:);
%     X=FX(4:5,:);
%     
%     F_tot1 = [F_tot1 F];
%     X_tot1 = [X_tot1 X];
% end
% F_tot1 = [F_tot1 [Fx1_d(end);Fy1_d(end);Fz1_d(end)]];
% X_tot1 = [X_tot1 [xp1_d(end);yp1_d(end)]];
% 
% close all
% figure()
% grid on
% hold on
% plot(time_r,[Jmat.step_1.Ftot_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fx1_d F_tot1(1,:)'],'LineWidth',2)
% hold off
% 
% figure()
% grid on
% hold on
% plot(time_r,[Jmat.step_1.Ftot_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fy1_d F_tot1(2,:)'],'LineWidth',2)
% hold off
% 
% figure()
% grid on
% hold on
% plot(time_r,[Jmat.step_1.Ftot_tot(3,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' Fz1_d F_tot1(3,:)'],'LineWidth',2)
% hold off
% 
% figure()
% grid on
% hold on
% plot(time_r,[Jmat.step_1.Z_tot(1,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' xp1_d X_tot1(1,:)'],'LineWidth',2)
% hold off
% 
% figure()
% grid on
% hold on
% plot(time_r,[Jmat.step_1.Z_tot(2,1:sum(wpg_param.discretization(1:wpg_param.nbpolypi)))' yp1_d X_tot1(2,:)'],'LineWidth',2)
% hold off