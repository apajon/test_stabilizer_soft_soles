% duration=0.8;%in second
duration=1.2;%in second
f=wpg_param.frequency;
duration_disc=duration*f;
time=[0:duration_disc-1]'*1/f;


w=wpg_param.w;

q=[-13,-3,-w];
k1=-q(1)*q(2)-1;
k2=(q(1)+q(2))/-q(3);
k3=w;

perturb=zeros(duration_disc,1);
perturb=[time perturb];
%%
time_r=time;


% k=diag(zmp.k_diag);
% k(sum(wpg_param.discretization(1:wpg_param.nbpolypi))+1:end)=[];
xp_d=repmat(trajectories_zmp.xpzmp(1),duration_disc,1);
yp_d=repmat(trajectories_zmp.ypzmp(1),duration_disc,1);

xp1_d=repmat(trajectories_zmp.xpzmp1(1),duration_disc,1);
yp1_d=repmat(trajectories_zmp.ypzmp1(1),duration_disc,1);

xp2_d=repmat(trajectories_zmp.xpzmp2(1),duration_disc,1);
yp2_d=repmat(trajectories_zmp.ypzmp2(1),duration_disc,1);

Lfooty=wpg_param.inttoankle+wpg_param.exttoankle;
traslx_1=wpg_param.step_number_pankle_fixed(1,2)-wpg_param.backtoankle;
trasly_1=wpg_param.step_number_pankle_fixed(1,3)-(Lfooty/2-wpg_param.inttoankle);
xp1_d_sole=xp1_d-traslx_1;
yp1_d_sole=yp1_d-trasly_1;

traslx_2=wpg_param.step_number_pankle_fixed(1,2)-wpg_param.backtoankle;
trasly_2=wpg_param.step_number_pankle_fixed(2,3)+(Lfooty/2-wpg_param.inttoankle);
xp2_d_sole=xp2_d-traslx_2;
yp2_d_sole=yp2_d-trasly_2;
%%
Fx_d=repmat(-zmp.A_xfcom(1,:)*wpg_param.psa_abcdDSP(1:size(zmp.A_xfcom,2))-zmp.B_xfcom(1),duration_disc,1);
Fy_d=repmat(-zmp.A_yfcom(1,:)*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_yfcom,2))-zmp.B_yfcom(1),duration_disc,1);
Fz_d=repmat(wpg_param.mg,duration_disc,1);

Fx1_d=0.5.*Fx_d;
Fy1_d=0.5.*Fy_d;
Fz1_d=0.5.*Fz_d;

Fx2_d=0.5.*Fx_d;
Fy2_d=0.5.*Fy_d;
Fz2_d=0.5.*Fz_d;
%%
xx_d=repmat(trajectories_zmp.xpcom(1),duration_disc,1);
yx_d=repmat(trajectories_zmp.ypcom(1),duration_disc,1);

xsx_d=repmat(zmp.A_xcom_spd(1,:)*wpg_param.psa_abcdDSP(1:size(zmp.A_xcom_spd,2))+zmp.B_xcom_spd(1),duration_disc,1);
ysx_d=repmat(zmp.A_ycom_spd(1,:)*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_xcom_spd,2))+zmp.B_ycom_spd(1),duration_disc,1);
%%
xp=[xp_d(1)];
xp1=[xp1_d(1)];
xp2=[xp2_d(1)];

yp=[yp_d(1)];
yp1=[yp1_d(1)];
yp2=[yp2_d(1)];

xp_=[xp_d(1)];
xp1_=[xp1_d(1)];
xp2_=[xp2_d(1)];

yp_=[yp_d(1)];
yp1_=[yp1_d(1)];
yp2_=[yp2_d(1)];

xp1_sole=[xp1_d_sole(1)];
xp2_sole=[xp2_d_sole(1)];

yp1_sole=[yp1_d_sole(1)];
yp2_sole=[yp2_d_sole(1)];

Fx=[Fx_d(1)];
Fy=[Fy_d(1)];
Fz=[Fz_d(1)];

Fx_=[Fx_d(1)];
Fy_=[Fy_d(1)];
Fz_=[Fz_d(1)];

Fx1=[Fx1_d(1)];
Fy1=[Fy1_d(1)];
Fz1=[Fz1_d(1)];

Fx2=[Fx2_d(1)];
Fy2=[Fy2_d(1)];
Fz2=[Fz2_d(1)];

xx=[xx_d(1)];
yx=[yx_d(1)];

xsx=[xsx_d(1)];
ysx=[ysx_d(1)];

xax=[0];
yax=[0];

xp_=[xp_d(1)];
xp1_=[xp1_d(1)];
xp2_=[xp2_d(1)];

yp_=[yp_d(1)];
yp1_=[yp1_d(1)];
yp2_=[yp2_d(1)];

xx_=[xx_d(1)];
yx_=[yx_d(1)];

xsx_=[xsx_d(1)];
ysx_=[ysx_d(1)];

Fx_=[Fx_d(1)];
Fy_=[Fy_d(1)];
Fz_=[Fz_d];

Fx1_=[Fx1_d(1)];
Fy1_=[Fy1_d(1)];
Fz1_=[Fz1_d(1)];

Fx2_=[Fx2_d(1)];
Fy2_=[Fy2_d(1)];
Fz2_=[Fz2_d(1)];
%%
load('Simulation_3D/results/J_0.mat');
angleact1 = Jmat.step_1.angleact_tot(:,1);
displ1 = Jmat.step_1.displ_tot(:,1);
angleact_tot1 = [angleact1];
displ_tot1 = [displ1];

angleact_tot1_ = [angleact1];
displ_tot1_ = [displ1];

angleact2 = Jmat.step_2.angleact_tot(:,1);
displ2 = Jmat.step_2.displ_tot(:,1);
angleact_tot2 = [angleact2];
displ_tot2 = [displ2];
%%

xp_d_init=xp_d(1);
xp_d=[time xp_d-xp_d_init];
yp_d_init=yp_d(1);
yp_d=[time yp_d-yp_d_init];

xp1_d=[time xp1_d];
yp1_d=[time yp1_d];

xp2_d=[time xp2_d];
yp2_d=[time yp2_d];

xx_d_init=xx_d(1);
xx_d=[time xx_d-xx_d_init];
yx_d_init=yx_d(1);
yx_d=[time yx_d-yx_d_init];

xsx_d_init=xsx_d(1);
xsx_d=[time xsx_d-xsx_d_init];
ysx_d_init=ysx_d(1);
ysx_d=[time ysx_d-ysx_d_init];

Fx1_d=[time Fx1_d];
Fy1_d=[time Fy1_d];
Fz1_d=[time Fz1_d];

Fx2_d=[time Fx2_d];
Fy2_d=[time Fy2_d];
Fz2_d=[time Fz2_d];

J1=[];
J1_inv=[];
for i=2:80
    J1=[J1;Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)];
    J1_inv=[J1_inv;J1(end-5:end,:)^-1];
end
%%
% addpath ./test' stabilizer'/
[sole_1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_init(0.65,80000*4,0.3);
[sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,1,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);

[sole_2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_init(0.65,80000*4,0.3);
[sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,angleact2,displ2,1,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);

xp1=[time(1) xp1];
yp1=[time(1) yp1];
Fx1=[time(1) Fx1];
Fy1=[time(1) Fy1];
Fz1=[time(1) Fz1];

%%
for i=2:duration_disc
    i
    %%
    simOut=sim('test stabilizer/simu_control_robot_4','StartTime',sprintf('%d',time(2)),'StopTime',sprintf('%d',time(i)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f),'srcworkspace','current');
    displangleact=simOut.get('displangleact');
    displ1=displangleact.Data(1:3,:,end);
    angleact1=displangleact.Data(4:6,:,end);
    
    [sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,1,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);
    Fx1=[Fx1;time(i) Ftot1(1,end)];
    Fy1=[Fy1;time(i) Ftot1(2,end)];
    Fz1=[Fz1;time(i) Ftot1(3,end)];
    xp1_sole=[xp1_sole;Z1(1,end)];
    yp1_sole=[yp1_sole;Z1(2,end)];
    
    xp1=[xp1;time(i) xp1_sole(end)+traslx_1];
    yp1=[yp1;time(i) yp1_sole(end)+traslx_1];
end