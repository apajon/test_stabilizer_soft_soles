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
% %%
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
% %%
Fx_d=repmat(-zmp.A_xfcom(1,:)*wpg_param.psa_abcdDSP(1:size(zmp.A_xfcom,2))-zmp.B_xfcom(1),duration_disc,1);
Fy_d=repmat(-zmp.A_yfcom(1,:)*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_yfcom,2))-zmp.B_yfcom(1),duration_disc,1);
Fz_d=repmat(wpg_param.mg,duration_disc,1);

Fx1_d=0.5.*Fx_d;
Fy1_d=0.5.*Fy_d;
Fz1_d=0.5.*Fz_d;

Fx2_d=0.5.*Fx_d;
Fy2_d=0.5.*Fy_d;
Fz2_d=0.5.*Fz_d;
% %%
xx_d=repmat(trajectories_zmp.xpcom(1),duration_disc,1);
yx_d=repmat(trajectories_zmp.ypcom(1),duration_disc,1);

xsx_d=repmat(zmp.A_xcom_spd(1,:)*wpg_param.psa_abcdDSP(1:size(zmp.A_xcom_spd,2))+zmp.B_xcom_spd(1),duration_disc,1);
ysx_d=repmat(zmp.A_ycom_spd(1,:)*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_xcom_spd,2))+zmp.B_ycom_spd(1),duration_disc,1);
% %%
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
% %%
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
clear sole_1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1 sole_2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2
% addpath ./test' stabilizer'/
[sole_1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_init(0.65,80000*4,0.3);
[sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,1,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);

[sole_2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_init(0.65,80000*4,0.3);
[sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,angleact2,displ2,1,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);
%% no noise
for i=2:duration_disc
    i
    %%
    delta_F1=[Fx1(end);Fy1(end);Fz1(end);xp1(end);yp1(end);0]-[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0];
    displangleact1=[Jmat.step_1.displ_tot(:,i-1);Jmat.step_1.angleact_tot(:,i-1)]-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F1;
    displ1=displangleact1(1:3,:);
    angleact1=displangleact1(4:6,:);
    
    angleact_tot1 = [angleact_tot1 angleact1];
    displ_tot1 = [displ_tot1 displ1];
    
    [sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,i,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);
    Fx1=[Fx1;Ftot1(1,end)];
    Fy1=[Fy1;Ftot1(2,end)];
    Fz1=[Fz1;Ftot1(3,end)];
    xp1_sole=[xp1_sole;Z1(1,end)];
    yp1_sole=[yp1_sole;Z1(2,end)];
    
    delta_F2=[Fx2(end);Fy2(end);Fz2(end);xp2(end);yp2(end);0]-[Fx2_d(i);Fy2_d(i);Fz2_d(i);xp2_d(i);yp2_d(i);0];
    displangleact2=[Jmat.step_2.displ_tot(:,i-1);Jmat.step_2.angleact_tot(:,i-1)]-Jmat.step_2.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F2;
    displ2=displangleact2(1:3,:);
    angleact2=displangleact2(4:6,:);
    
    angleact_tot2 = [angleact_tot2 angleact2];
    displ_tot2 = [displ_tot2 displ2];
    
    [sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,angleact2,displ2,i,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);
    Fx2=[Fx2;Ftot2(1,end)];
    Fy2=[Fy2;Ftot2(2,end)];
    Fz2=[Fz2;Ftot2(3,end)];
    xp2_sole=[xp2_sole;Z2(1,end)];
    yp2_sole=[yp2_sole;Z2(2,end)];
    
    %%
    xp1=xp1_sole+traslx_1;
    yp1=yp1_sole+trasly_1;
    
    xp2=xp2_sole+traslx_2;
    yp2=yp2_sole+trasly_2;
    
    Fx=[Fx;Fx1(end)+Fx2(end)];
    Fy=[Fy;Fy1(end)+Fy2(end)];
    Fz=[Fz;Fz1(end)+Fz2(end)];
    
    xp=[xp; (xp2(end)+Fz1(end)/Fz2(end)*xp1(end))/(1+Fz1(end)/Fz2(end))];
    yp=[yp; (yp2(end)+Fz1(end)/Fz2(end)*yp1(end))/(1+Fz1(end)/Fz2(end))];
    
    xax=[xax;-Fx(end)/wpg_param.m];
    yax=[yax;-Fy(end)/wpg_param.m];
    
    toto=cumtrapz(xax);
    xsx=[xsx;xsx(1)+toto(end)/f];
    toto=cumtrapz(yax);
    ysx=[ysx;ysx(1)+toto(end)/f];
    
    toto=cumtrapz(xsx);
    xx=[xx;xx(1)+toto(end)/f];
    toto=cumtrapz(ysx);
    yx=[yx;yx(1)+toto(end)/f];
    
%     xx=[xx;xp(end)+1/w^2*(-Fx(end))/wpg_param.m];
%     yx=[yx;yp(end)+1/w^2*(-Fy(end))/wpg_param.m];
%     
% %     xsx=[xsx;(xx(end)-xx(end-1))*f];
% %     ysx=[ysx;(yx(end)-yx(end-1))*f];
%     xsx=[xsx;xsx(end)+(-Fx(end))/wpg_param.m/f];
%     ysx=[ysx;ysx(end)+(-Fy(end))/wpg_param.m/f];
end
% rmpath ./test' stabilizer'/
%% with perturbation
for i=2:duration_disc
    i
    %%
%     delta_xp=xp_d(i)-xp(end);
%     delta_xx=xx_d(i)-xx(end);
%     delta_xsx=xsx_d(i)-xsx(end);
%     delta_corr_x=k1*delta_xx+k2*delta_xsx+k3*delta_xp;
%     xp_=[xp_;xp_d(i)+delta_corr_x];
%     xp1_=[xp1_;xp1_d(i)+delta_corr_x];
%     xp2_=[xp2_;xp2_d(i)+delta_corr_x];
%     Fx_=[Fx_;-(xx(end)-xp_(end))*w^2*wpg_param.m];
%     Fx1_=[Fx1_;k(i)*Fx_(end)];
%     Fx2_=[Fx2_;(1-k(i))*Fx_(end)];
%     
%     delta_yp=yp_d(i)-yp(end);
%     delta_yx=yx_d(i)-yx(end);
%     delta_ysx=ysx_d(i)-ysx(end);
%     delta_corr_y=k1*delta_yx+k2*delta_ysx+k3*delta_yp;
%     yp_=[yp_;yp_d(i)+delta_corr_y];
%     yp1_=[yp1_;yp1_d(i)+delta_corr_y];
%     yp2_=[yp2_;yp2_d(i)+delta_corr_y];
%     Fy_=[Fy_;-(yx(end)-yp_(end))*w^2*wpg_param.m];
%     Fy1_=[Fy1_;k(i)*Fy_(end)];
%     Fy2_=[Fy2_;(1-k(i))*Fy_(end)];
    %%
%     if i>=40 && i<61
%         Fx1_d(i)=Fx1_d(i)-10;
%         xp1_d(i)=xp1_d(i)-0.04;
%     end
    delta_F1=[Fx1(end);Fy1(end);Fz1(end);xp1(end);yp1(end);0]-[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0];
%     if i>=80
%         displangleact1=[Jmat.step_1.displ_tot(:,80);Jmat.step_1.angleact_tot(:,80)]-Jmat.step_1.J_tot((80-1)*6+1:(80-1)*6+6,:)^-1*delta_F1;
%     else
%         displangleact1=[Jmat.step_1.displ_tot(:,i-1);Jmat.step_1.angleact_tot(:,i-1)]-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F1;
%         Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F1
%     end
    if i>=80
        displangleact1=[displ_tot1(:,end);angleact_tot1(:,end)]-Jmat.step_1.J_tot((80-1)*6+1:(80-1)*6+6,:)^-1*delta_F1;
    else
        displangleact1=[displ_tot1(:,end);angleact_tot1(:,end)]-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F1;
        Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F1;
    end
    displ1=displangleact1(1:3,:);
    angleact1=displangleact1(4:6,:);
    
    angleact_tot1 = [angleact_tot1 angleact1];
    displ_tot1 = [displ_tot1 displ1];

    [sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,i,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);
    Fx1=[Fx1;Ftot1(1,end)];
    Fy1=[Fy1;Ftot1(2,end)];
    Fz1=[Fz1;Ftot1(3,end)];
    xp1_sole=[xp1_sole;Z1(1,end)];
    yp1_sole=[yp1_sole;Z1(2,end)];
    
%     if i>=40 && i<61
%         Fx2_d(i)=Fx2_d(i)-10;
%         xp2_d(i)=xp2_d(i)-0.04;
%     end
    delta_F2=[Fx2(end);Fy2(end);Fz2(end);xp2(end);yp2(end);0]-[Fx2_d(i);Fy2_d(i);Fz2_d(i);xp2_d(i);yp2_d(i);0];
%     if i>=80
%         displangleact2=[Jmat.step_2.displ_tot(:,80-1);Jmat.step_2.angleact_tot(:,80-1)]-Jmat.step_2.J_tot((80-1)*6+1:(80-1)*6+6,:)^-1*delta_F2;
%     else
%         displangleact2=[Jmat.step_2.displ_tot(:,i-1);Jmat.step_2.angleact_tot(:,i-1)]-Jmat.step_2.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F2;
%     end
    if i>=80
        displangleact2=[displ_tot2(:,end);angleact_tot2(:,end)]-Jmat.step_1.J_tot((80-1)*6+1:(80-1)*6+6,:)^-1*delta_F1;
    else
        displangleact2=[displ_tot2(:,end);angleact_tot2(:,end)]-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F1;
    end
    displ2=displangleact2(1:3,:);
    angleact2=displangleact2(4:6,:);
    
    angleact_tot2 = [angleact_tot2 angleact2];
    displ_tot2 = [displ_tot2 displ2];
    
    [sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,angleact2,displ2,i,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);
    Fx2=[Fx2;Ftot2(1,end)];
    Fy2=[Fy2;Ftot2(2,end)];
    Fz2=[Fz2;Ftot2(3,end)];
    xp2_sole=[xp2_sole;Z2(1,end)];
    yp2_sole=[yp2_sole;Z2(2,end)];
    
    %%
    xp1=xp1_sole+traslx_1;
    yp1=yp1_sole+trasly_1;
    
    xp2=xp2_sole+traslx_2;
    yp2=yp2_sole+trasly_2;
    
    Fx=[Fx;Fx1(end)+Fx2(end)];
    Fy=[Fy;Fy1(end)+Fy2(end)];
    Fz=[Fz;Fz1(end)+Fz2(end)];
    
    xp=[xp; (xp2(end)+Fz1(end)/Fz2(end)*xp1(end))/(1+Fz1(end)/Fz2(end))];
    yp=[yp; (yp2(end)+Fz1(end)/Fz2(end)*yp1(end))/(1+Fz1(end)/Fz2(end))];
    
%     xax=[xax;-Fx(end)/wpg_param.m];
%     yax=[yax;-Fy(end)/wpg_param.m];
%     
%     toto=cumtrapz(xax);
%     xsx=[xsx;xsx(1)+toto(end)/f];
%     toto=cumtrapz(yax);
%     ysx=[ysx;ysx(1)+toto(end)/f];
%     
%     toto=cumtrapz(xsx);
%     xx=[xx;xx(1)+toto(end)/f];
%     toto=cumtrapz(ysx);
%     yx=[yx;yx(1)+toto(end)/f];
    
%     xx=[xx;xp(end)+1/w^2*(-Fx(end))/wpg_param.m];
%     yx=[yx;yp(end)+1/w^2*(-Fy(end))/wpg_param.m];
%     
% %     xsx=[xsx;(xx(end)-xx(end-1))*f];
% %     ysx=[ysx;(yx(end)-yx(end-1))*f];
%     xsx=[xsx;xsx(end)+(-Fx(end))/wpg_param.m/f];
%     ysx=[ysx;ysx(end)+(-Fy(end))/wpg_param.m/f];
end
% rmpath ./test' stabilizer'/
%% with noise
for i=2:sum(wpg_param.discretization(1:wpg_param.nbpolypi))
    i
    %%
    delta_F1=[Fx1(end);Fy1(end);Fz1(end);xp1(end);yp1(end);0]-[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0];
    displangleact1=[Jmat.step_1.displ_tot(:,i-1);Jmat.step_1.angleact_tot(:,i-1)]-Jmat.step_1.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F1;
    displ1=displangleact1(1:3,:);
    angleact1=displangleact1(4:6,:);
    
    angleact_tot1 = [angleact_tot1 angleact1];
    displ_tot1 = [displ_tot1 displ1];
    
    [sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,i,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);
    Fx1=[Fx1;awgn(Ftot1(1,end),20,'measured')];
    Fy1=[Fy1;awgn(Ftot1(2,end),20,'measured')];
    Fz1=[Fz1;awgn(Ftot1(3,end),20,'measured')];
    xp1_sole=[xp1_sole;awgn(Z1(1,end),20,'measured')];
    yp1_sole=[yp1_sole;awgn(Z1(2,end),20,'measured')];
    
    delta_F2=[Fx2(end);Fy2(end);Fz2(end);xp2(end);yp2(end);0]-[Fx2_d(i);Fy2_d(i);Fz2_d(i);xp2_d(i);yp2_d(i);0];
    displangleact2=[Jmat.step_2.displ_tot(:,i-1);Jmat.step_2.angleact_tot(:,i-1)]-Jmat.step_2.J_tot((i-1)*6+1:(i-1)*6+6,:)^-1*delta_F2;
    displ2=displangleact2(1:3,:);
    angleact2=displangleact2(4:6,:);
    
    angleact_tot2 = [angleact_tot2 angleact2];
    displ_tot2 = [displ_tot2 displ2];
    
    [sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,angleact2,displ2,i,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);
    Fx2=[Fx2;awgn(Ftot2(1,end),20,'measured')];
    Fy2=[Fy2;awgn(Ftot2(2,end),20,'measured')];
    Fz2=[Fz2;awgn(Ftot2(3,end),20,'measured')];
    xp2_sole=[xp2_sole;awgn(Z2(1,end),20,'measured')];
    yp2_sole=[yp2_sole;awgn(Z2(2,end),20,'measured')];
    
    %%
    xp1=xp1_sole+traslx_1;
    yp1=yp1_sole+trasly_1;
    
    xp2=xp2_sole+traslx_2;
    yp2=yp2_sole+trasly_2;
    
    Fx=[Fx;Fx1(end)+Fx2(end)];
    Fy=[Fy;Fy1(end)+Fy2(end)];
    Fz=[Fz;Fz1(end)+Fz2(end)];
    
    xp=[xp; (xp2(end)+Fz1(end)/Fz2(end)*xp1(end))/(1+Fz1(end)/Fz2(end))];
    yp=[yp; (yp2(end)+Fz1(end)/Fz2(end)*yp1(end))/(1+Fz1(end)/Fz2(end))];
    
    xax=[xax;-Fx(end)/wpg_param.m];
    yax=[yax;-Fy(end)/wpg_param.m];
    
    toto=cumtrapz(xax);
    xsx=[xsx;xsx(1)+toto(end)/f];
    toto=cumtrapz(yax);
    ysx=[ysx;ysx(1)+toto(end)/f];
    
    toto=cumtrapz(xsx);
    xx=[xx;xx(1)+toto(end)/f];
    toto=cumtrapz(ysx);
    yx=[yx;yx(1)+toto(end)/f];
    
%     xx=[xx;xp(end)+1/w^2*(-Fx(end))/wpg_param.m];
%     yx=[yx;yp(end)+1/w^2*(-Fy(end))/wpg_param.m];
%     
% %     xsx=[xsx;(xx(end)-xx(end-1))*f];
% %     ysx=[ysx;(yx(end)-yx(end-1))*f];
%     xsx=[xsx;xsx(end)+(-Fx(end))/wpg_param.m/f];
%     ysx=[ysx;ysx(end)+(-Fy(end))/wpg_param.m/f];
end
%%
close all
figure()
hold on
plot(time_r,[xp1_d xp1])
hold off

figure()
hold on
plot(time_r,[yp1_d yp1])
hold off

close all
figure()
hold on
plot(time_r,[xp2_d xp2])
hold off

figure()
hold on
plot(time_r,[yp2_d yp2])
hold off

close all
figure()
hold on
plot(time_r,[xp_d xp])
hold off

figure()
hold on
plot(time_r,[yp_d yp])
hold off

close all
figure()
hold on
plot(time_r,[Fx_d Fx ])
hold off

figure()
hold on
plot(time_r,[Fy_d Fy ])
hold off

close all
figure()
hold on
plot(time_r,[Fx1_d Fx1 ])
hold off

figure()
hold on
plot(time_r,[Fy1_d Fy1 ])
hold off

figure()
hold on
plot(time_r,[Fx2_d Fx2 ])
hold off

figure()
hold on
plot(time_r,[Fy2_d Fy2 ])
hold off

close all
figure()
hold on
plot(time_r,[Fx1_d+Fx2_d Fx ])
hold off

figure()
hold on
plot(time_r,[Fy1_d+Fy2_d Fy ])
hold off

figure()
hold on
plot(time_r,[Fz1_d+Fz2_d Fz])
hold off

figure()
hold on
plot(time_r,[xp_d xp])
hold off

figure()
hold on
plot(time_r,[yp_d yp])
hold off
%%
delta_xp=xp_d(2:end)-xp(1:end-1);
delta_xx=xx_d(2:end)-xx(1:end-1);
delta_xsx=xsx_d(2:end)-xsx(1:end-1);
delta_corr_x=k1*delta_xx+k2*delta_xsx+k3*delta_xp;
xp_=[xp_;xp_d(i)+delta_corr_x];
xp1_=[xp1_;xp1_d(i)+delta_corr_x];
xp2_=[xp2_;xp2_d(i)+delta_corr_x];
Fx_=[Fx_;-(xx(end)-xp_(end))*w^2*wpg_param.m];
Fx1_=[Fx1_;0.5*Fx_(end)];
Fx2_=[Fx2_;(1-0.5)*Fx_(end)];

delta_yp=yp_d(i)-yp(end);
delta_yx=yx_d(i)-yx(end);
delta_ysx=ysx_d(i)-ysx(end);
delta_corr_y=k1*delta_yx+k2*delta_ysx+k3*delta_yp;
yp_=[yp_;yp_d(i)+delta_corr_y];
yp1_=[yp1_;yp1_d(i)+delta_corr_y];
yp2_=[yp2_;yp2_d(i)+delta_corr_y];
Fy_=[Fy_;-(yx(end)-yp_(end))*w^2*wpg_param.m];
Fy1_=[Fy1_;0.5*Fy_(end)];
Fy2_=[Fy2_;(1-0.5)*Fy_(end)];

Fz1_=[Fz1_;Fz1_d(i)];
Fz2_=[Fz2_;Fz2_d(i)];