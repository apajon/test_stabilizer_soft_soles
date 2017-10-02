test_movement=false;
if test_movement
%     duration=0.8;%in second
    duration=7.7;%in second
    duration=0.2;
else
    duration=2;%in second
%     duration=1.4;%in second
    wpg_param.dt_type_phase=zeros(duration*200+1,1);
end
f=wpg_param.frequency;
duration_disc=duration*f;
time=[0:duration_disc-1]'*1/f;

Tp=0.05;

w=wpg_param.w;

poles=[-13,-3,-w];
% poles=[-w,-w,-w];
% poles=[-13,-13,-w];
% poles=[-20,-3,-w];
% poles=[-20,-20,-w];
% poles=[-100,-3,-w];
% poles=[-100,-100,-w];

% k1=-poles(1)*poles(2)-1;
% k2=(poles(1)+poles(2))/-poles(3);
% k3=-w;
k3=-(sum(poles))*Tp-1;
k1=prod(poles)*Tp/w^2-(1+k3);
k2=-Tp*(1+(poles(1)*poles(2)+poles(1)*poles(3)+poles(3)*poles(2))/w^2);

Tgain=0.001;


perturb=zeros(duration_disc,1);
perturb=[time perturb];
%%
time_r=time;

if test_movement
    xp_d=trajectories_zmp.xpzmp(1:duration_disc);
    yp_d=trajectories_zmp.ypzmp(1:duration_disc);

    xp1_d=trajectories_zmp.xpzmp1;
    yp1_d=trajectories_zmp.ypzmp1;

    xp2_d=trajectories_zmp.xpzmp2;
    yp2_d=trajectories_zmp.ypzmp2;
else
    % k=diag(zmp.k_diag);
    % k(sum(wpg_param.discretization(1:wpg_param.nbpolypi))+1:end)=[];
    xp_d=repmat(trajectories_zmp.xpzmp(1),duration_disc,1);
    yp_d=repmat(0,duration_disc,1);

    xp1_d=repmat(trajectories_zmp.xpzmp1(1),duration_disc,1);
    yp1_d=repmat(trajectories_zmp.ypzmp1(1),duration_disc,1);

    xp2_d=repmat(trajectories_zmp.xpzmp2(1),duration_disc,1);
    yp2_d=repmat(trajectories_zmp.ypzmp2(1),duration_disc,1);
end

footy=wpg_param.inttoankle+wpg_param.exttoankle;
translx_R=wpg_param.step_number_pankle_fixed(1,2)-wpg_param.backtoankle;
transly_R=wpg_param.step_number_pankle_fixed(1,3)-(footy/2-wpg_param.inttoankle);
xpR_d_sole=xp1_d-translx_R;
ypR_d_sole=yp1_d-transly_R;

translx_L=wpg_param.step_number_pankle_fixed(1,2)-wpg_param.backtoankle;
transly_L=wpg_param.step_number_pankle_fixed(2,3)+(footy/2-wpg_param.inttoankle);
xpL_d_sole=xp2_d-translx_L;
ypL_d_sole=yp2_d-transly_L;
%%
% if test_movement
%     Fx_d=repmat(0,duration_disc,1);
%     Fy_d=repmat(0,duration_disc,1);
%     Fz_d=repmat(wpg_param.mg,duration_disc,1);
% 
%     Fx1_d=0.5.*Fx_d;
%     Fy1_d=0.5.*Fy_d;
%     Fz1_d=0.5.*Fz_d;
% 
%     Fx2_d=0.5.*Fx_d;
%     Fy2_d=0.5.*Fy_d;
%     Fz2_d=0.5.*Fz_d;
% else
    Fx_d=repmat(0,duration_disc,1);
    Fy_d=repmat(0,duration_disc,1);
    Fz_d=repmat(wpg_param.mg,duration_disc,1);

    Fx1_d=0.5.*Fx_d;
    Fy1_d=0.5.*Fy_d;
    Fz1_d=0.5.*Fz_d;

    Fx2_d=0.5.*Fx_d;
    Fy2_d=0.5.*Fy_d;
    Fz2_d=0.5.*Fz_d;
% end
%%
if test_movement
    xx_d=trajectories_zmp.xpcom(1:duration_disc);
    yx_d=trajectories_zmp.ypcom(1:duration_disc);
    zx_d=repmat(0,duration_disc,1);

    xsx_d=zmp.A_xcom_spd(1:duration_disc,:)*wpg_param.psa_abcdDSP(1:size(zmp.A_xcom_spd,2))+zmp.B_xcom_spd(1:duration_disc,:);
    ysx_d=zmp.A_ycom_spd(1:duration_disc,:)*wpg_param.psa_abcdDSP(end/2+1:end/2+size(zmp.A_ycom_spd,2))+zmp.B_ycom_spd(1:duration_disc,:);
else
    xx_d=repmat(trajectories_zmp.xpcom(1),duration_disc,1);
    yx_d=repmat(0,duration_disc,1);
    zx_d=repmat(0,duration_disc,1);

    xsx_d=repmat(0,duration_disc,1);
    ysx_d=repmat(0,duration_disc,1);
end
%%
xp=[xp_d(1)];
xpR=[xp1_d(1)];
xpL=[xp2_d(1)];

yp=[yp_d(1)];
ypR=[yp1_d(1)];
ypL=[yp2_d(1)];

xp_=[xp_d(1)];
xpR_=[xp1_d(1)];
xpL_=[xp2_d(1)];

yp_=[yp_d(1)];
ypR_=[yp1_d(1)];
ypL_=[yp2_d(1)];

xpR_sole=[xpR_d_sole(1)];
xpL_sole=[xpL_d_sole(1)];

ypR_sole=[ypR_d_sole(1)];
ypL_sole=[ypL_d_sole(1)];

Fx=[Fx_d(1)];
Fy=[Fy_d(1)];
Fz=[Fz_d(1)];

Fx_=[Fx_d(1)];
Fy_=[Fy_d(1)];
Fz_=[Fz_d(1)];

FxR=[Fx1_d(1)];
FyR=[Fy1_d(1)];
FzR=[Fz1_d(1)];

FxL=[Fx2_d(1)];
FyL=[Fy2_d(1)];
FzL=[Fz2_d(1)];

xx=[xx_d(1)];
yx=[yx_d(1)];
zx=[0];

xsx=[xsx_d(1)];
ysx=[ysx_d(1)];
zsx=[0];

xax=[0];
yax=[0];

xp_=[xp_d(1)];
xpR_=[xp1_d(1)];
xpL_=[xp2_d(1)];

yp_=[yp_d(1)];
ypR_=[yp1_d(1)];
ypL_=[yp2_d(1)];

xx_=[xx_d(1)];
yx_=[yx_d(1)];

xsx_=[xsx_d(1)];
ysx_=[ysx_d(1)];

Fx_=[Fx_d(1)];
Fy_=[Fy_d(1)];
Fz_=[Fz_d];

FxR_=[Fx1_d(1)];
FyR_=[Fy1_d(1)];
FzR_=[Fz1_d(1)];

FxL_=[Fx2_d(1)];
FyL_=[Fy2_d(1)];
FzL_=[Fz2_d(1)];

xp_filtered=[xp_d(1)];
yp_filtered=[yp_d(1)];    
    
xp_satured=[xp_d(1)];
yp_satured=[yp_d(1)];

xp_control=[xp_d(1)];
yp_control=[yp_d(1)];

xp_ED=[xp_d(1)];
yp_ED=[yp_d(1)];
%%
load('Simulation_3D/results/J_0.mat');
angleactR = Jmat.step_1.angleact_tot(:,1);
displR = Jmat.step_1.displ_tot(:,1);
angleact_totR = [angleactR];
displ_totR = [displR];

angleactL = Jmat.step_2.angleact_tot(:,1);
displL = Jmat.step_2.displ_tot(:,1);
angleact_totL = [angleactL];
displ_totL = [displL];
%%

% xp_d_init=xp_d(1);
% xp_d=[time xp_d-xp_d_init];
% yp_d_init=yp_d(1);
% yp_d=[time yp_d-yp_d_init];
% 
% xp1_d=[time xp1_d];
% yp1_d=[time yp1_d];
% 
% xp2_d=[time xp2_d];
% yp2_d=[time yp2_d];
% 
% xx_d_init=xx_d(1);
% xx_d=[time xx_d-xx_d_init];
% yx_d_init=yx_d(1);
% yx_d=[time yx_d-yx_d_init];
% 
% xsx_d_init=xsx_d(1);
% xsx_d=[time xsx_d-xsx_d_init];
% ysx_d_init=ysx_d(1);
% ysx_d=[time ysx_d-ysx_d_init];
% 
% Fx1_d=[time Fx1_d];
% Fy1_d=[time Fy1_d];
% Fz1_d=[time Fz1_d];
% 
% Fx2_d=[time Fx2_d];
% Fy2_d=[time Fy2_d];
% Fz2_d=[time Fz2_d];
% 
% xp_=[time(1) xp_d(1)];
% xp1_=[time(1) xp1_d(1)];
% xp2_=[time(1) xp2_d(1)];
% 
% yp_=[time(1) yp_d(1)];
% yp1_=[time(1) yp1_d(1)];
% yp2_=[time(1) yp2_d(1)];
% 
% Fx_=[time(1) Fx_d(1)];
% Fy_=[time(1) Fy_d(1)];
% Fz_=[time(1) Fz_d];
% 
% Fx1_=[time(1) Fx1_d(1)];
% Fy1_=[time(1) Fy1_d(1)];
% Fz1_=[time(1) Fz1_d(1)];
% 
% Fx2_=[time(1) Fx2_d(1)];
% Fy2_=[time(1) Fy2_d(1)];
% Fz2_=[time(1) Fz2_d(1)];
% 
% xp1=[time(1) xp1];
% yp1=[time(1) yp1];
% Fx1=[time(1) Fx1];
% Fy1=[time(1) Fy1];
% Fz1=[time(1) Fz1];

xp=[time(1) xp_d(1)];
yp=[time(1) yp_d(1)];

xax=[time(1) 0];
yax=[time(1) 0];
zax=[time(1) 0];

pbc_tot=[0 0 0 0 0 0];
%%
ActiveNoise=0;
%noise on FR and FL
noise1=wgn(duration_disc,1,1,2)*ActiveNoise;
noise2=wgn(duration_disc,1,1,2)*ActiveNoise;
noise3=wgn(duration_disc,1,1,2)*ActiveNoise;
noise4=wgn(duration_disc,1,1,2)*ActiveNoise;
noise5=wgn(duration_disc,1,1,2)*ActiveNoise;
noise6=wgn(duration_disc,1,1,2)*ActiveNoise;

%noise on ZMP1&2
noise7=wgn(duration_disc,1,1,1e-4)*ActiveNoise;
noise8=wgn(duration_disc,1,1,1e-4)*ActiveNoise;
noise9=wgn(duration_disc,1,1,1e-4)*ActiveNoise;
noise10=wgn(duration_disc,1,1,1e-4)*ActiveNoise;
%noise on ZMP
noise11=wgn(duration_disc,1,1,1e-4)*ActiveNoise;
noise12=wgn(duration_disc,1,1,1e-4)*ActiveNoise;

toto3=angleact_totR(:,1)';
toto4=displ_totR(:,1)';
    
toto5=angleact_totL(:,1)';
toto6=displ_totL(:,1)';
%%
% addpath ./test' stabilizer'/
[sole_R param_sopt_R ABe_R Pg_R FcFreeSurf_R D_R]=flexible_sole_model_init(0.65,80000,0.3);
[sole_R ZR FtotR param_sopt_R ABe_R Pg_R FcFreeSurf_R D_R]=flexible_sole_model_update(sole_R,angleactR,displR,1,param_sopt_R,ABe_R,Pg_R,FcFreeSurf_R,D_R,false);

[sole_L param_sopt_L ABe_L Pg_L FcFreeSurf_L D_L]=flexible_sole_model_init(0.65,80000,0.3);
[sole_L ZL FtotL param_sopt_L ABe_L Pg_L FcFreeSurf_L D_L]=flexible_sole_model_update(sole_L,angleactL,displL,1,param_sopt_L,ABe_L,Pg_L,FcFreeSurf_L,D_L,false);

%%
RorL=wpg_param.firstSS;%-1 left, +1 right
RFootIndex=1;
LFootIndex=2;

[xv yv xvR_reduced yvR_reduced xvL_reduced yvL_reduced xvR yvR xvL yvL]=compute_convexhull_vertices(wpg_param.step_number_pankle_fixed(RFootIndex,:),wpg_param.step_number_pankle_fixed(LFootIndex,:),wpg_param);
%%
% i=1;
% if true
%     i=i+1;
% for i=2:353
for i=2:duration_disc
% for i=341:duration_disc
    i
%     size(yp1_)
%     size(yp2_)
%     size(Fz1_)
%     size(Fz2_)
%     keyboard
    %% Perturbation
%     if i>=40 && i<151
% %         xx_d(i)=xx_d(i)-0.001;
%         xp_d(i)=xp_d(i)-0.04;
%     end
    if ~test_movement
        if i>=20 && i<1200
            xx_d(i)=xx_d(i)+0.00;
            xp_d(i)=xp_d(i)+0.00;

            yx_d(i)=yx_d(i)-0.05;
            yp_d(i)=yp_d(i)-0.05;
        end
    end
    %% high level control
    delta_xp=xp_d(i)-xp(end,2);
    delta_xx=xx_d(i)-xx(end);
    delta_xsx=xsx_d(i)-xsx(end);
    delta_corr_x=k1*delta_xx+k2*delta_xsx+k3*delta_xp;
    xp_=[xp_;xp_d(i)+delta_corr_x];
    
    delta_yp=yp_d(i)-yp(end,2);
    delta_yx=yx_d(i)-yx(end);
    delta_ysx=ysx_d(i)-ysx(end);
    delta_corr_y=k1*delta_yx+k2*delta_ysx+k3*delta_yp;
    yp_=[yp_;yp_d(i)+delta_corr_y];
    
    %% compute the orthogonal projection of ZMP on ankle segment
    if wpg_param.dt_type_phase(i)==0
        if wpg_param.dt_type_phase(i-1)~=wpg_param.dt_type_phase(i)
            if RorL==-1
                RFootIndex=RFootIndex+2;
                
                translx_R=wpg_param.step_number_pankle_fixed(RFootIndex,2)-wpg_param.backtoankle;
                transly_R=wpg_param.step_number_pankle_fixed(RFootIndex,3)-(footy/2-wpg_param.inttoankle);
            else
                LFootIndex=LFootIndex+2;
                
                translx_L=wpg_param.step_number_pankle_fixed(LFootIndex,2)-wpg_param.backtoankle;
                transly_L=wpg_param.step_number_pankle_fixed(LFootIndex,3)+(footy/2-wpg_param.inttoankle);
            end
            RorL=-1*RorL;
            [xv yv xvR_reduced yvR_reduced xvL_reduced yvL_reduced xvR yvR xvL yvL]=compute_convexhull_vertices(wpg_param.step_number_pankle_fixed(RFootIndex,:),wpg_param.step_number_pankle_fixed(LFootIndex,:),wpg_param);
        end
        
        a=wpg_param.step_number_pankle_fixed(LFootIndex,2:3)-wpg_param.step_number_pankle_fixed(RFootIndex,2:3);
        a=a(2)/a(1);
        [xp_ortho yp_ortho]=intersec_line(a,wpg_param.step_number_pankle_fixed(RFootIndex,2),wpg_param.step_number_pankle_fixed(RFootIndex,3),-1/a,xp_(end),yp_(end));

        kn=[wpg_param.step_number_pankle_fixed(LFootIndex,2)-wpg_param.step_number_pankle_fixed(RFootIndex,2);wpg_param.step_number_pankle_fixed(LFootIndex,3)-wpg_param.step_number_pankle_fixed(RFootIndex,3)]'...
            *[xp_ortho-wpg_param.step_number_pankle_fixed(RFootIndex,2);yp_ortho-wpg_param.step_number_pankle_fixed(RFootIndex,3)];
        kd=norm([wpg_param.step_number_pankle_fixed(RFootIndex,2)-wpg_param.step_number_pankle_fixed(LFootIndex,2);wpg_param.step_number_pankle_fixed(RFootIndex,3)-wpg_param.step_number_pankle_fixed(LFootIndex,3)])^2;
        k_2=kn/kd;
    else
        if RorL==+1
            k_2=0;
        else
            k_2=1;
        end
    end
    %% Bring back the ZMP into the Convex hull
    if wpg_param.dt_type_phase(i)==0
        if k_2<=1 && k_2>=0
            %bring the zmp into the convex hull defined  by the 2 feet if
            %needed
            %we search the closest point to ZMP projected on the convex
            %hull edges along the perpendicular to ankle segment
            [xp_0 yp_0]=projection_convex(xp_(end),yp_(end),xp_ortho,yp_ortho,xv,yv);
        elseif k_2>1
            k_2=1;
            %bring the zmp into the convex hull defined by 1 foot if needed
            %we search the closest point to ZMP projected on the convex
            %hull edges in the direction of the ankle
            [xp_0 yp_0]=projection_convex(xp_(end),yp_(end),...
                                        wpg_param.step_number_pankle_fixed(LFootIndex,2),wpg_param.step_number_pankle_fixed(LFootIndex,3),...
                                        xvL,yvL);
        else
            k_2=0;
            %bring the zmp into the convex hull deined by 1 foot if needed
            %we search the closest point to ZMP projected on the convex
            %hull edges in the direction of the ankle
            [xp_0 yp_0]=projection_convex(xp_(end),yp_(end),...
                                        wpg_param.step_number_pankle_fixed(RFootIndex,2),wpg_param.step_number_pankle_fixed(RFootIndex,3),...
                                        xvR,yvR);
        end
    else
        if k_2==1
            %bring the zmp into the convex hull defined by 1 foot if needed
            %we search the closest point to ZMP projected on the convex
            %hull edges in the direction of the ankle
            [xp_0 yp_0]=projection_convex(xp_(end),yp_(end),...
                                        wpg_param.step_number_pankle_fixed(LFootIndex,2),wpg_param.step_number_pankle_fixed(LFootIndex,3),...
                                        xvL,yvL);
        else
            %bring the zmp into the convex hull defined by 1 foot if needed
            %we search the closest point to ZMP projected on the convex
            %hull edges in the direction of the ankle
            [xp_0 yp_0]=projection_convex(xp_(end),yp_(end),...
                                        wpg_param.step_number_pankle_fixed(RFootIndex,2),wpg_param.step_number_pankle_fixed(RFootIndex,3),...
                                        xvR,yvR);
        end
    end
    
    %% compute the position of ZMP1&2
    if k_2<1 && k_2>0
        N=sum(wpg_param.step_number_pankle_fixed([RFootIndex LFootIndex],2:3))/2;

        a1=[xp_0 yp_0]-wpg_param.step_number_pankle_fixed(RFootIndex,2:3);
        a1=a1(2)/a1(1);
        [xp_L yp_L]=intersec_line(-1/a,wpg_param.step_number_pankle_fixed(LFootIndex,2),wpg_param.step_number_pankle_fixed(LFootIndex,3),a1,N(1),N(2));

        a2=[xp_0 yp_0]-wpg_param.step_number_pankle_fixed(LFootIndex,2:3);
        a2=a2(2)/a2(1);
        [xp_R yp_R]=intersec_line(-1/a,wpg_param.step_number_pankle_fixed(RFootIndex,2),wpg_param.step_number_pankle_fixed(RFootIndex,3),a2,N(1),N(2));
        %TODO saturation of ZMP1&2
        if k_2>0.5
            [xp_R yp_R]=projection_convex(xp_R,yp_R,...
                                            wpg_param.step_number_pankle_fixed(RFootIndex,2),wpg_param.step_number_pankle_fixed(RFootIndex,3),...
                                            xvR,yvR);
            xp_L=xp_R-1/k_2*(xp_R-xp_0);
            yp_L=yp_R-1/k_2*(yp_R-yp_0);
        else
            [xp_L yp_L]=projection_convex(xp_L,yp_L,...
                                            wpg_param.step_number_pankle_fixed(LFootIndex,2),wpg_param.step_number_pankle_fixed(LFootIndex,3),...
                                            xvL,yvL);
            xp_R=1/(1-k_2)*(xp_0-k_2*xp_L);
            yp_R=1/(1-k_2)*(yp_0-k_2*yp_L);

        end
    elseif k_2==1
        xp_L=xp_0;yp_L=yp_0;
        xp_R=wpg_param.step_number_pankle_fixed(RFootIndex,2);yp_R=wpg_param.step_number_pankle_fixed(RFootIndex,3);
    else%if k_2==0
        xp_R=xp_0;yp_R=yp_0;
        xp_L=wpg_param.step_number_pankle_fixed(LFootIndex,2);yp_L=wpg_param.step_number_pankle_fixed(LFootIndex,3);
    end
    
    %% Force repartition
    xp_0_satured=xp_0;
    yp_0_satured=yp_0;
    if wpg_param.dt_type_phase(i)==0
        fz_R=min(max(Fz_(1)*(1-k_2),30),wpg_param.mg-30);
        fz_L=min(max(Fz_(1)*k_2,30),wpg_param.mg-30);
        
        xp_0_satured=fz_R/Fz_(1)*xp_R+fz_L/Fz_(1)*xp_L;
        yp_0_satured=fz_R/Fz_(1)*yp_R+fz_L/Fz_(1)*yp_L;
    elseif RorL==+1
        fz_R=min(Fz_(1),wpg_param.mg);
        fz_L=30;
    elseif RorL==-1
        fz_R=30;
        fz_L=min(Fz_(1),wpg_param.mg);
    end
        
    fx_=(xx(end)-xp_0_satured)*w^2*wpg_param.m;
    fy_=(yx(end)-yp_0_satured)*w^2*wpg_param.m;
    
    fx_R=fx_*(1-k_2);
    fx_L=fx_*k_2;
    fy_R=fy_*(1-k_2);
    fy_L=fy_*k_2;
    
    %%
    xpR_=[xpR_;xp_R];
    xpL_=[xpL_;xp_L];

    ypR_=[ypR_;yp_R];
    ypL_=[ypL_;yp_L];
    
    FxR_=[FxR_;fx_R];
    FxL_=[FxL_;fx_L];
    FyR_=[FyR_;fy_R];
    FyL_=[FyL_;fy_L];
    FzR_=[FzR_;fz_R];
    FzL_=[FzL_;fz_L];   
    
    
    Fx_=[Fx_;fx_];
    Fy_=[Fy_;fy_];
%     Fz_=[];
    
    xp_filtered=[xp_filtered;xp_0];
    yp_filtered=[yp_filtered;yp_0];    
    
    xp_satured=[xp_satured;xp_0_satured];
    yp_satured=[yp_satured;yp_0_satured];
    
    xp_control=[xp_control;xp_(end)];
    yp_control=[yp_control;yp_(end)];
    
    %% mid level control
    if wpg_param.dt_type_phase(i)==0
        if wpg_param.dt_type_phase(i-1)~=wpg_param.dt_type_phase(i)
            if RorL==-1
                delta_FR=[FxR(end)+noise1(i);FyR(end)+noise2(i);FzR(end)+noise3(i);xpR(end);ypR(end);0]-[FxR_(i);FyR_(i);FzR_(i);xpR_(i);ypR_(i);0];
                delta_FL=-[FxL_(i);FyL_(i);FzL_(i);0;0;0];
            else
                delta_FR=-[FxR_(i);FyR_(i);FzR_(i);0;0;0];
                delta_FL=[FxL(end)+noise4(i);FyL(end)+noise5(i);FzL(end)+noise6(i);xpL(end);ypL(end);0]-[FxL_(i);FyL_(i);FzL_(i);xpL_(i);ypL_(i);0];
            end
        else
            delta_FR=[FxR(end)+noise1(i);FyR(end)+noise2(i);FzR(end)+noise3(i);xpR(end);ypR(end);0]-[FxR_(i);FyR_(i);FzR_(i);xpR_(i);ypR_(i);0];
            delta_FL=[FxL(end)+noise4(i);FyL(end)+noise5(i);FzL(end)+noise6(i);xpL(end);ypL(end);0]-[FxL_(i);FyL_(i);FzL_(i);xpL_(i);ypL_(i);0];
        end
    else
        if RorL==-1
            delta_FR=-[0;0;30;0;0;0];
            delta_FL=[FxL(end)+noise4(i);FyL(end)+noise5(i);FzL(end)+noise6(i);xpL(end);ypL(end);0]-[FxL_(i);FyL_(i);FzL_(i);xpL_(i);ypL_(i);0];
        else
            delta_FR=[FxR(end)+noise1(i);FyR(end)+noise2(i);FzR(end)+noise3(i);xpR(end);ypR(end);0]-[FxR_(i);FyR_(i);FzR_(i);xpR_(i);ypR_(i);0];
            delta_FL=-[0;0;30;0;0;0];
        end
    end
    displangleactR=[toto4(end,:)';toto3(end,:)']-Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1*delta_FR*(1/f*1/Tgain);
    displangleactL=[toto6(end,:)';toto5(end,:)']-Jmat.step_2.J_tot((40)*6+1:(40)*6+6,:)^-1*delta_FL*(1/f*1/Tgain);


    displR=displangleactR(1:3,:);
    angleactR=displangleactR(4:6,:);

    angleact_totR = [angleact_totR angleactR];
    displ_totR = [displ_totR displR];

    displL=displangleactL(1:3,:);
    angleactL=displangleactL(4:6,:);

    angleact_totL = [angleact_totL angleactL];
    displ_totL = [displ_totL displL];
    %% Lags on foot movement
%     toto3=toto3+(angleact_totR(:,end)'-toto3)*1/f/Tp;
%     toto4=toto4+(displ_totR(:,end)'-toto4)*1/f/Tp;
%     
%     toto5=toto5+(angleact_totL(:,end)'-toto5)*1/f/Tp;
%     toto6=toto6+(displ_totL(:,end)'-toto6)*1/f/Tp;

%     if wpg_param.dt_type_phase(i)==0
%         toto3=angleact_totR(:,end)';
%         toto4=displ_totR(:,end)';
% 
%         toto5=angleact_totL(:,end)';
%         toto6=displ_totL(:,end)';
%     else
%         if RorL==-1
%             toto3=(angleact_totR(:,end)');
%             toto4=(displ_totR(:,end)');
% 
%             toto5=(angleact_totL(:,end)');
%             toto6=(displ_totL(:,end)');
%         else
%             toto3=(angleact_totR(:,end)');
%             toto4=(displ_totR(:,end)');
% 
%             toto5=(angleact_totL(:,end)');
%             toto6=(displ_totL(:,end)');
%         end
%     end

    
    toto3=angleact_totR(:,end)';
    toto4=displ_totR(:,end)';
    
    toto5=angleact_totL(:,end)';
    toto6=displ_totL(:,end)';
    %% robot simulation
%     [sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,i,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);
    [sole_R ZR FtotR param_sopt_R ABe_R Pg_R FcFreeSurf_R D_R]=flexible_sole_model_update(sole_R,toto3,toto4,i,param_sopt_R,ABe_R,Pg_R,FcFreeSurf_R,D_R,false);

    FxR=[FxR;FtotR(1,end)];
    FyR=[FyR;FtotR(2,end)];
    FzR=[FzR;FtotR(3,end)];
    xpR_sole=[xpR_sole;ZR(1,end)];
    ypR_sole=[ypR_sole;ZR(2,end)];
    
%     [sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,angleact2,displ2,i,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);
    [sole_L ZL FtotL param_sopt_L ABe_L Pg_L FcFreeSurf_L D_L]=flexible_sole_model_update(sole_L,toto5,toto6,i,param_sopt_L,ABe_L,Pg_L,FcFreeSurf_L,D_L,false);

    FxL=[FxL;FtotL(1,end)];
    FyL=[FyL;FtotL(2,end)];
    FzL=[FzL;FtotL(3,end)];
    xpL_sole=[xpL_sole;ZL(1,end)];
    ypL_sole=[ypL_sole;ZL(2,end)];
    
    xpR=[xpR;xpR_sole(i)+translx_R+noise7(i)];
    ypR=[ypR;ypR_sole(i)+transly_R+noise8(i)];
    
    xpL=[xpL;xpL_sole(i)+translx_L+noise9(i)];
    ypL=[ypL;ypL_sole(i)+transly_L+noise10(i)];
    
    if wpg_param.dt_type_phase(i)==0
        Fx=[Fx;FxR(end)+FxL(end)];
        Fy=[Fy;FyR(end)+FyL(end)];
        Fz=[Fz;FzR(end)+FzL(end)];
        
        xp_ED=[xp_ED;FzR(end)/Fz(end)*xpR(end)+FzL(end)/Fz(end)*xpL(end)];
        yp_ED=[yp_ED;FzR(end)/Fz(end)*ypR(end)+FzL(end)/Fz(end)*ypL(end)];
    elseif RorL==+1
        Fx=[Fx;FxR(end)];
        Fy=[Fy;FyR(end)];
        Fz=[Fz;FzR(end)];
        
        xp_ED=[xp_ED;xpR(end)];
        yp_ED=[yp_ED;ypR(end)];
    elseif RorL==-1
        Fx=[Fx;FxL(end)];
        Fy=[Fy;FyL(end)];
        Fz=[Fz;FzL(end)];
        
        xp_ED=[xp_ED;xpL(end)];
        yp_ED=[yp_ED;ypL(end)];
    end

    
%     xax=[xax;time(i) -Fx(end)/wpg_param.m];
%     yax=[yax;time(i) -Fy(end)/wpg_param.m];
    xax=[xax;time(i) Fx(end)/wpg_param.m];
    yax=[yax;time(i) Fy(end)/wpg_param.m];
    zax=[zax;time(i) Fz(end)/wpg_param.m-wpg_param.g];
       
    %%
    simOut=sim('test stabilizer/simu_control_robot_5','StartTime','0','StopTime',sprintf('%d',time(i)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f));
%     simOut=sim('test stabilizer/simu_control_robot_6','StartTime','0','StopTime',sprintf('%d',time(i)),'SolverType','Fixed-step','FixedStep',sprintf('%d',1/f));
    simu.xx=simOut.get('xx');
    simu.yx=simOut.get('yx');
    simu.zx=simOut.get('zx');
    simu.xsx=simOut.get('xsx');
    simu.ysx=simOut.get('ysx');
    simu.zsx=simOut.get('zsx');
%     simu.xax=simOut.get('xax');
%     simu.yax=simOut.get('yax');
    xx=[xx;simu.xx.Data(end)];
    yx=[yx;simu.yx.Data(end)];
    zx=[zx;simu.zx.Data(end)];
    xsx=[xsx;simu.xsx.Data(end)];
    ysx=[ysx;simu.ysx.Data(end)];
    zsx=[zsx;simu.zsx.Data(end)];
%     xax=[xax;time(i) simu.xax.Data(end)];
%     yax=[yax;time(i) simu.yax.Data(end)];

    xp=[xp;time(i) xx(end-1)-xax(end,2)/w^2+noise11(i)];
    yp=[yp;time(i) yx(end-1)-yax(end,2)/w^2+noise12(i)];
end
%%
close all
figure()
clf
title('Reaction force in x-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fx_d Fx Fx_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
clf
title('Reaction force in y-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fy_d Fy Fy_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('Reaction force in z-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fz_d Fz Fz_d])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
clf
title('Reaction force1 in x-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fx1_d FxR FxR_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
clf
title('Reaction force1 in y-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fy1_d FyR FyR_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('Reaction force1 in z-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fz1_d FzR FzR_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
clf
title('Reaction force2 in x-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fx2_d FxL FxL_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
clf
title('Reaction force2 in y-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fy2_d FyL FyL_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('Reaction force2 in z-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fz2_d FzL FzL_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMP in x-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[xp_d xp(:,2) xp_smooth (FzR./Fz.*xpR+FzL./Fz.*xpL)])
hleg = legend('WPG desired','robot (simulation)','WPG modified*','robot (ED)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMP in y-axis')
xlabel('t(s)')
ylabel('y(m)')
hold on
plot(time_r,[yp_d yp(:,2) yp_smooth (FzR./Fz.*ypR+FzL./Fz.*ypL)])
hleg = legend('WPG desired','robot (simulation)','WPG modified*','robot (ED)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMPs in x-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[xpR xpL (FzR./Fz.*xpR+FzL./Fz.*xpL)])
plot(time_r,[xpR_ xpL_ (FzR_./Fz_.*xpR_+FzL_./Fz_.*xpL_)],'--')
plot(time_r,xp_d,'k')
hleg = legend('ZMPR (ED)','ZMPL (ED)','ZMP (ED)','ZMPR* (control)','ZMPL* (control)','ZMP* (control)','ZMP WPG','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMPs in y-axis')
xlabel('t(s)')
ylabel('y(m)')
hold on
plot(time_r,[ypR ypL (FzR./Fz.*ypR+FzL./Fz.*ypL)])
plot(time_r,[ypR_ ypL_ ((FzR_./Fz_).*ypR_+FzL_./Fz_.*ypL_)],'--')
plot(time_r,yp_d,'k')
hleg = legend('ZMPR (ED)','ZMPL (ED)','ZMP (ED)','ZMPR* (control)','ZMPL* (control)','ZMP* (control)','ZMP WPG','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

% figure()
% title('ZMP1 in x-axis')
% xlabel('t(s)')
% ylabel('x(m)')
% hold on
% plot(time_r,[xp1_d xpR xpR_])
% hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% % Make the text of the legend italic and color it brown
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
% hold off
% 
% figure()
% title('ZMP1 in y-axis')
% xlabel('t(s)')
% ylabel('y(m)')
% hold on
% plot(time_r,[yp1_d ypR ypR_])
% hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% % Make the text of the legend italic and color it brown
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
% hold off
% 
% figure()
% title('ZMP2 in x-axis')
% xlabel('t(s)')
% ylabel('x(m)')
% hold on
% plot(time_r,[xp2_d xpL xpL_])
% hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% % Make the text of the legend italic and color it brown
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
% hold off
% 
% figure()
% title('ZMP2 in y-axis')
% xlabel('t(s)')
% ylabel('y(m)')
% hold on
% plot(time_r,[yp2_d ypL ypL_])
% hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% % Make the text of the legend italic and color it brown
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
% hold off

figure()
title('COM in x-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[xx_d xx])
hleg = legend('WPG desired','robot (simulation)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('COM in y-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[yx_d yx])
hleg = legend('WPG desired','robot (simulation)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('COM in z-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[zx_d zx])
hleg = legend('WPG desired','robot (simulation)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('COM error')
xlabel('t(s)')
ylabel('y(m)')
hold on
plot(time_r,[xx-xx_d yx-yx_d])
hleg = legend('COM error x-axis','COM error y-axis','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off