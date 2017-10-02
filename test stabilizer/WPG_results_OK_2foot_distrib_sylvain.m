test_movement=true;
if test_movement
    duration=1.7;%in second
else
    duration=1.2;%in second
%     duration=1.4;%in second
end
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

if test_movement
    xp_d=trajectories_zmp.xpzmp(1:duration_disc);
    yp_d=trajectories_zmp.ypzmp(1:duration_disc);

    xp1_d=trajectories_zmp.xpzmp1(1:duration_disc);
    yp1_d=trajectories_zmp.ypzmp1(1:duration_disc);

    xp2_d=trajectories_zmp.xpzmp2(1:duration_disc);
    yp2_d=trajectories_zmp.ypzmp2(1:duration_disc);
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
zx=[0];

xsx=[xsx_d(1)];
ysx=[ysx_d(1)];
zsx=[0];

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
noise1=wgn(duration_disc,1,1,2)*0;
noise2=wgn(duration_disc,1,1,2)*0;
noise3=wgn(duration_disc,1,1,2)*0;
noise4=wgn(duration_disc,1,1,2)*0;
noise5=wgn(duration_disc,1,1,2)*0;
noise6=wgn(duration_disc,1,1,2)*0;

noise7=wgn(duration_disc,1,1,1e-4)*0;
noise8=wgn(duration_disc,1,1,1e-4)*0;
noise9=wgn(duration_disc,1,1,1e-4)*0;
noise10=wgn(duration_disc,1,1,1e-4)*0;
noise11=wgn(duration_disc,1,1,1e-4)*0;
noise12=wgn(duration_disc,1,1,1e-4)*0;

toto3=angleact_tot1(:,1)';
toto4=displ_tot1(:,1)';
    
toto5=angleact_tot2(:,1)';
toto6=displ_tot2(:,1)';
%%
% addpath ./test' stabilizer'/
[sole_1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_init(0.65,80000*4,0.3);
[sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,1,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);

[sole_2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_init(0.65,80000*4,0.3);
[sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,angleact2,displ2,1,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);

%%
figure(1);clf;axis equal;

xv2 = [wpg_param.step_number_pankle_fixed(2,2)+(wpg_param.fronttoankle-wpg_param.sole_margin);...
        wpg_param.step_number_pankle_fixed(2,2)+(wpg_param.fronttoankle-wpg_param.sole_margin); ...
        wpg_param.step_number_pankle_fixed(2,2)-(wpg_param.backtoankle-wpg_param.sole_margin);...
        wpg_param.step_number_pankle_fixed(2,2)-(wpg_param.backtoankle-wpg_param.sole_margin); ...
        ];
yv2 = [wpg_param.step_number_pankle_fixed(2,3)-(wpg_param.inttoankle-wpg_param.sole_margin);...
        wpg_param.step_number_pankle_fixed(2,3)+(wpg_param.exttoankle-wpg_param.sole_margin); ...
        wpg_param.step_number_pankle_fixed(2,3)+(wpg_param.exttoankle-wpg_param.sole_margin);...
        wpg_param.step_number_pankle_fixed(2,3)-(wpg_param.inttoankle-wpg_param.sole_margin); ...
        ];
xv2 = [xv2 ; xv2(1)]; yv2 = [yv2 ; yv2(1)];
hold on;plot(xv2,yv2);hold off

xv1 = [wpg_param.step_number_pankle_fixed(1,2)-(wpg_param.backtoankle-wpg_param.sole_margin); ...
        wpg_param.step_number_pankle_fixed(1,2)-(wpg_param.backtoankle-wpg_param.sole_margin);...
        wpg_param.step_number_pankle_fixed(1,2)+(wpg_param.fronttoankle-wpg_param.sole_margin); ...
        wpg_param.step_number_pankle_fixed(1,2)+(wpg_param.fronttoankle-wpg_param.sole_margin);...
        ];
yv1 = [wpg_param.step_number_pankle_fixed(1,3)+(wpg_param.inttoankle-wpg_param.sole_margin); ...
        wpg_param.step_number_pankle_fixed(1,3)-(wpg_param.exttoankle-wpg_param.sole_margin);...
        wpg_param.step_number_pankle_fixed(1,3)-(wpg_param.exttoankle-wpg_param.sole_margin); ...
        wpg_param.step_number_pankle_fixed(1,3)+(wpg_param.inttoankle-wpg_param.sole_margin);...
        ];
xv1 = [xv1 ; xv1(1)]; yv1 = [yv1 ; yv1(1)];
hold on;plot(xv1,yv1);hold off

xv=[xv1(1:end-1);xv2(1:end-1)];
yv=[yv1(1:end-1);yv2(1:end-1)];
% xv = [wpg_param.step_number_pankle_fixed(1,2)-(wpg_param.backtoankle-wpg_param.sole_margin); ...
%         wpg_param.step_number_pankle_fixed(1,2)+(wpg_param.fronttoankle-wpg_param.sole_margin);...
%         wpg_param.step_number_pankle_fixed(2,2)+(wpg_param.fronttoankle-wpg_param.sole_margin); ...
%         wpg_param.step_number_pankle_fixed(2,2)-(wpg_param.backtoankle-wpg_param.sole_margin);...
%         ];
% yv = [wpg_param.step_number_pankle_fixed(1,3)-(wpg_param.exttoankle-wpg_param.sole_margin); ...
%         wpg_param.step_number_pankle_fixed(1,3)-(wpg_param.exttoankle-wpg_param.sole_margin);...
%         wpg_param.step_number_pankle_fixed(2,3)+(wpg_param.exttoankle-wpg_param.sole_margin); ...
%         wpg_param.step_number_pankle_fixed(2,3)+(wpg_param.exttoankle-wpg_param.sole_margin);...
%         ];
xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
hold on;plot(xv,yv);hold off

hold on;plot(wpg_param.step_number_pankle_fixed(1:2,2),wpg_param.step_number_pankle_fixed(1:2,3),'+');hold off

%%
for i=2:duration_disc
    i
%     keyboard
    %% Perturbation
%     if i>=40 && i<151
% %         xx_d(i)=xx_d(i)-0.001;
%         xp_d(i)=xp_d(i)-0.04;
%     end
    if ~test_movement
        if i>=20 && i<300
            xx_d(i)=xx_d(i)+0.03;
            xp_d(i)=xp_d(i)+0.03;

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
    a=wpg_param.step_number_pankle_fixed(2,2:3)-wpg_param.step_number_pankle_fixed(1,2:3);
    a=a(2)/a(1);
    [xp_ortho yp_ortho]=intersec_line(a,wpg_param.step_number_pankle_fixed(1,2),wpg_param.step_number_pankle_fixed(1,3),-1/a,xp_(end),yp_(end));
        
    kn=[wpg_param.step_number_pankle_fixed(2,2)-wpg_param.step_number_pankle_fixed(1,2);wpg_param.step_number_pankle_fixed(2,3)-wpg_param.step_number_pankle_fixed(1,3)]'...
        *[xp_ortho-wpg_param.step_number_pankle_fixed(1,2);yp_ortho-wpg_param.step_number_pankle_fixed(1,3)];
    kd=norm([wpg_param.step_number_pankle_fixed(1,2)-wpg_param.step_number_pankle_fixed(2,2);wpg_param.step_number_pankle_fixed(1,3)-wpg_param.step_number_pankle_fixed(2,3)])^2;
    k_2=kn/kd;
    
    %% Bring back the ZMP into the Convex hull
    if k_2<=1 && k_2>=0
        %compute the convex hull vertices
%         xv = [wpg_param.step_number_pankle_fixed(1,2)-(wpg_param.backtoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(1,2)+(wpg_param.fronttoankle-wpg_param.sole_margin);...
%                 wpg_param.step_number_pankle_fixed(2,2)+(wpg_param.fronttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(2,2)-(wpg_param.backtoankle-wpg_param.sole_margin);...
%                 ];
%         yv = [wpg_param.step_number_pankle_fixed(1,3)-(wpg_param.exttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(1,3)-(wpg_param.exttoankle-wpg_param.sole_margin);...
%                 wpg_param.step_number_pankle_fixed(2,3)+(wpg_param.exttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(2,3)+(wpg_param.exttoankle-wpg_param.sole_margin);...
%                 ];
%         xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
        
        %bring the zmp into the convex hull defined  by the 2 feet if
        %needed
        %we search the closest point to ZMP projected on the convex
        %hull edges along the perpendicular to ankle segment
        [xp_0 yp_0]=projection_convex(xp_(end),yp_(end),xp_ortho,yp_ortho,xv,yv);
    elseif k_2>1
        k_2=1;
        
%         xv = [wpg_param.step_number_pankle_fixed(2,2)-(wpg_param.backtoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(2,2)+(wpg_param.fronttoankle-wpg_param.sole_margin);...
%                 wpg_param.step_number_pankle_fixed(2,2)+(wpg_param.fronttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(2,2)-(wpg_param.backtoankle-wpg_param.sole_margin);...
%                 ];
%         yv = [wpg_param.step_number_pankle_fixed(2,3)-(wpg_param.inttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(2,3)-(wpg_param.inttoankle-wpg_param.sole_margin);...
%                 wpg_param.step_number_pankle_fixed(2,3)+(wpg_param.exttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(2,3)+(wpg_param.exttoankle-wpg_param.sole_margin);...
%                 ];
%         xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
        
        %bring the zmp into the convex hull defined 1 foot if needed
        %we search the closest point to ZMP projected on the convex
        %hull edges in the direction of the ankle
        [xp_0 yp_0]=projection_convex(xp_(end),yp_(end),...
                                    wpg_param.step_number_pankle_fixed(2,2),wpg_param.step_number_pankle_fixed(2,3),...
                                    xv2,yv2);
    else
        k_2=0;
%         xv = [wpg_param.step_number_pankle_fixed(1,2)-(wpg_param.backtoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(1,2)+(wpg_param.fronttoankle-wpg_param.sole_margin);...
%                 wpg_param.step_number_pankle_fixed(1,2)+(wpg_param.fronttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(1,2)-(wpg_param.backtoankle-wpg_param.sole_margin);...
%                 ];
%         yv = [wpg_param.step_number_pankle_fixed(1,3)+(wpg_param.inttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(1,3)+(wpg_param.inttoankle-wpg_param.sole_margin);...
%                 wpg_param.step_number_pankle_fixed(1,3)-(wpg_param.exttoankle-wpg_param.sole_margin); ...
%                 wpg_param.step_number_pankle_fixed(1,3)-(wpg_param.exttoankle-wpg_param.sole_margin);...
%                 ];
%         xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
        
        %bring the zmp into the convex hull defined 1 foot if needed
        %we search the closest point to ZMP projected on the convex
        %hull edges in the direction of the ankle
        [xp_0 yp_0]=projection_convex(xp_(end),yp_(end),...
                                    wpg_param.step_number_pankle_fixed(1,2),wpg_param.step_number_pankle_fixed(1,3),...
                                    xv1,yv1);
    end
    
    %% compute the position of ZMP1&2
    if k_2<1 && k_2>0
        N=sum(wpg_param.step_number_pankle_fixed(1:2,2:3))/2;
        
        a1=[xp_0 yp_0]-wpg_param.step_number_pankle_fixed(1,2:3);
        a1=a1(2)/a1(1);
        [xp_2 yp_2]=intersec_line(-1/a,wpg_param.step_number_pankle_fixed(2,2),wpg_param.step_number_pankle_fixed(2,3),a1,N(1),N(2));
        
        a2=[xp_0 yp_0]-wpg_param.step_number_pankle_fixed(2,2:3);
        a2=a2(2)/a2(1);
        [xp_1 yp_1]=intersec_line(-1/a,wpg_param.step_number_pankle_fixed(1,2),wpg_param.step_number_pankle_fixed(1,3),a2,N(1),N(2));
        %TODO saturation of ZMP1&2
        if k_2>0.5
            [xp_1 yp_1]=projection_convex(xp_1,yp_1,...
                                            wpg_param.step_number_pankle_fixed(1,2),wpg_param.step_number_pankle_fixed(1,3),...
                                            xv1,yv1);
            xp_2=xp_1-1/k_2*(xp_1-xp_0);
            yp_2=yp_1-1/k_2*(yp_1-yp_0);
        else
            [xp_2 yp_2]=projection_convex(xp_2,yp_2,...
                                            wpg_param.step_number_pankle_fixed(2,2),wpg_param.step_number_pankle_fixed(2,3),...
                                            xv2,yv2);
            xp_1=1/(1-k_2)*(xp_0-k_2*xp_2);
            yp_1=1/(1-k_2)*(yp_0-k_2*yp_2);

        end
    elseif k_2==1
        xp_2=xp_0;yp_2=yp_0;
        xp_1=wpg_param.step_number_pankle_fixed(1,2);yp_1=wpg_param.step_number_pankle_fixed(1,3);
    else
        xp_1=xp_0;yp_1=yp_0;
        xp_2=wpg_param.step_number_pankle_fixed(2,2);yp_2=wpg_param.step_number_pankle_fixed(2,3);
    end
    %%
    fx_=-(xx(end)-xp_0)*w^2*wpg_param.m;
    fy_=-(yx(end)-yp_0)*w^2*wpg_param.m;
    
    fx_1=fx_*(1-k_2);
    fx_2=fx_*k_2;
    fy_1=fy_*(1-k_2);
    fy_2=fy_*k_2;

    fz_1=min(max(Fz_(1)*(1-k_2),10),wpg_param.mg-10);
    fz_2=min(max(Fz_(1)*k_2,10),wpg_param.mg-10);
    %%
    figure(1);
    hold on;
    plot(xp_(end),yp_(end),'o')
    plot(xp_0,yp_0,'ro')
    plot(xp_0,yp_0,'r*')
    plot(xp_ortho,yp_ortho,'go')
    plot([xp_(end);xp_0],[yp_(end);yp_0],'k')
    
    plot(xx(end),yx(end),'*')
    
    plot(xp_1,yp_1,'mo')
    plot(xp_2,yp_2,'mo')
    hold off
    
    %%
    xp1_=[xp1_;xp_1];
    xp2_=[xp2_;xp_2];
    xp_filtered=xp_0;
    Fx_=[Fx_;fx_];

    yp1_=[yp1_;yp_1];
    yp2_=[yp2_;yp_2];
    
    Fz1_=[Fz1_;fz_1];
    Fz2_=[Fz2_;fz_2];   
    
    yp_filtered=yp_0;    
    Fy_=[Fy_;fy_];
    
    Fx1_=[Fx1_;fx_1];
    Fx2_=[Fx2_;fx_2];
    Fy1_=[Fy1_;fy_1];
    Fy2_=[Fy2_;fy_2];
    
    xp_smooth=xp_;
    yp_smooth=yp_;
    
    %% mid level control
%     delta_F1=[Fx1(end);Fy1(end);Fz1(end);xp1(end);yp1(end);0]-[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0];
    delta_F1=[Fx1(end)+noise1(i);Fy1(end)+noise2(i);Fz1(end)+noise3(i);xp1(end);yp1(end);0]-[Fx1_(i);Fy1_(i);Fz1_(i);xp1_(i);yp1_(i);0];
%     delta_F1=[Fx1_d(i);Fy1_d(i);Fz1_d(i);xp1_d(i);yp1_d(i);0]-[Fx1_(end);Fy1_(end);Fz1_(end);xp1_(end);yp1_(end);0];

%     displangleact1=[displ_tot1(:,end);angleact_tot1(:,end)]-Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1*delta_F1*(1/f*1/0.05);
    displangleact1=[toto4(end,:)';toto3(end,:)']-Jmat.step_1.J_tot((40)*6+1:(40)*6+6,:)^-1*delta_F1*(1/f*1/0.005);


    displ1=displangleact1(1:3,:);
    angleact1=displangleact1(4:6,:);
    
    angleact_tot1 = [angleact_tot1 angleact1];
    displ_tot1 = [displ_tot1 displ1];
    
%     delta_F2=[Fx2(end);Fy2(end);Fz2(end);xp2(end);yp2(end);0]-[Fx2_d(i);Fy2_d(i);Fz2_d(i);xp2_d(i);yp2_d(i);0];
    delta_F2=[Fx2(end)+noise4(i);Fy2(end)+noise5(i);Fz2(end)+noise6(i);xp2(end);yp2(end);0]-[Fx2_(i);Fy2_(i);Fz2_(i);xp2_(i);yp2_(i);0];
%     delta_F2=[Fx2_d(i);Fy2_d(i);Fz2_d(i);xp2_d(i);yp2_d(i);0]-[Fx2_(end);Fy2_(end);Fz2_(end);xp2_(end);yp2_(end);0];

    
%     displangleact2=[displ_tot2(:,end);angleact_tot2(:,end)]-Jmat.step_2.J_tot((40)*6+1:(40)*6+6,:)^-1*delta_F2*(1/f*1/0.05);
    displangleact2=[toto6(end,:)';toto5(end,:)']-Jmat.step_2.J_tot((40)*6+1:(40)*6+6,:)^-1*delta_F2*(1/f*1/0.005);


    displ2=displangleact2(1:3,:);
    angleact2=displangleact2(4:6,:);
    
    angleact_tot2 = [angleact_tot2 angleact2];
    displ_tot2 = [displ_tot2 displ2];
    %% Lags on foot movement
    toto3=toto3+(angleact_tot1(:,end)'-toto3)*1/f/0.05;
    toto4=toto4+(displ_tot1(:,end)'-toto4)*1/f/0.05;
    
    toto5=toto5+(angleact_tot2(:,end)'-toto5)*1/f/0.05;
    toto6=toto6+(displ_tot2(:,end)'-toto6)*1/f/0.05;
    
%     toto3=angleact_tot1(:,end)';
%     toto4=displ_tot1(:,end)';
%     
%     toto5=angleact_tot2(:,end)';
%     toto6=displ_tot2(:,end)';
    %% robot simulation
%     [sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,angleact1,displ1,i,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);
    [sole_1 Z1 Ftot1 param_sopt_1 ABe_1 Pg_1 FcFreeSurf_1 D_1]=flexible_sole_model_update(sole_1,toto3,toto4,i,param_sopt_1,ABe_1,Pg_1,FcFreeSurf_1,D_1,false);

    Fx1=[Fx1;Ftot1(1,end)];
    Fy1=[Fy1;Ftot1(2,end)];
    Fz1=[Fz1;Ftot1(3,end)];
    xp1_sole=[xp1_sole;Z1(1,end)];
    yp1_sole=[yp1_sole;Z1(2,end)];
    
%     [sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,angleact2,displ2,i,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);
    [sole_2 Z2 Ftot2 param_sopt_2 ABe_2 Pg_2 FcFreeSurf_2 D_2]=flexible_sole_model_update(sole_2,toto5,toto6,i,param_sopt_2,ABe_2,Pg_2,FcFreeSurf_2,D_2,false);

    Fx2=[Fx2;Ftot2(1,end)];
    Fy2=[Fy2;Ftot2(2,end)];
    Fz2=[Fz2;Ftot2(3,end)];
    xp2_sole=[xp2_sole;Z2(1,end)];
    yp2_sole=[yp2_sole;Z2(2,end)];
    
    xp1=xp1_sole+traslx_1+noise7(i);
    yp1=yp1_sole+trasly_1+noise8(i);
    
    xp2=xp2_sole+traslx_2+noise9(i);
    yp2=yp2_sole+trasly_2+noise10(i);
    
    Fx=[Fx;Fx1(end)+Fx2(end)];
    Fy=[Fy;Fy1(end)+Fy2(end)];
    Fz=[Fz;Fz1(end)+Fz2(end)];
    
    xax=[xax;time(i) -Fx(end)/wpg_param.m];
    yax=[yax;time(i) -Fy(end)/wpg_param.m];
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
plot(time_r,[Fz_d Fz Fz1_+Fz2_])
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
plot(time_r,[Fx1_d Fx1 Fx1_])
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
plot(time_r,[Fy1_d Fy1 Fy1_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('Reaction force1 in z-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fz1_d Fz1 Fz1_])
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
plot(time_r,[Fx2_d Fx2 Fx2_])
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
plot(time_r,[Fy2_d Fy2 Fy2_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('Reaction force2 in z-axis')
xlabel('t(s)')
ylabel('F(N)')
hold on
plot(time_r,[Fz2_d Fz2 Fz2_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMP in x-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[xp_d xp(:,2) xp_smooth (Fz1./Fz.*xp1+Fz2./Fz.*xp2)])
hleg = legend('WPG desired','robot (simulation)','WPG modified*','robot (ED)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMP in y-axis')
xlabel('t(s)')
ylabel('y(m)')
hold on
plot(time_r,[yp_d yp(:,2) yp_smooth (Fz1./Fz.*yp1+Fz2./Fz.*yp2)])
hleg = legend('WPG desired','robot (simulation)','WPG modified*','robot (ED)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMPs in x-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[xp1 xp2 (Fz1./Fz.*xp1+Fz2./Fz.*xp2)])
plot(time_r,[xp1_ xp2_ (Fz1_./Fz_.*xp1_+Fz2_./Fz_.*xp2_)],'--')
hleg = legend('ZMP1 (ED)','ZMP2 (ED)','ZMP (ED)','ZMP1* (control)','ZMP2* (control)','ZMP* (control)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMPs in y-axis')
xlabel('t(s)')
ylabel('y(m)')
hold on
plot(time_r,[yp1 yp2 (Fz1./Fz.*yp1+Fz2./Fz.*yp2)])
plot(time_r,[yp1_ yp2_ (Fz1_./Fz_.*yp1_+Fz2_./Fz_.*yp2_)],'--')
hleg = legend('ZMP1 (ED)','ZMP2 (ED)','ZMP (ED)','ZMP1* (control)','ZMP2* (control)','ZMP* (control)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMP1 in x-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[xp1_d xp1 xp1_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMP1 in y-axis')
xlabel('t(s)')
ylabel('y(m)')
hold on
plot(time_r,[yp1_d yp1 yp1_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMP2 in x-axis')
xlabel('t(s)')
ylabel('x(m)')
hold on
plot(time_r,[xp2_d xp2 xp2_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('ZMP2 in y-axis')
xlabel('t(s)')
ylabel('y(m)')
hold on
plot(time_r,[yp2_d yp2 yp2_])
hleg = legend('WPG desired','robot (ED)','WPG modified*','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

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
ylabel('y(m)')
hold on
plot(time_r,[yx_d yx])
hleg = legend('WPG desired','robot (simulation)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

figure()
title('COM in z-axis')
xlabel('t(s)')
ylabel('z(m)')
hold on
plot(time_r,[zx_d zx]+0.782)
hleg = legend('WPG desired','robot (simulation)','Location','East');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
hold off

% figure()
% hold on
% plot(time_r,[xp(:,2) xx (xx+xsx/w)])
% hold off
% 
% figure()
% hold on
% plot(time_r,[yp(:,2) yx (yx+ysx/w)])
% hold off
%%
close all
figure()
hold on
plot(time_r(1:length(Fy)),[Fx1_d(1:length(Fy))+Fx2_d(1:length(Fy)) Fx Fx_])
hold off

figure()
hold on
plot(time_r(1:length(Fy)),[Fy1_d(1:length(Fy))+Fy2_d(1:length(Fy)) Fy Fy_])
hold off

figure()
hold on
plot(time_r(1:length(Fy)),[Fz1_d(1:length(Fy))+Fz2_d(1:length(Fy)) Fz])
hold off

figure()
hold on
plot(time_r(1:length(Fy)),[xp_d(1:length(Fy)) xp(:,2) xp_])
hold off

figure()
hold on
plot(time_r(1:length(Fy)),[yp_d(1:length(Fy)) yp(:,2) yp_])
hold off

figure()
hold on
plot(time_r(1:length(xx)),[xx_d(1:length(xx)) xx])
hold off

figure()
hold on
plot(time_r(1:length(yx)),[yx_d(1:length(yx)) yx])
hold off
%%
close all
figure()
hold on
plot(time_r(1:length(Fy)),[Fz1_d Fz1_ Fz1_d+pbc_tot(:,1)])
hold off

figure()
hold on
plot(time_r(1:length(Fy)),[Fz2_d Fz2_ Fz2_d+pbc_tot(:,2)])
hold off

figure()
hold on
plot(time_r(1:length(Fy)),[yp1_d yp1_ yp1_d+pbc_tot(:,3)])
plot(time_r(1:length(Fy)),[yp1_],'*g')
hold off

figure()
hold on
plot(time_r(1:length(Fy)),[yp2_d yp2_ yp2_d+pbc_tot(:,4)])
plot(time_r(1:length(Fy)),[yp2_],'*g')
hold off