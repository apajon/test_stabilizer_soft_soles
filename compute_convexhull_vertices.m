function [xv yv xvR_reduced yvR_reduced xvL_reduced yvL_reduced xvR yvR xvL yvL]=compute_convexhull_vertices(footR,footL,wpg_param)

%%foot R and L borders, used in SSP
% figure(1);clf;axis equal;
xvL = [footL(2)+(wpg_param.fronttoankle-wpg_param.sole_margin);...
        footL(2)+(wpg_param.fronttoankle-wpg_param.sole_margin); ...
        footL(2)-(wpg_param.backtoankle-wpg_param.sole_margin);...
        footL(2)-(wpg_param.backtoankle-wpg_param.sole_margin); ...
        ];
yvL = [footL(3)-(wpg_param.inttoankle-wpg_param.sole_margin);...
        footL(3)+(wpg_param.exttoankle-wpg_param.sole_margin); ...
        footL(3)+(wpg_param.exttoankle-wpg_param.sole_margin);...
        footL(3)-(wpg_param.inttoankle-wpg_param.sole_margin); ...
        ];
xvL = [xvL ; xvL(1)]; yvL = [yvL ; yvL(1)];
% hold on;plot(xvL,yvL);hold off

xvR = [footR(2)-(wpg_param.backtoankle-wpg_param.sole_margin); ...
        footR(2)-(wpg_param.backtoankle-wpg_param.sole_margin);...
        footR(2)+(wpg_param.fronttoankle-wpg_param.sole_margin); ...
        footR(2)+(wpg_param.fronttoankle-wpg_param.sole_margin);...
        ];
yvR = [footR(3)+(wpg_param.inttoankle-wpg_param.sole_margin); ...
        footR(3)-(wpg_param.exttoankle-wpg_param.sole_margin);...
        footR(3)-(wpg_param.exttoankle-wpg_param.sole_margin); ...
        footR(3)+(wpg_param.inttoankle-wpg_param.sole_margin);...
        ];
xvR = [xvR ; xvR(1)]; yvR = [yvR ; yvR(1)];
% hold on;plot(xvR,yvR);hold off

%%foot R and L intersection with perpendicular to ankle segment.
a=footL(2:3)-footR(2:3);
a_ortho=[-a(2);a(1)];
[xvL1 yvL1 vL1]=projection_convex(footL(2)+10*a_ortho(1),footL(3)+10*a_ortho(2),...
                                footL(2),footL(3),...
                                xvL,yvL);
[xvL2 yvL2 vL2]=projection_convex(footL(2)-10*a_ortho(1),footL(3)-10*a_ortho(2),...
                                footL(2),footL(3),...
                                xvL,yvL);
                            
[xvR1 yvR1 vR1]=projection_convex(footR(2)+10*a_ortho(1),footR(3)+10*a_ortho(2),...
                                footR(2),footR(3),...
                                xvR,yvR);
[xvR2 yvR2 vR2]=projection_convex(footR(2)-10*a_ortho(1),footR(3)-10*a_ortho(2),...
                                footR(2),footR(3),...
                                xvR,yvR);
                            
% vL1_norm=(xvL(1:end-1)-xvL1).^2+(yvL(1:end-1)-yvL1).^2;
% vL2_norm=(xvL(1:end-1)-xvL2).^2+(yvL(1:end-1)-yvL2).^2;
% vR1_norm=(xvR(1:end-1)-xvR1).^2+(yvR(1:end-1)-yvR1).^2;
% vR2_norm=(xvR(1:end-1)-xvR2).^2+(yvR(1:end-1)-yvR2).^2;
% 
% [B vL1]=sort(vL1_norm);
% [B vL2]=sort(vL2_norm);
% [B vR1]=sort(vR1_norm);
% [B vR2]=sort(vR2_norm);
vL1=vL1(2);
vL2=vL2(2);
vR1=vR1(2);
vR2=vR2(2);

if vL1>vL2
    xuL=xvL2;
    xlL=xvL1;

    yuL=yvL2;
    ylL=yvL1;
else
    xuL=xvL1;
    xlL=xvL2;

    yuL=yvL1;
    ylL=yvL2;
end


if vL1==1 || vL2==1
    switch (vL1*vL2)
         case 2
             xvL_reduced=[xuL;xvL(2:end-3);xlL];
             yvL_reduced=[yuL;yvL(2:end-3);ylL];
         case 3
             xvL_reduced=[xuL;xvL(2:end-2);xlL];
             yvL_reduced=[yuL;yvL(2:end-2);ylL];
         case 4
             xvL_reduced=[xuL;xvL(2:end-1);xlL];
             yvL_reduced=[yuL;yvL(2:end-1);ylL];
    end
elseif vL1==2 || vL2==2
    switch (vL1*vL2)
         case 6
             xvL_reduced=[xuL;xvL(3:end-2);xlL];
             yvL_reduced=[yuL;yvL(3:end-2);ylL];
         case 8
             xvL_reduced=[xuL;xvL(3:end-1);xlL];
             yvL_reduced=[yuL;yvL(3:end-1);ylL];
    end
else
    xvL_reduced=[xuL;xvL(4:end-1);xlL];
    yvL_reduced=[yuL;yvL(4:end-1);ylL];
end

if vR1<vR2
    xuR=xvR2;
    xlR=xvR1;

    yuR=yvR2;
    ylR=yvR1;
else
    xuR=xvR1;
    xlR=xvR2;

    yuR=yvR1;
    ylR=yvR2;
end


if vR1==1 || vR2==1
    switch (vR1*vR2)
         case 2
             xvR_reduced=[xlR;xvR(2:end-3);xuR];
             yvR_reduced=[ylR;yvR(2:end-3);yuR];
         case 3
             xvR_reduced=[xlR;xvR(2:end-2);xuR];
             yvR_reduced=[ylR;yvR(2:end-2);yuR];
         case 4
             xvR_reduced=[xlR;xvR(2:end-1);xuR];
             yvR_reduced=[ylR;yvR(2:end-1);yuR];
    end
elseif vR1==2 || vR2==2
    switch (vR1*vR2)
         case 6
             xvR_reduced=[xlR;xvR(3:end-2);xuR];
             yvR_reduced=[ylR;yvR(3:end-2);yuR];
         case 8
             xvR_reduced=[xlR;xvR(3:end-1);xuR];
             yvR_reduced=[ylR;yvR(3:end-1);yuR];
    end
else
    xvR_reduced=[xlR;xvR(4:end-1);xuR];
    yvR_reduced=[ylR;yvR(4:end-1);yuR];
end

% xv=[xlR;xuR;xuL;xlL];
% yv=[ylR;yuR;yuL;ylL];
xv=[xvR_reduced;xvL_reduced];
yv=[yvR_reduced;yvL_reduced];
xv = [xv ; xv(1)]; yv = [yv ; yv(1)];

figure(1);clf;axis equal;
hold on;plot(xvL,yvL);hold off
hold on;plot(xvR,yvR);hold off
hold on;plot(xv,yv,'red');hold off
hold on;plot(xvR_reduced,yvR_reduced,'green');hold off
hold on;plot(xvL_reduced,yvL_reduced,'green');hold off