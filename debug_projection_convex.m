close all
figure()
xv = [wpg_param.step_number_pankle_fixed(2,2)-wpg_param.backtoankle; ...
        wpg_param.step_number_pankle_fixed(2,2)+wpg_param.fronttoankle;...
        wpg_param.step_number_pankle_fixed(2,2)+wpg_param.fronttoankle; ...
        wpg_param.step_number_pankle_fixed(2,2)-wpg_param.backtoankle;...
        ];
yv = [wpg_param.step_number_pankle_fixed(2,3)-wpg_param.inttoankle; ...
        wpg_param.step_number_pankle_fixed(2,3)-wpg_param.inttoankle;...
        wpg_param.step_number_pankle_fixed(2,3)+wpg_param.exttoankle; ...
        wpg_param.step_number_pankle_fixed(2,3)+wpg_param.exttoankle;...
        ];
xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
hold on;
plot(xv,yv)
plot(wpg_param.step_number_pankle_fixed(2,2),wpg_param.step_number_pankle_fixed(2,3),'+')
hold off
for i=0:12 
    toto=[i/50+wpg_param.step_number_pankle_fixed(2,2) 0];
%     toto=[wpg_param.step_number_pankle_fixed(2,2)+wpg_param.fronttoankle+0.1 i/50+wpg_param.step_number_pankle_fixed(2,3)];
    [x0 y0]=projection_convex(toto(1),toto(2),wpg_param.step_number_pankle_fixed(2,2),wpg_param.step_number_pankle_fixed(2,3),xv,yv);
   
    hold on;
    plot(toto(1),toto(2),'o')
    plot(x0,y0,'ro')
    plot([toto(1);x0],[toto(2);y0],'k')
    hold off
end