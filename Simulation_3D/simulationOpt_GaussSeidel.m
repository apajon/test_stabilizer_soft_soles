function cost=simulationOpt_GaussSeidel(sole,coornew,friction,Fdes,Zdes,connecext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Stiffness                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /*Position of heel compare to ankle position 
% -0.0749569
% 0
% -0.031668
% */
% /*Position of halx compare to ankle position 
% 0.174899
% 0
% -0.031668
% */
eSole = 0.06;
Lankle = [-0.0749569;0.;-0.031668];
% stressVM0 = zeros(sole.analysis.NN,1);
% plotsole(1,sole.connecext,coornew,stressVM0,[],[],[],[]);
% 
% /* Dimension of the foot
%  * length = 0.24985640 m
%  * widht = 0.1914 m
%  * e = 0.031668 m
%  */
%e = 0.031668;

% displ_in=0.755568
% 0
% -0.0512123
% 
% angleact_in=0
% 0.510589
% 0
coornew(:,1)= coornew(:,1) + Lankle(1);
coornew(:,2) = coornew(:,2) + Lankle(2);
coornew(:,3) = coornew(:,3)  + Lankle(3) -eSole;

stressVM0 = zeros(sole.analysis.NN,1);
% plotsole(1,sole.connecext,coornew,stressVM0,[],[],[],[]);

% angleact = [0.;0.5;0.];
% R = Rot(angleact(1),angleact(2),angleact(3));
% coornew = (R*coornew')';
% plotsole(4,sole.connecext,coornew,stressVM0,[],[],[],[]);

Ankle = [0; 0; 0];
% coornew(:,1)= coornew(:,1);
% coornew(:,2) = coornew(:,2);
% coornew(:,3) = coornew(:,3)+0.12;
% 
% 
% 
% plotsole(5,sole.connecext,coornew,stressVM0,[],[],[],[]);
% PosAnkle = [-0.0749569;0.0;0.1713+efoot];
% 
% angleact = [degtorad(0);1.1216;degtorad(0)];
% R = Rot(angleact(1),angleact(2),angleact(3));
% stressVM0 = zeros(sole.analysis.NN,1);
% plotsole(4,sole.connecext, (R*coornew')',stressVM0,[],[],[],[]);
% coornew = (R*coornew')';
% PosAnklevec = repmat(PosAnkle, 1, sole.nTot);
% coornew = coornew' + PosAnklevec;
% stressVM0 = zeros(sole.analysis.NN,1);
% plotsole(5,sole.connecext,coornew',stressVM0,[],[],[],[]);
% 
% 


x_center = 0;
y_center = 0;
z_center = 0;
center = repmat([x_center; y_center; z_center], 1, sole.nTot);
% % center = repmat([-0.0749569;0.0;0.1713-efoot], 1, sole.nTot);
Pg = [];
FcSurf = [];
P0 = coornew' - center;
% angleact = [degtorad(0);1.1216;degtorad(0)];
% R = Rot(angleact(1),angleact(2),angleact(3));
% plotsole(3,sole.connecext, (R*coornew')',stressVM0,[],[],[],[]);

angle_X_degres_ini = 2 + zeros(1,11);
angle_X_degres_middle = (2.:-0.2:-2.);
angle_X_degres_end = -2 + zeros(1,11);
angle_X_degres = [angle_X_degres_ini angle_X_degres_middle angle_X_degres_end];

angle_X = degtorad(angle_X_degres);

angle_Y_degres_ini = -2 + zeros(1,11);
angle_Y_degres_middle1 = (-2:0.2:0);
angle_Y_degres_middle2 = (-0.2:-0.2:-2);
angle_Y_degres_end = -2 + zeros(1,11);
angle_Y_degres = [angle_Y_degres_ini angle_Y_degres_middle1 angle_Y_degres_middle2 angle_Y_degres_end];

angle_Y = degtorad(angle_Y_degres);
% angle_Y = zeros(size(angle_X_degres));
% 
% angle_Z_degres_ini = 2 + zeros(1,11);
% angle_Z_degres_middle = (2:-0.2:-2);
% angle_Z_degres_end = -2 + zeros(1,11);
% angle_Z_degres = [angle_Z_degres_ini angle_Z_degres_middle angle_Z_degres_end];
% 
% angle_Z = degtorad(angle_Z_degres);
angle_Z = zeros(size(angle_X_degres));

center = repmat([x_center; y_center; z_center], 1, sole.nTot);
displZ_ini = (0:-0.001:-0.01)-0.03;
displZ_middle = (-0.03:0.00145:-0.001)-0.01;
displZ_end = (-0.01:0.003:0.0205)+0.0005;
displZ = [displZ_ini displZ_middle displZ_end];
displZ = zeros(size(angle_Z));
%DispY = 0:0.27:22.3;
% DispX = 0:0.27:22.3;
DispX = zeros(size(displZ));
DispY = zeros(size(displZ));
% disp('prep stress comput')
% ABe = prep_stressVonMises(sole);
% disp('finished')

% displ = [0;0;-3];
% angleact = [0.17;0;0];
MatForce=zeros(3,length(angle_X));
OZ = zeros(2,length(angle_X));
for i=1:(length(angle_X))
%     displ = [DispX(i);DispY(i);displZ(i)];
%     angleact = [angle_X(i);angle_Y(i);angle_Z(i)];
%     displ = [-0.0749569;0.0;0.1713-efoot];
    displ = [0.;0.;0.12];
    angleact = [0.;0.5;0.];
    [Pg,Pc,dP,contact,FcSurf,Ftot,Z,displ,angleact,Kcart] = contacts_GaussSeidel(sole,P0,center,displ,Pg,friction,angleact,i,FcSurf,Fdes(:,1),Zdes(:,1));
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    Ankle = R*(Ankle)+displ;    
    Fc = FcSurf(contact,:);
%     stressVM = stressVonMises(sole,dP,sole.nodesDirichlet3,ABe); % VonMises' Stress
%     pause(3);
    MatForce(:,i) = Ftot;
    OZ(:,i) = [Z(1);Z(2)];

    plotsoleAnkle(2,sole.connecext,Pg,stressVM0,Pc,Fc,Z,Ftot,Ankle)
    Vector = [Z(1)-Ankle(1);Z(2)-Ankle(2);0-Ankle(3)]
    cost = 1;
    %Kcart = StiffCart(angle,contact,sole,Pc);
%     if i==1
%         cost = 0;
%     else
%         if i<(length(Fdes)-26)
             %cost = cost + norm(Kcart(1:3,1:3),2);
%         end
%         cost=cost-norm(Kcart(4:6,4:6),2);
%     end
end
figure; clf;plot(OZ(1,:),OZ(2,:))
hold on
plot(-(0.13/2-0.055),0.1,'+')
hold on
%%% foot
sidex1 = zeros(12) + 0.13/2;
sidex2 = zeros(12) - 0.13/2;
sidey = 0:0.02:0.23; 
plot(sidex1,sidey)
hold on
plot(sidex2,sidey)
end
