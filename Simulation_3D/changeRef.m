function [ZMP,F] = changeRef(ZMPx, ZMPy, Fx, Fy, Fz,backtoankle,fronttoankle,exttoankle,inttoankle,xpankle,ypankle,rightorleft)
%%%%%%%%%% Foot %%%%%%%%%%
% ---------------------
% |                    |
% |                    |   
% |                    |  0.13
% |                    |    
% ---------------------
%          0.23
% taille du pied :
% distance en l’arrière du pied et la cheville : 0.1
% distance en l’avant du pied et la cheville : 0.13
% distance en l’extérieur du pied et la cheville : 0.075 (bord gauche dans notre cas)
% distance en l’intérieur du pied et la cheville : 0.055 (bord droit dans notre cas)
% Je t’ai mis en plus le tracé des courbes.
% backtoankle=0.098; %from back to ankle of foot
% fronttoankle=0.128; %from  front to ankle of foot
% exttoankle=0.076; %from exterior to ankle of foot
% inttoankle=0.054; %from interior to ankle of foot
% xpankle=1.257728354176543; %x coordinate of ankle position
% ypankle=-0.045000000000611; %y coordinate of ankle position
%rightorleft : 1 for right and -1 for left
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lfooty = 0.13;
% Lfootx = 0.23;
Lfooty = exttoankle+inttoankle;
Lfootx = backtoankle+fronttoankle;
%%% translation of trajectory
% traslx = min(ZMPx)-0.023;
% trasly = Lfooty/2-0.01;
traslx=xpankle-backtoankle;
trasly=ypankle-rightorleft*(Lfooty/2-inttoankle);
% traslx=xpankle;
% trasly=ypankle;
%%% foot
% sidex1 = zeros(1,14);
% sidex2 = zeros(1,14) + Lfootx;
sidey = -(Lfooty/2):0.001:(Lfooty/2);
sidex1 = zeros(1,size(sidey,2));
sidex2 = zeros(1,size(sidey,2)) + Lfootx;
figure (3); clf; axis equal; 
hold on
plot(sidex1,sidey)
plot(sidex2,sidey)
% hold on
sidex = 0:0.001:Lfootx;
sidey1 = zeros(1,size(sidex,2)) - Lfooty/2;
sidey2 = zeros(1,size(sidex,2)) + Lfooty/2;
plot(sidex,sidey1)
% hold on
plot(sidex,sidey2)
% ZMPxNew = ZMPx-traslx;
% ZMPyNew = ZMPy+trasly;
ZMPxNew = ZMPx-traslx;
ZMPyNew = ZMPy-trasly;

plot(ZMPxNew,ZMPyNew,'-')
% plot(0.1,(Lfootx/2-(2*0.055)),'+')
plot(backtoankle,rightorleft*(Lfooty/2-inttoankle),'+')
hold off

figure(4);plot(Fz)

% % theta=pi/2;
% % R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
ZMP = [ZMPxNew';ZMPyNew'];
% %ZMP = R*ZMP;
% % R * ZMPxNew
% % figure; clf; plot(ZMP(1,:),ZMP(2,:))
% % hold on
% % plot(-(0.13/2-0.055),0.1,'+')
% % hold on
% % %%% foot
% % sidex1 = zeros(12) + 0.13/2;
% % sidex2 = zeros(12) - 0.13/2;
% % sidey = 0:0.02:0.23; 
% % plot(sidex1,sidey)
% % hold on
% % plot(sidex2,sidey)
% % [R,Rtheta,Rphi,Rpsi] = Rot(0,0,-pi/2);
% % F1 = [Fx';Fy';Fz'];
% % % Rtheta * 
F = [Fx';Fy';Fz'];
% % F = R*F1;
% % figure; clf; stem(F(1,:))
% % figure; clf; stem(F(2,:))
% % figure; clf; stem(F(3,:))