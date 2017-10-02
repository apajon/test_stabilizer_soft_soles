function plotsoleAnkle(number,connecext,coor,stressVM,Pc,Fc,Z,Ftot,Ankle)
figure(number)
clf();
patch('faces',connecext,'vertices',coor,'FaceVertexCData',stressVM(:,1),'FaceColor','interp');
grid on
axis([-0.1  0.3 -0.24  0.24  -0.1  0.2])
view(3);
hold on;
plot3(Ankle(1),Ankle(2),Ankle(3),'+r');
text(Ankle(1),Ankle(2),Ankle(3),'Ankle','VerticalAlignment','bottom', 'HorizontalAlignment','right');
hold on;
%view(90,0)
%In case of contact
if ~isempty(Fc)%----------------------------------------------------
    % Plot the contact forces
    quiver3(Pc(:,1),Pc(:,2),Pc(:,3),Fc(:,1)*0.25,Fc(:,2)*0.25,Fc(:,3)*0.25,0)

    %-----------------------calcul ZMP position--------------------------------
    %Fsol = sum(Fc,1);
    Mo = sum(cross(Pc',Fc'),2);
    Xzmp = Z(1);
    Yzmp = Z(2);
    Zzmp = 0;

    % Track the Zmp (green line)
    quiver3(Xzmp,Yzmp,Zzmp,Ftot(1)*0.2,Ftot(2)*0.2,Ftot(3)*0.2,0)
    
    drawnow();
end
end
