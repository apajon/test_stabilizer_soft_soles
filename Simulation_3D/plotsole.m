function plotsole(number,connecext,coor,stressVM,Pc,Fc,Z,Ftot)
figure(number)
clf();
patch('faces',connecext,'vertices',coor,'FaceVertexCData',stressVM(:,1),'FaceColor','interp');
% grid on
axis([-0.1  0.3 -0.24  0.24  -0.2  0.05])
%axis([-0.3  0.3 -0.3  0.3  -0.1  0.2])
view(3);
hold on;
% pointA = [-10,10,0];
% pointB = [-10,-10,0];
% pointC = [10,-10,0];pointD = [10,10,0];
% points = [pointA' pointB' pointC' pointD']; % using the data given in the question
% fill3(points(1,:),points(2,:),points(3,:),'r')
% alpha(1)
% hold on;
%view(90,0)

% [X,Y] = meshgrid(-0.3:.05:0.3);
% R = sqrt(X.^2 + Y.^2) + eps;
% Z = sin(0)./R;
% % patch(surf2patch(X,Y,Z,Z),'FaceColor','g','EdgeColor','flat',...
% %       'Marker','o','MarkerFaceColor','flat');
% colorPlane = zeros(size(X,1)*size(X,1),3);
% colorPlane(:,1) = 0.9;
% colorPlane(:,2) = 0.9;
% colorPlane(:,3) = 0.9;
% patch(surf2patch(X,Y,Z,Z),'FaceVertexCData',colorPlane);
% shading faceted;
% hold on;
%view(3)
%In case of contact
if ~isempty(Fc)%----------------------------------------------------
    alpha(1)
    % Plot the contact forces
    quiver3(Pc(:,1),Pc(:,2),Pc(:,3),-Fc(:,1)*0.005,-Fc(:,2)*0.005,-Fc(:,3)*0.005,0)

    %-----------------------calcul ZMP position--------------------------------
    %Fsol = sum(Fc,1);
    Mo = sum(cross(Pc',Fc'),2);
    Xzmp = Z(1);
    Yzmp = Z(2);
    Zzmp = 0;

    % Track the Zmp (green line)
    %quiver3(Xzmp,Yzmp,Zzmp,-Ftot(1)*0.0003,-Ftot(2)*0.0003,-Ftot(3)*0.0003,0)
    
    drawnow();
end
end
