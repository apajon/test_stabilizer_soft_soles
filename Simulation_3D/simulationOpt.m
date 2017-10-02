function cost=simulationOpt(sole,friction,Fdes,Zdes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Stiffness                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cost = 0;
cost_lin = 0;
cost_rot = 0;

stressVM0 = zeros(sole.nTot,1);
% plotsole(1,sole.elements_surf,sole.coor,stressVM0,[],[],[],[]);

x_center = 0;
y_center = 0;
z_center = 0;
center = repmat([x_center; y_center; z_center], 1, sole.nTot);
Pg = [];
FcFreeSurf = [];
P0 = sole.coor' - center;

% disp('prep stress comput')
%ABe = prep_stressVonMises(sole);
% disp('finished')

% displ = [0.;0.;0];
angleact = [0.;0.;0];
% angleact = [sign(Fdes(2,1))*0.0000001;-sign(Fdes(1,1))*0.0000001;0.];
% displ = 1e-05 * [0.247743419624653;-0.309586755764408;-0.008896995069604];
% angleact = [0.17;0.;0.];
%angleact = [-0.034906585039887;0.034906585039887;0];
% angleact = 1e-03 * [-0.305682622062386;-0.203397104826419;-0.005494340064746];
% angleact = 1e-03 * [-sign(Zdes(2,1))*0.305682622062386;sign(Zdes(1,1)-max(sole.coor(sole.nodesDirichlet,1))/2)*0.203397104826419;-sign(Zdes(2,1))*0.005494340064746];
% angleact = [-sign(Zdes(2,1))*0.0001  sign(Zdes(1,1)-max(sole.coor(sole.nodesDirichlet,1))/2)*0.0001   -sign(Zdes(2,1))*0.0001]';
% displ = [sign(Zdes(2,1))*0.000001;sign(Zdes(1,1)-max(sole.coor(sole.nodesDirichlet,1))/2)*0.000001;-sign(Zdes(2,1))*0.000001];
displ = [-sign(Fdes(1,1))*0.0000001;-sign(Fdes(2,1))*0.0000001;-0.0000001];
% displ = [-sign(Fdes(1,1))*0.0001;-sign(Fdes(2,1))*0.0001;-0.0001];
% displ = [0;0;-0.005];

%angleact = [-0.034906585039887;0.034906585039887;0];
%angleact = [0.049243281471969;-0.056533641644207;0];
%angleact = [0.054906585039887;0.034906585039887;0];
angleact_tot = zeros(3,size(Fdes,2));
displ_tot = zeros(3,size(Fdes,2));
% if sole.IteOpt~=3
% a = load('angleact_tot.mat');
% angleact_tot = a.angleact_tot;
% a = load('displ_tot.mat');
% displ_tot = a.displ_tot;        
% end
cost_rot_tot = zeros(1,size(Fdes,2));
cost_lin_tot = zeros(1,size(Fdes,2));
no_conv_tot = [];
cont_no_conv = 0;

fet1=[];
fet2=[];
fet3=[];
fet4=[];
tpp=[];

% fet1_=(sum(sole.coor(:,:)==repmat([0    0.0650    0.0150],[size(sole.coor(:,:),1),1]),2)==3);
% fet2_=(sum(sole.coor(:,:)==repmat([0    -0.0650    0.0150],[size(sole.coor(:,:),1),1]),2)==3);
% fet3_=(sum(sole.coor(:,:)==repmat([0.2300    0.0650    0.0150],[size(sole.coor(:,:),1),1]),2)==3);
% fet4_=(sum(sole.coor(:,:)==repmat([0.2300    -0.0650    0.0150],[size(sole.coor(:,:),1),1]),2)==3);
fet1_=(sum(sole.coor(:,:)==repmat([0    0.0650    0.021],[size(sole.coor(:,:),1),1]),2)==3);
fet2_=(sum(sole.coor(:,:)==repmat([0    -0.0650    0.021],[size(sole.coor(:,:),1),1]),2)==3);
fet3_=(sum(sole.coor(:,:)==repmat([0.2300    0.0650    0.021],[size(sole.coor(:,:),1),1]),2)==3);
fet4_=(sum(sole.coor(:,:)==repmat([0.2300    -0.0650    0.021],[size(sole.coor(:,:),1),1]),2)==3);
% fet1_=(sum(sole.coor(:,:)==repmat([0    0.0650    0.006],[size(sole.coor(:,:),1),1]),2)==3);
% fet2_=(sum(sole.coor(:,:)==repmat([0    -0.0650    0.006],[size(sole.coor(:,:),1),1]),2)==3);
% fet3_=(sum(sole.coor(:,:)==repmat([0.2300    0.0650    0.006],[size(sole.coor(:,:),1),1]),2)==3);
% fet4_=(sum(sole.coor(:,:)==repmat([0.2300    -0.0650    0.006],[size(sole.coor(:,:),1),1]),2)==3);


angleact_init=angleact;
Pg_init=Pg;
FcFreeSurf_init=FcFreeSurf;
i=1;

% while i<=(size(Fdes,2))
%         i
%     %disp = [0;DispY(i);displZ(i)];
%     %disp(['Step ' num2str(i)])
% %     if sole.IteOpt~=3
% %         angleact = angleact_tot(:,i);
% %         displ = displ_tot(:,i);
% %     end
%     [Pg,contact,FcFreeSurf,Ftot,Z,displ,angleact,Kcart,no_conv,dP] = contacts(sole,P0,center,displ,Pg,friction,angleact,i,FcFreeSurf,Fdes(:,i),Zdes(:,i));
%     Fc = FcFreeSurf(contact,:);
%     Pc = Pg(sole.nodesFreeSurf(contact),:);
% 
%     %stressVM = stressVonMises(sole,dP,ABe); % VonMises' Stress
%     if ~isempty(no_conv)
%         cont_no_conv = cont_no_conv + 1;
%         no_conv_tot(cont_no_conv) = no_conv;
%     end
% 
% %     figure(3)
% %     hold on
% %     plot(Zdes(1,i),Zdes(2,i),'o')
% %     hold off
%     
%     plotsole(2,sole.elements_surf,Pg,stressVM0,Pc,Fc,Z,Ftot)
%     
%     fet1=[fet1;Pg(sole.nodesDirichlet(1),:)];
%     fet2=[fet2;Pg(sole.nodesDirichlet(2),:)];
% %     fet3=[fet3;Pg(sole.nodesDirichlet(7),:)];
% %     fet4=[fet4;Pg(sole.nodesDirichlet(8),:)];
%     fet3=[fet3;Pg(sole.nodesDirichlet(3),:)];
%     fet4=[fet4;Pg(sole.nodesDirichlet(4),:)];
% 
%     tpp=[tpp;angleact(1:3)'];
% 
%     if i<(length(Fdes)/10) %minimize linear stiffness for the first 10% of the step
%         %disp('min linear stiffness')
%         cost_lin = cost_lin + norm(Kcart(1:3,1:3),2);
%     end
%     %disp('max rotational stiffness')
%     %%% 04/02/2015 
%     % cost_rot = cost_rot - sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6))));
%     cost_rot = cost_rot + sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6)))); % maximize the minimum stiffness. Be carefull, this is not differentiable when changing of minimum eigenvalue !!!!!!
%     %cost_rot = cost_rot + norm(Kcart(4:6,4:6),2); % maximize the minimum stiffness
%     cost_rot_tot(i) = sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6))));
%     cost_lin_tot(i) = norm(Kcart(1:3,1:3),2);
%     angleact_tot(:,i)=angleact;
%     displ_tot(:,i)=displ;
% %     angleact
%         i=i+1;
%     if i>=3
% % %         angleact_tot(:,1)==angleact_tot(:,2) & angleact_tot(:,1)~=angleact_tot(:,i-1)
% %         if angleact_tot(:,1)==angleact_tot(:,2) & angleact_tot(:,1)~=angleact_tot(:,i-1)
% % %             i
% % %             angleact_tot(:,1)==angleact_tot(:,2) & angleact_tot(:,1)~=angleact_tot(:,i-1)
% %             angleact=angleact_tot(:,i-1);
% %             displ = displ_tot(:,i-1);
% %             i=1;
% % %             displ = [0.;0.;0];
% %             Pg=Pg_init;
% %             FcFreeSurf=FcFreeSurf_init;
% %             fet1=[];
% %             fet2=[];
% %             fet3=[];
% %             fet4=[];
% %             tpp=[];
% %         elseif angleact_tot(:,1)~=angleact_tot(:,2)
% %             i=size(Fdes,2)+1;
% %             angleact=angleact_tot(:,2);
% %             displ = displ_tot(:,i-1);
% % %             i=1;
% % %             displ = [0.;0.;0];
% %             Pg=Pg_init;
% %             FcFreeSurf=FcFreeSurf_init;
% %             fet1=[];
% %             fet2=[];
% %             fet3=[];
% %             fet4=[];
% %             tpp=[];
% %         end
%         if sum(angleact_tot(:,1)==angleact_init) && sum(angleact_tot(:,1)~=angleact_tot(:,i-1))
%             angleact=angleact_tot(:,i-1);
%             angleact_init=angleact;
%             displ = displ_tot(:,i-1);
%             i=1;
%             Pg=Pg_init;
%             FcFreeSurf=FcFreeSurf_init;
%             fet1=[];
%             fet2=[];
%             fet3=[];
%             fet4=[];
%             tpp=[];
%         end
%         if(angleact_tot(:,1)~=angleact_tot(:,2))
%             i=size(Fdes,2)+1;
%         end
%     end
% end

% angleact=angleact_tot(:,2);
Kcart_tot=[];
J_tot=[];
Ftot_tot=[];
Z_tot=[];
Mo_tot=[];
for i=1:(size(Fdes,2))
    if (mod(i/(size(Fdes,2)),0.1)>mod((i+1)/(size(Fdes,2)),0.1))
%         {'Foot step' i 'from' (size(Fdes,2))}
        {'Achieve' round(100*(i+1)/(size(Fdes,2))) '%'}
    end
%     i
    %disp = [0;DispY(i);displZ(i)];
    %disp(['Step ' num2str(i)])
%     if sole.IteOpt~=3
%         angleact = angleact_tot(:,i);
%         displ = displ_tot(:,i);
%     end
    [Pg,contact,FcFreeSurf,Ftot,Z,displ,angleact,Kcart,no_conv,dP,J,Mo] = contacts(sole,P0,center,displ,Pg,friction,angleact,i,FcFreeSurf,Fdes(:,i),Zdes(:,i));
    Fc = FcFreeSurf(contact,:);
    Pc = Pg(sole.nodesFreeSurf(contact),:);

    %stressVM = stressVonMises(sole,dP,ABe); % VonMises' Stress
    if ~isempty(no_conv)
        cont_no_conv = cont_no_conv + 1;
        no_conv_tot(cont_no_conv) = no_conv;
    end

%     figure(3)
%     hold on
%     plot(Zdes(1,i),Zdes(2,i),'o')
%     plot(Z(1),Z(2),'*')
%     hold off
%     sqrt((Zdes(1,i)-Z(1))^2+(Zdes(2,i)-Z(2))^2)
%     
%     plotsole(2,sole.elements_surf,Pg,stressVM0,Pc,Fc,Z,Ftot)
    
    fet1=[fet1;Pg(fet1_,:)];
    fet2=[fet2;Pg(fet2_,:)];
    fet3=[fet3;Pg(fet3_,:)];
    fet4=[fet4;Pg(fet4_,:)];

    tpp=[tpp;angleact(1:3)'];

    if i<(length(Fdes)/10) %minimize linear stiffness for the first 10% of the step
        %disp('min linear stiffness')
        cost_lin = cost_lin + norm(Kcart(1:3,1:3),2);
    end
    %disp('max rotational stiffness')
    %%% 04/02/2015 
    % cost_rot = cost_rot - sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6))));
    cost_rot = cost_rot + sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6)))); % maximize the minimum stiffness. Be carefull, this is not differentiable when changing of minimum eigenvalue !!!!!!
    %cost_rot = cost_rot + norm(Kcart(4:6,4:6),2); % maximize the minimum stiffness
    cost_rot_tot(i) = sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6))));
    cost_lin_tot(i) = norm(Kcart(1:3,1:3),2);
    angleact_tot(:,i)=angleact;
    displ_tot(:,i)=displ;
    Kcart_tot=[Kcart_tot;Kcart];
    J_tot=[J_tot;J];
    Ftot_tot=[Ftot_tot Ftot];
    Z_tot=[Z_tot Z];
    Mo_tot=[Mo_tot Mo];
end

% detJ = zeros(size(sole.elements_vol,1),1);
% vol = 0;
% for j = 1:size(sole.elements_vol,1)
%     vertices = sole.coor(sole.elements_vol(j,:),:);
%     detJ(j) = det([1,1,1,1;vertices']);
%     vol = vol + abs(detJ(j)/6);
% end
% volini = 8.969999999999988e-04;
% c = vol-volini;
% penalties = 0;
% if c > 0
%     penalties = 1000*vol;
% end

for i=1:size(no_conv_tot)
    if (no_conv_tot(i) ~= 1 && no_conv_tot(i) ~= length(Fdes))
        cost_rot_tot(no_conv_tot(i))=(cost_rot_tot(no_conv_tot(i)-1)+cost_rot_tot(no_conv_tot(i)+1))/2;
        cost_lin_tot(no_conv_tot(i))=(cost_lin_tot(no_conv_tot(i)-1)+cost_lin_tot(no_conv_tot(i)+1))/2;
    end
end
cost_lin_fin = 0;
% for i=1:(length(Fdes)/10)
for i=1:(length(Fdes)/10)
    cost_lin_fin = cost_lin_fin + cost_lin_tot(i);
end
% cost = cost_lin_fin - sum(cost_rot_tot);

save('Simulation_3D/results/trajectory.mat','fet1','fet2','fet3','fet4','tpp')
save('Simulation_3D/results/J_matrix.mat','J_tot','Kcart_tot','angleact_tot','displ_tot','Ftot_tot','Z_tot','Mo_tot')
% 
% 
% save 'results/cost_lin_fin.mat' cost_lin_tot;
% save 'results/cost_rot_fin.mat' cost_rot_tot;
% % if sole.IteOpt==3
% %      save 'angleact_tot.mat' angleact_tot;
% %      save 'displ_tot.mat' displ_tot;
% % end
% a = load('results/cost_rot_fin.mat');
% cost_rot_tot_fin = a.cost_rot_tot;
% a = load('results/cost_rot_ini.mat');
% cost_rot_tot_ini = a.cost_rot_tot;
% figure(6); clf;
% plot(cost_rot_tot_fin,'b');
% hold on;
% plot(cost_rot_tot_ini,'r')
% 
% a = load('results/cost_lin_fin.mat');
% cost_lin_tot_fin = a.cost_lin_tot;
% a = load('results/cost_lin_ini.mat');
% cost_lin_tot_ini = a.cost_lin_tot;
% figure(7); clf;
% plot(cost_lin_tot_fin,'b');
% hold on;
% plot(cost_lin_tot_ini,'r')
end
