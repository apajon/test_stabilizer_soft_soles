function [displ,angleact,Fc_mc3,PcSurf,ind_Cont,ind_Cont3,Ftot,Z,Kcart,no_conv,J,Mo] = PsoleAngle(m,Ccc,sole,P0,Pg,fric,contAngle,FcSurf,displ,angleact,Fdes,Zdes)
D = 3;
no_conv = [];
if contAngle==1
    Pfree = P0(:,sole.nodesFreeSurf);
%     PabsOld = zeros(D,sole.nFreeSurf);
    PabsOld(1,:) = P0(1,sole.nodesFreeSurf);%-0.1;
    PabsOld(2,:) = P0(2,sole.nodesFreeSurf);%-0.1;
    PabsOld(3,:) = P0(3,sole.nodesFreeSurf);    
    FcNew3 = zeros(D*m,1);
else
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld = Pg(sole.nodesFreeSurf,:)';
    FcNew = FcSurf'; 
    %FcNew3 = zeros(D*m,1);
    FcNew3 = reshape(FcNew,D*m,1);
end
diplini = displ;
angleactini = angleact;

% if contAngle==1
%     FtotZMPdes = [Fdes;Zdes];    
%     angleact = [angleact(1);angleact(2)];
%     [displ,angleact,FtotZMP,Fc_mc3_out,ind_cont_out,Kcart] = GaussFtotZMP2Angles(contAngle,fric,m,displ,angleact,FtotZMPdes,FcNew3,Pfree,PabsOld,Ccc);    angleact = [angleact(1);angleact(2);0];
%     PcSurf = [];
% else
    %%% 3 angles %%%
    FtotZMPdes = [Fdes;Zdes;0];    
    [displ,angleact,FtotZMP,Fc_mc3_out,PcSurf,ind_cont_out,Kcart,J] = GaussFtotZMP(contAngle,fric,m,diplini,angleactini,FtotZMPdes,FcNew3,Pfree,PabsOld,Ccc);
%     disp(['Mzmp_z ' num2str(FtotZMP(6))])
    %%% 2 angles %%%
%     FtotZMPdes = [Fdes;Zdes];    
%     angleact = [angleact(1);angleact(2)];
%     [displ,angleact,FtotZMP,Fc_mc3_out,ind_cont_out,Kcart]  = GaussFtotZMP2Angles(contAngle,fric,m,displ,angleact,FtotZMPdes,FcNew3,Pfree,PabsOld,Ccc);
%     angleact = [angleact(1);angleact(2);0];
%     PcSurf = [];    
% end

if norm(displ)==0
    displ = diplini;
    angleact = angleactini;
    a = load('Kcart.mat');
    Kcart = a.Kcart;
    no_conv = contAngle;
end
if ind_cont_out(1)==0
    a = find(ind_cont_out==0,2);
    Fc_mc3 = Fc_mc3_out(1:D*(a(2)-1));
    ind_Cont = (sort(ind_cont_out(1:(a(2)-1))+1))';        
else
    a = find(ind_cont_out==0,1);
    Fc_mc3 = Fc_mc3_out(1:D*(a-1));
    ind_Cont = (sort(ind_cont_out(1:(a-1))+1))';              
end

ind_Cont3 = sort([D*ind_Cont-2 D*ind_Cont-1 D*ind_Cont]);
Ftot = FtotZMP(1:3);
Z = FtotZMP(4:5);
Mo = FtotZMP(6);

end

