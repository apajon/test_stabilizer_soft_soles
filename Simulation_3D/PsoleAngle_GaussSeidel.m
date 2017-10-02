function [displ,angleact,Fc_mc3,ind_Cont,ind_Cont3,Ftot,Z,Kcart] = PsoleAngle_GaussSeidel(m,Ccc,sole,P0,Pg,fric,contAngle,FcSurf,displ,angleact,Fdes,Zdes)
D = 3;
if contAngle==1
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld(1,:) = P0(1,sole.nodesFreeSurf);%-0.1;
    PabsOld(2,:) = P0(2,sole.nodesFreeSurf);%-0.1;
    PabsOld(3,:) = P0(3,sole.nodesFreeSurf);    
    FcNew3 = zeros(D*m,1);
else
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld = Pg(sole.nodesFreeSurf,:)';
    FcNew = FcSurf'; 
    FcNew3 = reshape(FcNew,D*m,1);
end

FtotZMPdes = [Fdes;Zdes];

[displ,angleact,FtotZMP,Fc_mc3_out,ind_cont_out] = Gauss(contAngle,fric,m,displ,angleact,[FtotZMPdes;0],FcNew3,Pfree,PabsOld,Ccc);

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
Kcart = [];

end

