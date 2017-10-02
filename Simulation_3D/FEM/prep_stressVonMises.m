function ABe = prep_stressVonMises(sole)
%--------------------------------------------------------------------------
%         Computation of the VonMises' Stress
%--------------------------------------------------------------------------
mu=sole.E/(2*(1+sole.nu));
lambda=sole.E*sole.nu/((1+sole.nu)*(1-2*sole.nu));
for j = 1:size(sole.elements_vol,1)                         % computes nodal stresses
    vertices = sole.coor(sole.elements_vol(j,:),:);
    PhiGrad = [1,1,1,1;vertices']\[zeros(1,3);eye(3)];
    Be = zeros(6,12);
    Be([1,4,5],1:3:10) = PhiGrad';
    Be([4,2,6],2:3:11) = PhiGrad';
    Be([5,6,3],3:3:12) = PhiGrad';
    C(1:3,1:3) = lambda*ones(3,3)+2*mu*eye(3);
    C(4:6,4:6) = mu*eye(3);

    ABe{j} = C*Be;    
end
