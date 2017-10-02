function Lambda = gradient_simu_Lambda_contact(ind_cont,der_C_cc,ind_slip,Fc_mc3,RBlock_c)
    der_W_c = zeros(3*length(ind_cont),3*length(ind_cont));
%     FcOld3_mc(3:3:end) = 0;
    for i=1:3*length(ind_cont)
        der_W_c(:,i) = RBlock_c * der_C_cc{i} * RBlock_c' * Fc_mc3;
    end
    %RI = RBlock * repmat([1 0 0;0 1 0;0 0 1],length(ind_cont),length(ind_cont));
    %dPold = (RI_old + der_W_c_old);
%     dPold(3:3:end,:)=0;
%     dPold(:,3:3:end)=0;
    
    Lambda = [der_W_c;zeros(2*length(ind_slip),3*length(ind_cont))];
end