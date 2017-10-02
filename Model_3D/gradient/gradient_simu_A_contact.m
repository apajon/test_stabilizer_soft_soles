function A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric,contAngle,displ_ini,psi_ini,angleact)
    A11 = -Wc;
    A12 = zeros(3*length(ind_cont),2*length(ind_slip));
    A21 = zeros(2*length(ind_slip),3*length(ind_cont));
    A22 = zeros(2*length(ind_slip),2*length(ind_slip));
    if ~isempty(ind_slip)
        I = [];
        for i=1:length(ind_slip)
            I(i) = find(ind_cont==ind_slip(i));
        end
        Fc_mc_splip = Fc_mc(:,I);
        for i=1:length(ind_slip)
            if contAngle==1
                Rpsiini = [cos(psi_ini), -sin(psi_ini), 0; sin(psi_ini), cos(psi_ini), 0; 0, 0, 1];
                theta = angleact(1);
                phi = angleact(2);
                psi = angleact(3);
                [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi);
                Rini = Rtheta * Rphi * Rpsiini;
                delta_t = P_mc(:,I(i)) - displ_ini - Rini * PabsOld_mc(:,I(i));
            else
                delta_t = P_mc(:,I(i))-PabsOld_mc(:,I(i));
            end
            delta_t = delta_t(1:2,1);
            A21((2*i)-1:2*i,(3*I(i))-2:3*I(i)) = [norm(delta_t).*eye(2) fric.*delta_t];
            A12((3*I(i))-2:3*I(i),(2*i)-1:2*i) = [1 0; 0 1; 0 0];
            A22((2*i)-1:2*i,(2*i)-1:2*i) =  Fc_mc_splip(1:2,i)*(delta_t'./norm(delta_t))+fric.*Fc_mc_splip(3,i).*eye(2);
        end
    end
    A = [A11 A12;A21 A22];
end