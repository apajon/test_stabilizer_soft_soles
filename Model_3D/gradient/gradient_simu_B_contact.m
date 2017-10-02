function B = gradient_simu_B_contact(angleact,P_mc3,Fc_mc3,ind_cont,Ccc,ind_slip,psi_first,contAngle)
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);    
    [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi);
    for i=1:length(ind_cont)
        if contAngle==1
            Rpsiini = [cos(psi_first), -sin(psi_first), 0; sin(psi_first), cos(psi_first), 0; 0, 0, 1];            
            dRinidphi = Rtheta*dR_dphi*Rpsiini;
            dRinidtheta = dR_dtheta*Rphi*Rpsiini;
            btheta = dR_dtheta*P_mc3((3*i-2):3*i);
            bphi = dR_dphi*P_mc3((3*i-2):3*i);
            bpsi = dR_dpsi*P_mc3((3*i-2):3*i);
            btheta(1) = btheta(1) - dRinidtheta(1,:)*P_mc3((3*i-2):3*i);
            btheta(2) = btheta(2) - dRinidtheta(2,:)*P_mc3((3*i-2):3*i);
            bphi(1) = bphi(1) - dRinidphi(1,:)*P_mc3((3*i-2):3*i);
            bphi(2) = bphi(2) - dRinidphi(2,:)*P_mc3((3*i-2):3*i);
        else
            btheta = dR_dtheta*P_mc3((3*i-2):3*i);
            bphi = dR_dphi*P_mc3((3*i-2):3*i);
            bpsi = dR_dpsi*P_mc3((3*i-2):3*i);
        end
        for j=1:length(ind_cont)
            btheta = btheta + (dR_dtheta * Ccc((3*i-2):3*i,(3*j-2):3*j) * R' + R * Ccc((3*i-2):3*i,(3*j-2):3*j) * dR_dtheta') * Fc_mc3((3*j-2):3*j);
            bphi = bphi + (dR_dphi * Ccc((3*i-2):3*i,(3*j-2):3*j) * R' + R * Ccc((3*i-2):3*i,(3*j-2):3*j) * dR_dphi') * Fc_mc3((3*j-2):3*j);
            bpsi = bpsi + (dR_dpsi * Ccc((3*i-2):3*i,(3*j-2):3*j) * R' + R * Ccc((3*i-2):3*i,(3*j-2):3*j) * dR_dpsi') * Fc_mc3((3*j-2):3*j);
        end
        B12((3*i-2):3*i,:) = [btheta bphi bpsi];
    end
    B11 = repmat(eye(3),length(ind_cont),1);
    B = [B11 B12;zeros(2*length(ind_slip),3) zeros(2*length(ind_slip),3)];
end