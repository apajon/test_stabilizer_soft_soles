function [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angle)
theta = angle(1);
phi = angle(2);
psi = angle(3);
Rtheta = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Rphi = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
Rpsi = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R = Rtheta*Rphi*Rpsi;
dR_dtheta = [0 0 0; 0 -sin(theta) -cos(theta); 0 cos(theta) -sin(theta)]* Rphi * Rpsi;
dR_dphi = Rtheta * [-sin(phi) 0 cos(phi); 0 0 0; -cos(phi) 0 -sin(phi)] * Rpsi;
dR_dpsi = Rtheta * Rphi * [-sin(psi) -cos(psi) 0; cos(psi) -sin(psi) 0; 0 0 0];
end