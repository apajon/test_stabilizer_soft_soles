function [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi)
Rtheta = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Rphi = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
Rpsi = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R = Rtheta*Rphi*Rpsi;