function [c, ceq] = confun(x,Fz_d,xp_smooth,yp_smooth,Fz1_d,xp1_d,yp1_d,Fz2_d,xp2_d,yp2_d)
% Nonlinear inequality constraints
c = [0];
% Nonlinear equality constraints
% ceq = [-Fz_d*yp_smooth+(Fz1_d+x(1))*(yp1_d+x(3))+(Fz2_d+x(2))*(yp2_d+x(4))];
ceq = [-Fz_d*xp_smooth+(Fz1_d+x(1))*(xp1_d+x(3))+(Fz2_d+x(2))*(xp2_d+x(4)); ...
        -Fz_d*yp_smooth+(Fz1_d+x(1))*(yp1_d+x(5))+(Fz2_d+x(2))*(yp2_d+x(6)); ...
        ((xp1_d+x(3))-xp_smooth)*((yp2_d+x(6))-yp_smooth)-((xp2_d+x(4))-xp_smooth)*((yp1_d+x(5))-yp_smooth)];
end