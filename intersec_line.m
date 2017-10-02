function [x y]=intersec_line(a1,x1,y1,a2,x2,y2)
%give the intersection point of two straight line
% assuming we have the components of lines functions as:
%(L1) : y=a1(x-x1)+y1
%(L2) : y=a2(x-x2)+y2
%if L1 and L2 are parallel, return (x2,y2) coordinate
if a1==a2
    x=x2;
    y=y2;
elseif isinf(a1)
    x=x1;
    y=a2*(x-x2)+y2;
elseif isinf(a2)
    x=x2;
    y=a1*(x-x1)+y1;
else
    x=(a1*x1-a2*x2-y1+y2)/(a1-a2);
    y=a1*(x-x1)+y1;
end
