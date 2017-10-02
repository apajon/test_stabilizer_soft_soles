function [x0 y0 ii]=projection_convex(xp,yp,xi,yi,xv,yv)
%we search the projection of the point p=(xp,yp) in direction of i=(xi,yi)
%on the convex polygon define by the vertices (xv,yv)
%i has to be inside the convex polygon
%if p is iside the convex polygon, the function return p
if ~inpolygon(xi,yi,xv,yv)
    error(message('i=(xi,yi) is not inside the polygon=(xv,yv)'));
end

if inpolygon(xp,yp,xv,yv)
    x0=xp;
    y0=yp;
    ii=0;
else
%     a0=[xp yp]-[xi yi];
%     a0=a0(2)/a0(1);
%     [x0 y0]=search_projection_convex(xp,yp,xi,yi,a0,xv,yv);
    [x0 y0,ii]=polyxpoly([xp;xi],[yp;yi],xv,yv);
end

% function [x0 y0]=search_projection_convex(xp,yp,xi,yi,a,xv,yv)
% %we search the closest point to ZMP projected on the convex
% %hull edges along the perpendicular to ankle segment
% x0=[];y0=[];
% norm0=Inf;
% for j=1:length(xv)-1
%     a0=(yv(j+1)-yv(j))/(xv(j+1)-xv(j));
%     [x0_ y0_]=intersec_line(a0,xv(j),yv(j),a,xp,yp);
%     
%     v1=[x0_;y0_]-[xv(j);yv(j)];
%     v2=[x0_;y0_]-[xv(j+1);yv(j+1)];
%     if sign(v1(1)/v2(1))==-1 || sign(v1(2)/v2(2))==-1%inpolygon(x0_,y0_,xv,yv)
%         if norm0>norm([xp-x0_;yp-y0_])+norm([xi-x0_;yi-y0_])
%             norm0=norm([xp-x0_;yp-y0_])+norm([xi-x0_;yi-y0_]);
%             x0=x0_;
%             y0=y0_;
%         end
%     else
%         if norm([xv(j+1)-x0_;yv(j+1)-y0_])>norm([xv(j)-x0_;yv(j)-y0_])
%             x0_=xv(j);
%             y0_=yv(j);
%         else
%             x0_=xv(j+1);
%             y0_=yv(j+1);
%         end
%         if norm0>norm([xp-x0_;yp-y0_])+norm([xi-x0_;yi-y0_])
%             norm0=norm([xp-x0_;yp-y0_])+norm([xi-x0_;yi-y0_]);
%             x0=x0_;
%             y0=y0_;
%         end
%     end
% end