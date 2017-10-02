%%% Spline Class
classdef SplineClass<handle
    properties
        idx_u
        idx_v
        basisValues
        u
        v
        r
        spline_res
        l_spli
    end
    methods
        function [obj]=SplineClass(sole,spline_res1)
            obj.spline_res = spline_res1;
            obj.l_spli = length(obj.spline_res);
            coor = zeros(sole.nTot,3);
            coor(:,1)= sole.coor(:,1) - sole.trasl(1);
            coor(:,2) = sole.coor(:,2) - sole.trasl(2);
            coor(:,3) = sole.coor(:,3) - sole.trasl(3);
            azimuth = zeros(1,size(coor,1));
            elevation = zeros(1,size(coor,1));
            r1 = zeros(1,size(coor,1));
            for i=1:size(sole.coor,1);
                [azimuth(1,i),elevation(1,i),r1(1,i)] = cart2sph(-coor(i,3),coor(i,2),coor(i,1));
            end
            spl1 = splineBasis([azimuth;elevation;r1], obj.spline_res, obj.spline_res);
            %%% Convert spline to Matlab Class
            splCell = struct2cell(spl1);
            obj.idx_u = splCell(1,:);
            obj.idx_v = splCell(2,:);
            obj.basisValues = splCell(3,:);
            obj.u = splCell(4,:);
            obj.v = splCell(5,:);
            obj.r = splCell(6,:);
        end
        function [obj]=update(obj,sole)            
            coor = zeros(sole.nTot,3);
            coor(:,1)= sole.coor(:,1) - sole.trasl(1);
            coor(:,2) = sole.coor(:,2) - sole.trasl(2);
            coor(:,3) = sole.coor(:,3) - sole.trasl(3);
            azimuth = zeros(1,size(coor,1));
            elevation = zeros(1,size(coor,1));
            r1 = zeros(1,size(coor,1));
            for i=1:size(sole.coor,1);
                [azimuth(1,i),elevation(1,i),r1(1,i)] = cart2sph(-coor(i,3),coor(i,2),coor(i,1));
            end
            spl1 = splineBasis([azimuth;elevation;r1], obj.spline_res, obj.spline_res);
            %%% Convert spline to Matlab Class
            splCell = struct2cell(spl1);
            obj.idx_u = splCell(1,:);
            obj.idx_v = splCell(2,:);
            obj.basisValues = splCell(3,:);
            obj.u = splCell(4,:);
            obj.v = splCell(5,:);
            obj.r = splCell(6,:);           
        end        
    end
end