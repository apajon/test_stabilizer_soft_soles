function coornew = remesh(sole,pname,coornew)
elename = 'elements1.txt';
%delete 'newmesh.msh';
delete 'new.msh';
fname = 'newmesh.msh';
coorNodesnew(:,1) = 1:sole.nTot;
coorNodesnew(:,2) = coornew(:,1);
coorNodesnew(:,3) = coornew(:,2);
coorNodesnew(:,4) = coornew(:,3);
coorNodesnew = coorNodesnew';
fid = fopen(fname,'w');
fprintf(fid,'$MeshFormat \n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n',sole.nTot);
fprintf(fid,'%d %f %f %f\n',coorNodesnew);
fprintf(fid,'$EndNodes\n');
fclose(fid);

fid = fopen(fname,'a');
fmid = fopen([pname elename],'r');

for i=1:(sole.nEle+3)
    tline = fgets(fmid);
    %h=sscanf(tline)';
    fprintf(fid,tline);
end

fclose(fid);
fclose(fmid);
% command = 'cd C:\Users\Giovanni\Documents\Giappone\Software\gmsh-svn-Windows64\gmsh-2.8.5-svn-Windows';
% status = dos(command)
%%% I put gmesh in the work directory, in other way it does not works and I
%%% don't know why
% command = 'C:\Users\Giovanni\Documents\Giappone\Software\gmsh-2.8.5-Windows\gmsh -refine C:\Users\Giovanni\Documents\Giappone\Code\sim3d\Simulation 3D - Desired ZMP - Position and Force - 3 Angles - meters\input\semelle1 L=0.23, l=0.13, e=0.04 m new centre\newmesh.msh';
command = 'gmsh newmesh.msh -0 -3 refineMesh optimize -o new.msh';

% status = dos(command,'-echo');
[status,cmdout] = dos(command);

fname = 'new.msh';
% %pname = '/'; 
% % pname = 'input/semelle1 L=100, l=40, e=20 rought/'; 
% % pname = 'input/semelle1 L=0.239, l=0.139, e=0.04 m/';
pname = '/';
[elements_surf,elements_vol,coornew] = input_mesh(pname,fname);     % reads nodes and elem. from GMSH file
% sole = soleFEM_newStiff(pname,fname);