function soleini = new_mesh(soleini,pname,coornew)
elename = 'elements1.txt';
delete 'newmesh.msh';
delete 'new.msh';
fname = 'newmesh.msh';
coorNodesnew(:,1) = 1:soleini.nTot;
coorNodesnew(:,2) = coornew(:,1);
coorNodesnew(:,3) = coornew(:,2);
coorNodesnew(:,4) = coornew(:,3);
coorNodesnew = coorNodesnew';
fid = fopen(fname,'w');
fprintf(fid,'$MeshFormat \n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n',soleini.nTot);
fprintf(fid,'%d %d %d %d\n',coorNodesnew);
fprintf(fid,'$EndNodes\n');
fclose(fid);

fid = fopen(fname,'a');
fmid = fopen([pname elename],'r');

for i=1:(soleini.nEle+3)
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
command = 'gmsh newmesh.msh -0 -3 optimize_ho -o new.msh';

% status = dos(command,'-echo');
status = dos(command);
% !C:/Users/Giovanni/Documents/Giappone/Software/gmsh-svn-Windows64/gmsh-2.8.5-svn-Windows/gmsh &
% fmid=fopen([pname fname],'r'); 
% tline = fgets(fmid);                        % sur la nouvelle version de GMSH, sauter 5l
% tline = fgets(fmid);                        % ancienne, sauter 2l
% %
% tline = fgets(fmid);
% tline = fgets(fmid);
% tline = fgets(fmid);
% %
% analysis.NN=sscanf(tline,'%d');
% temp=zeros(analysis.NN,1);
% coorNodes=zeros(analysis.NN,3);                   % matrix of cooordinates
% for i=1:analysis.NN,
%   tline = fgets(fmid);
%   h=sscanf(tline,'%d %f %f %f')';
%   temp(i)=h(1);
%   coorNodes(i,:)=h(2:length(h));
% end
fname = 'new.inp';
% %pname = '/'; 
% % pname = 'input/semelle1 L=100, l=40, e=20 rought/'; 
% % pname = 'input/semelle1 L=0.239, l=0.139, e=0.04 m/';
pname = '/';
soleini = soleFEM_newStiff(pname,fname);

% fname = 'new.msh';
% pname = '/';
% [analysis,elements,corrispnod,coorNodes]= input_mesh(pname,fname);  
% [connec,analysis]=input_solid(elements,analysis);
% [analysis,TD]=input_TD(elements,analysis);
% connecext = zeros(analysis.nTD,3);
% for i = 1:analysis.nTD
%    connecext(i,:) = TD(i).nodes;
% end
% Xmax = max(coorNodes(:,1));
% Xmin = min(coorNodes(:,1));
% Ymax = max(coorNodes(:,2));
% Ymin = min(coorNodes(:,2));
% Zmax = max(coorNodes(:,3));
% Zmin = min(coorNodes(:,3));
% nodesDirichlet = find(coorNodes(:,3)==Zmax);

