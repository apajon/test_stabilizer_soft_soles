function writemshfile(sole,pname)
elename = 'elements1.txt';
delete 'newmesh.msh';
fname = 'newmesh.msh';
coorNodesnew(:,1) = 1:sole.nTot;
coorNodesnew(:,2) = sole.coor(:,1);
coorNodesnew(:,3) = sole.coor(:,2);
coorNodesnew(:,4) = sole.coor(:,3);
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