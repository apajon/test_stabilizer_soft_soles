function [elements_surf,elements_vol,coor]=input_mesh(pname,mname)

  fmid=fopen([pname mname],'r');
  % Jump 4 lines %
  tline = fgets(fmid);
  tline = fgets(fmid);
  tline = fgets(fmid);
  tline = fgets(fmid);
  %%%%%%%%%%%%%%%
  tline = fgets(fmid);
  nNodes = sscanf(tline,'%d');
  
  temp = zeros(nNodes,1);
  coor = zeros(nNodes,3);                   % matrix of cooordinates
  for i=1:nNodes
      tline = fgets(fmid);
      h=sscanf(tline,'%d %f %f %f')';
      temp(i)=h(1);
      coor(i,:)=h(2:length(h));
  end
  
  % Jump 2 lines %
  tline = fgets(fmid);
  tline = fgets(fmid);
  %%%%%%%%%%%%%%%
  tline = fgets(fmid);
  nElem = sscanf(tline,'%d');
  cont_surf = 0;
  cont_vol = 0;
  for i=1:nElem
      tline = fgets(fmid);
      h=sscanf(tline,'%d')';
      %elements(i,:)=h(6:length(h));  %ancienne:6:fin | nouvelle : 7:fin
      if length(h)<9
          cont_surf = cont_surf + 1;
          elements_surf(cont_surf,:)=h(6:length(h));
      else
          cont_vol = cont_vol + 1;
          elements_vol(cont_vol,:)=h(6:length(h));
      end
  end
  
  fclose(fmid);
end

