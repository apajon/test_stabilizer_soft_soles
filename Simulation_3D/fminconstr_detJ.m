function [c,ceq] = fminconstr_detJ(p_mid,sole,spl,volini)
%%% New coor %%%
coornew = deformation(size(sole.coor),p_mid,spl);
%%% translate because spl is in coornew-sole.trasl
coornew(:,1)= coornew(:,1) + sole.trasl(1);
coornew(:,2) = coornew(:,2) + sole.trasl(2);
coornew(:,3) = coornew(:,3) + sole.trasl(3);

detJ = zeros(size(sole.elements_vol,1),1);
vol = 0;
for j = 1:size(sole.elements_vol,1)
    vertices = coornew(sole.elements_vol(j,:),:);
    detJ(j) = det([1,1,1,1;vertices']);
    vol = vol + abs(detJ(j)/6);
end
%volini = 8.969999999999988e-04;
volini = 8.970000000000016e-04;
%c = [vol-volini;-vol+(volini/2)];
c = [vol-volini;-detJ(j)];
figure(8); clf;
ite = 1;
stem(ite,vol,'b');
hold on;
stem(ite,volini,'r');
ceq = [];