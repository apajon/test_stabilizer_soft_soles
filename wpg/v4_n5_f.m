function [v4]=v4_n5_f(ti,i)
v4=[((-1).*ti(i-2).^3.*ti(i+1).^2.*((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i-2).^2.*ti(i-1).* ...
  ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*( ...
  (-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i-2).^2.*ti(i).*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*ti(i+2).^2.*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^3.*ti(i+1).* ...
  ((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^2.*ti(i).*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^( ...
  -1)+(-1).*ti(i-2).*ti(i-1).^2.*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).* ...
  ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^( ...
  -1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i+4).^2.*(ti(i+4)+( ...
  -1).*ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^(-1).*ti(i).^3.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  ti(i-1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)) ...
  .^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-1).*ti(i-2).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i)+ti(i+2)).^(-1).*ti(i+3).^2.*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).* ...
  ti(i)+ti(i+3)).^(-1)) (2.*ti(i-2).^3.*ti(i+1).*((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+3.*ti(i-2).^2.* ...
  ti(i+1).^2.*((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*(( ...
  -1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+ti(i-2).^2.*ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^( ...
  -1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i-2)+ti(i+3)).^(-1)+ti(i-2).^2.*ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).* ...
  ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+ ...
  ti(i-2).^2.*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+2.*ti(i-2).*ti(i-1).*ti(i+1).*((-1).* ...
  ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+2.*ti(i-2).^2.*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*(( ...
  -1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ ...
  ti(i+3)).^(-1)+ti(i-2).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).^2.*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1) ...
  .*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+2.*ti(i-2).*ti(i).*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).^2.*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  ti(i-1).^3.*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i-1)+ti(i+3)).^(-1)+3.*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^2.*ti(i+1).*((-1).*ti(i-1)+ti(i+1)) ...
  .^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+( ...
  ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^3.*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  ti(i-1).^2.*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^2.*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ ...
  2.*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).* ...
  ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(ti(i+4)+(-1).*ti(i-1)).^( ...
  -1).*ti(i-1).^2.*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*ti(i-1).^2.*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)) ...
  .^(-1)+ti(i-2).*ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).* ...
  ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*ti(i-1).^2.*((-1) ...
  .*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+2.*ti(i-2).*ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1)+ti(i-1).^2.*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*ti(i-1).*ti(i).*( ...
  (-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1) ...
  .*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*(( ...
  -1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1) ...
  .*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-1).*ti(i).*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*(( ...
  -1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+3.*ti(i+4).^2.*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  (ti(i+4)+(-1).*ti(i)).^(-1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*( ...
  (-1).*ti(i)+ti(i+3)).^(-1)+2.*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^(-1).* ...
  ti(i).^3.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+ ...
  ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+2.*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)) ...
  .^(-1).*ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i).^2.*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ ...
  ti(i+3)).^(-1)+(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).* ...
  ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+2.*ti(i-2).* ...
  ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+2.*ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).^2.*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^( ...
  -1).*((-1).*ti(i)+ti(i+3)).^(-1)+ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1) ...
  .*ti(i+3).^2.*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)) (( ...
  -1).*ti(i-2).^3.*((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-6).*ti(i-2).^2.*ti(i+1).*((-1).* ...
  ti(i-2)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^( ...
  -1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-3).*ti(i-2).*ti(i+1).^2.*((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^( ...
  -1)+(-1).*ti(i-2).^2.*ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i-2).^2.*ti(i+1).* ...
  ((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-2).*ti(i-2).*ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1) ...
  .*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).* ...
  ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i-2).^2.*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).* ...
  ((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-2).*ti(i-2).* ...
  ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*( ...
  (-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-2).*ti(i-2).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^( ...
  -1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).* ...
  ((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+( ...
  -1).*ti(i-2).^2.*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-2).*ti(i-2).^2.*((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-4).*ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1) ...
  .*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)) ...
  .^(-1)+(-2).*ti(i-2).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).^2.*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1) ...
  .*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i).*((-1) ...
  .*ti(i)+ti(i+1)).^(-1).*ti(i+2).^2.*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).* ...
  ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-3).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  ti(i-1).^2.*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^3.*((-1).*ti(i-1)+ti(i+1)).^(-1) ...
  .*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-3).* ...
  ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-3).*(ti(i+4)+(-1).*ti(i-1)).^(-1) ...
  .*ti(i-1).^2.*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^2.*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^( ...
  -1)+(-2).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1) ...
  .*ti(i-1)).^(-1).*ti(i-1).^2.*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1) ...
  .*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-2).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*( ...
  (-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^2.*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+( ...
  -1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-2).*(ti(i+4)+(-1).*ti(i-1)) ...
  .^(-1).*ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).*ti(i-1).^2.*((-1).*ti(i-1)+ti(i+1)).^(-1).*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1)+(-2).*ti(i-2).*ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*(( ...
  -1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).* ...
  ti(i-1).^2.*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).*ti(i-1).*ti(i).*((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1) ...
  +(-1).*ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-1).*ti(i).*((-1) ...
  .*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-2).*ti(i-2).*ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1)+(-1).*ti(i-1).^2.*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).* ...
  ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*( ...
  (-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-2).*ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^( ...
  -1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).* ...
  ((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^( ...
  -1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+( ...
  -1).*ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^( ...
  -1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-1).*ti(i).*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1) ...
  .*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1)+(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1) ...
  .*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i).*( ...
  (-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*(( ...
  -1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-3).*ti(i+4).^2.*(ti(i+4)+(-1).*ti(i-1)).^( ...
  -1).*(ti(i+4)+(-1).*ti(i)).^(-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i)+ti(i+3)).^(-1)+(-6).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^( ...
  -1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^( ...
  -1)+(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^(-1).*ti(i).^3.*((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-2).*ti(i+4).*(ti(i+4)+(-1) ...
  .*ti(i-1)).^(-1).*ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*(( ...
  -1).*ti(i)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+( ...
  -1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)) ...
  .^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-2).*ti(i+4).*(ti(i+4)+( ...
  -1).*ti(i-1)).^(-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1) ...
  .*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-2).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i).* ...
  ((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1) ...
  .*ti(i)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*( ...
  (-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-1).* ...
  ti(i-2).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-4).*ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^( ...
  -1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-2).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)) ...
  .^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^( ...
  -1)+(-1).*ti(i-2).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).^2.*((-1).* ...
  ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-2).*ti(i).*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).^2.*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)) (3.*ti(i-2).^2.*((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^( ...
  -1)+6.*ti(i-2).*ti(i+1).*((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+ti(i+1).^2.*((-1).*ti(i-2)+ti(i+1)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i-2)+ti(i+3)).^(-1)+ti(i-2).^2.*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+2.*ti(i-2).*ti(i-1).*(( ...
  -1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+2.*ti(i-2).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^( ...
  -1)+ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+2.*ti(i-2).*((-1).*ti(i-1)+ti(i+1)).^(-1) ...
  .*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i-2)+ti(i+3)).^(-1)+ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*(( ...
  -1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+ti(i+1).*((-1).* ...
  ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+ti(i-2).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+ ...
  2.*ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+4.*ti(i-2).*((-1).*ti(i)+ti(i+1)).^(-1).* ...
  ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1) ...
  .*ti(i-2)+ti(i+3)).^(-1)+2.*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*ti(i+2).^2.*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+3.*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).* ...
  ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^( ...
  -1)+3.*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^2.*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^( ...
  -1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).* ...
  ((-1).*ti(i-1)+ti(i+3)).^(-1)+3.*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1) ...
  .*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+2.*ti(i+4).* ...
  (ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).^2.*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1)+ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+2.*(ti(i+4)+(-1).*ti(i-1)).^( ...
  -1).*ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1) ...
  .*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+2.*( ...
  ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*( ...
  (-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i).*((-1) ...
  .*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1)+2.*ti(i-2).*ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-1).^2.*((-1).*ti(i-1)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).* ...
  ((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+2.*ti(i-1).* ...
  ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1) ...
  .*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^( ...
  -1)+ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^( ...
  -1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*(( ...
  -1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-1).*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1) ...
  +ti(i-2).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).* ...
  ((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+2.*ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+3)).^(-1)+ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)) ...
  .^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-2).*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*( ...
  (-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+ti(i).*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).* ...
  ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^( ...
  -1)+ti(i+4).^2.*(ti(i+4)+(-1).*ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^(-1).*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+6.*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^( ...
  -1).*(ti(i+4)+(-1).*ti(i)).^(-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i)+ti(i+3)).^(-1)+3.*(ti(i+4)+(-1).*ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^(-1).* ...
  ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+ ...
  ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+2.*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1) ...
  .*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*(( ...
  -1).*ti(i)+ti(i+3)).^(-1)+2.*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).* ...
  ((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(ti(i+4)+(-1) ...
  .*ti(i-1)).^(-1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^( ...
  -1)+(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).* ...
  ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+2.*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*( ...
  (-1).*ti(i)+ti(i+3)).^(-1)+2.*ti(i-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).* ...
  ((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+ti(i).^2.*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+2.*ti(i-2).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)) ...
  .^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^( ...
  -1)+4.*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+((-1).*ti(i)+ti(i+1)).^(-1).*(( ...
  -1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).^2.*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*(( ...
  -1).*ti(i)+ti(i+3)).^(-1)) ((-3).*ti(i-2).*((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*(( ...
  -1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-2).*ti(i+1).*(( ...
  -1).*ti(i-2)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-2).*ti(i-2).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1) ...
  +(-1).*ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1) ...
  .*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).* ...
  ti(i-2)+ti(i+3)).^(-1)+(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).* ...
  ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-2).*ti(i-2).*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^( ...
  -1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-2) ...
  .*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+(-1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).* ...
  ((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+3)).^(-1)+(-3).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1) ...
  .*ti(i-1)).^(-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1) ...
  +(-2).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^( ...
  -1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1) ...
  .*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*( ...
  (-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).*( ...
  (-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-2).*ti(i-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1) ...
  +(-1).*ti(i+1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-2).*((-1).*ti(i)+ti(i+1)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).* ...
  ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*ti(i).*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*ti(i+2).*((-1).*ti(i-1)+ti(i+2)).^( ...
  -1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1) ...
  .*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*ti(i+3).*(( ...
  -1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1)+(-2).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^(-1).*((-1).* ...
  ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-3).*(ti(i+4)+(-1).* ...
  ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^(-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-1).*ti(i+4).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1) ...
  +(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*ti(i-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-2).*(ti(i+4)+(-1).*ti(i-1)).^(-1) ...
  .*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*(( ...
  -1).*ti(i)+ti(i+3)).^(-1)+(-1).*(ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1) ...
  .*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-1).*ti(i-2).*( ...
  (-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ ...
  ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+(-2).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ ...
  ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1) ...
  +(-2).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*ti(i+3).*((-1).*ti(i-2)+ti(i+3)).^(-1) ...
  .*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)) (((-1).*ti(i-2)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^( ...
  -1)+((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)).^(-1).*((-1) ...
  .*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-2)+ti(i+2)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1)+( ...
  ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).* ...
  ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+((-1) ...
  .*ti(i-1)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)) ...
  .^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i-1)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1)+(ti(i+4)+(-1).* ...
  ti(i-1)).^(-1).*(ti(i+4)+(-1).*ti(i)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^( ...
  -1).*((-1).*ti(i)+ti(i+3)).^(-1)+(ti(i+4)+(-1).*ti(i-1)).^(-1).*((-1).*ti(i)+ti(i+1)).^(-1).*(( ...
  -1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)+((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i-2)+ti(i+3)).^(-1).*((-1).*ti(i-1)+ti(i+3)).^(-1) ...
  .*((-1).*ti(i)+ti(i+3)).^(-1))];
end