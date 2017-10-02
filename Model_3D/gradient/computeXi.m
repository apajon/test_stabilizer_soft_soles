function Xi = computeXi(displ,Zdes)
Xi = [1, 0, 0,      0,                  displ(3),           -displ(2)+Zdes(2);
      0, 1, 0,      -displ(3),          0,                  displ(1)-Zdes(1);
      0, 0, 1,      displ(2)-Zdes(2),  -displ(1)+Zdes(1), 0;
      0, 0, 0, 1, 0, 0;
      0, 0, 0, 0, 1, 0;
      0, 0, 0, 0, 0, 1];
end