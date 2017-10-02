lc=6;
l=0.139;
L=0.239;
e=0.02;

Point(1) = {0, 0, 0, lc};
Point(2) = {l, 0, 0, lc};
Point(3) = {l, 0, e, lc};
Point(4) = {0, 0, e, lc};
Point(5) = {0, L, 0, lc};
Point(6) = {l, L, 0, lc};
Point(7) = {l, L, e, lc};
Point(8) = {0, L, e, lc};


//Physical Point(2)={2};
//Physical Point(3)={3};
//Physical Point(6)={6};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {-5, -6, -7, -8};
Line Loop(3) = {-2, 10, 6, -11};
Line Loop(4) = {11, 7, -12, -3};
Line Loop(5) = {-4, 12, 8, -9};
Line Loop(6) = {-1, 9, 5, -10};


Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Physical Surface(30) = {3};
Physical Surface(50) = {5};
Physical Surface(10) = {1};
Physical Surface(40) = {4};
Physical Surface(20) = {2};
Physical Surface(60) = {6};


Surface Loop(100) = {1,2,3,4,5,6};

Volume(100)={100};

Physical Volume(100) = {100};

Field[1] = Cylinder;
Field[2] = Threshold;
