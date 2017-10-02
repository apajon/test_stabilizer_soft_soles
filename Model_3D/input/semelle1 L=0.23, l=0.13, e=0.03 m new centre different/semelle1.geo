lc=0.025;
lc1=0.005;
lc2=0.025;
l=0.13;
L=0.23;
e=0.03;

Point(1) = {0, -l/2, 0, lc1};
Point(2) = {0, l/2, 0, lc1};
Point(3) = {0, l/2, e, lc};
Point(4) = {0, -l/2, e, lc};
Point(5) = {L, -l/2, 0, lc1};
Point(6) = {L, l/2, 0, lc1};
Point(7) = {L, l/2, e, lc};
Point(8) = {L, -l/2, e, lc};
Point(9) = {0, 0, 0, lc2};
Point(10) = {L/2, l/2, 0, lc2};
Point(11) = {L/2, -l/2, 0, lc2};
Point(12) = {L, 0, 0, lc2};
//Physical Point(2)={2};
//Physical Point(3)={3};
//Physical Point(6)={6};


Line(1) = {1, 9};
Line(2) = {9, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 1};
Line(6) = {5, 12};
Line(7) = {12, 6};
Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 5};
Line(11) = {1, 11};
Line(12) = {11, 5};
Line(13) = {2, 10};
Line(14) = {10, 6};
Line(15) = {3, 7};
Line(16) = {4, 8};

Line Loop(1) = {1, 2, 3, 4, 5};
Line Loop(2) = {-6, -7, -8, -9, -10};
Line Loop(3) = {-3, 13, 14, 8, -15};
Line Loop(4) = {15, 9, -16, -4};
Line Loop(5) = {-5, 16, 10, -11, -12};
Line Loop(6) = {-1, -2, 11, 12, 6, 7, -13, -14};


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