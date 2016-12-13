Point(1) = {0, 0, 0, 500};
Point(2) = {3600, 0, 0, 500};
Point(3) = {4800, 0, 0, 500};
Point(4) = {6000, 0, 0, 500};
Point(5) = {13800, 0, 0, 500};
Point(6) = {13800, 1000, 0, 500};
Point(7) = {6000, 1000, 0, 500};
Point(8) = {4800, 1000, 0, 500};
Point(9) = {3600, 1000, 0, 500};
Point(10) = {0, 1000, 0, 500};

Line(1) = {10, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {2, 9};
Line(12) = {3, 8};
Line(13) = {4, 7};

Physical Line(20) = {1};
Physical Line(21) = {2, 3, 4, 5};
Physical Line(22) = {6};
Physical Line(23) = {7, 8, 9, 10};

Line Loop(24) = {10, 1, 2, 11};
Plane Surface(25) = {24};
Line Loop(26) = {11, -9, -12, -3};
Plane Surface(27) = {26};
Line Loop(28) = {12, -8, -13, -4};
Plane Surface(29) = {28};
Line Loop(30) = {13, -7, -6, -5};
Plane Surface(31) = {30};
Physical Surface(32) = {25, 27, 29, 31};
