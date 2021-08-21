// Gmsh project created on Sun Aug 08 12:22:16 2021
//+
Point(1) = {0, 0, -4000, 1.0};
//+
Point(2) = {-3000, 0, 0, 1.0};
//+
Point(3) = {0, 0, 4000, 1.0};
//+
Point(4) = {0, 5000, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {2, 4};
//+
Line(3) = {3, 4};
//+
Physical Point("fixed", 4) = {1, 2, 3};
//+
Physical Point("force", 5) = {4};
//+
Physical Curve("matprop1", 6) = {1, 3};
//+
Physical Curve("matprop2", 7) = {2};
