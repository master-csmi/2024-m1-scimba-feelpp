SetFactory("OpenCASCADE");
h=0.1;
dim=3;

Box(1) = {0, 0, 0, 1, 1, 1};
Characteristic Length{ PointsOf{ Volume{1}; } } = h;
Physical Surface("Gamma_D") = {1,2,3,4,5,6};
Physical Volume("Omega") = {1};