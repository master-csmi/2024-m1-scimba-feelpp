SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMax = 0.1;
dim=2;

Disk(1) = {0, 0, 0, 1.0};
Physical Curve("Gamma_D") = {1};
Physical Surface("Omega") = {1};