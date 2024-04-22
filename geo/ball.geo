SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMax = 0.1;
dim=3;

Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
Physical Surface("Gamma_D") = {1};
Physical Volume("Omega") = {1};