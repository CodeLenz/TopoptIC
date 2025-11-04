SetFactory("OpenCASCADE");

//
// Domínio quadrado 
//
//          ϕ0
//     ------------
//     |          |
//  ϕ0 |     0    |  ϕ0
//     |          | 
//     |__________|  
//          ϕ0

// Element size
lc = 0.01;

// Tamanho de elemento 
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// Corners
Point(1) = { 0, 0,   0};
Point(2) = { 1, 0,   0};
Point(3) = { 1, 1,   0};
Point(4) = { 0, 1,   0};

// Edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Circulo no meio 
Circle(5) = {0.5, 0.5, 0, 0.1, 0, 2*Pi};

// Loop para a superfície
Curve Loop(1) = {1, 2, 3, 4};

// Curve loop para o círculo 
Curve Loop(2) = {5};

// Gera a superfície
Plane Surface(1) = {1};

// Incorpora o círculo na malha
Curve{5} In Surface{1};

// Material
Physical Surface("Material,mat,1,1E-8") = {1};

// Boundary conditions - Open
Physical Curve("Open") = {1,2,3,4};

// Boundary conditions - bn
Physical Curve("bn,1.0") = {5};

// Convert triangles to quads
Recombine Surface{:};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "teste2.msh";
