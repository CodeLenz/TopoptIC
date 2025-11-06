//
// Domínio quadrado 
//
//          φm
//     ------------
//     |          |
//  φm |          |  hn
//     |          | 
//     |__________|  
//          φm

// Element size
lc = 0.01;

// Corners
Point(1) = { 0, 0,   0, lc};
Point(2) = { 1, 0,   0, lc};
Point(3) = { 1, 1,   0, lc};
Point(4) = { 0, 1,   0, lc};

// Edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Material
Physical Surface("Material,mat,1,1E-8") = {1};

// Boundary conditions - Open
Physical Curve("φm") = {1,3,4};

// Boundary conditions - bn
Physical Curve("hn,1.0") = {2};

// Convert triangles to quads
Recombine Surface{:};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "teste1.msh";
