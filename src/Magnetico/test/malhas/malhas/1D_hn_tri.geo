//
// Domínio retangular com dimensão L x H 
// L << H (1D)
//
// sem termo fonte e com hn = 1 na face da direita
// a face da esquerda tem potencial nulo prescrito
//
// A propriedade do material é \mu = 1 em todo o domínio
//
//        
//     ----------------------- ->
//     |                     | ->
//  φm |                     | -> hn
//     |                     | ->
//     |_____________________| -> 
//         
//

// Tamanho do elemento 
lc = 0.01;

// Dimensões 
L = 2.0;
H = 0.1;


// Cantos do domínio
Point(1) = { 0, 0,   0, lc};
Point(2) = { L, 0,   0, lc};
Point(3) = { L, H,   0, lc};
Point(4) = { 0, H,   0, lc};

// Arestas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Superfície
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Material
Physical Surface("Material,mat,1,1.0") = {1};

// Boundary conditions
Physical Curve("φm") = {4};

// Boundary conditions - hn
Physical Curve("hn,1.0") = {2};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "1D_hn_tri.msh";
