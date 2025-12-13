//
// Domínio retangular com dimensão L x H 
// L << H (1D)
//
// com  termo fonte 
// a face da esquerda tem potencial nulo prescrito
//
// A propriedade do material é \mu = 1 em todo o domínio
//
//
// φ(x) = -ρm*x² / 2 + ρm*L*x
// 
// φ(0) = 0
// φ(2) = -2ρm + 4ρm = 2ρm
//        
//     ----------------------- 
//     |                     | 
//  φm |       ρm            | 
//     |                     | 
//     |_____________________| 
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

// Termo fonte 
Physical Surface("ρm,1.0") = {1};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "1D_rhom_tri.msh";
