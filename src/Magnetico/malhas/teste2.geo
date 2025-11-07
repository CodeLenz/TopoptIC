SetFactory("OpenCASCADE");

//
// Retângulo 1 x 1 sem termo fonte unitário com mu=1.0 em todo 
// o domínio
//
//          φm        
//     ------------
//     |          |
//  φm |          |  φm
//     |          | 
//     |__________|  
//          φm
//
//

// Tamanho do elemento 
lc = 0.01;

// Cantos
Point(1) = { 0, 0,   0, lc};
Point(2) = { 1, 0,   0, lc};
Point(3) = { 1, 1,   0, lc};
Point(4) = { 0, 1,   0, lc};

// Arestas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Loop para a superfície
Curve Loop(1) = {1, 2, 3, 4};

// Gera a superfície
Plane Surface(1) = {1};

// Material
Physical Surface("Material,mat,1,1.0") = {1};

// Condição de contorno essencial
Physical Curve("φm") = {1,2,3,4};

// Termo fonte 
Physical Surface("ρm,1.0") ={1};

// Convert triangles to quads
Recombine Surface{:};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "teste2.msh";
