//
// Problema 8 x 5 engastado com força no meio da face lateral
//

// Comprimento da viga
L = 8.0;

// Altura da viga
H = 5.0;

// Espessura 
b = 0.01;

//
// Pontos
//
// Vamos dividir por 3 seções horizontais, para podermos 
// aplicar força distribuída no meio da face..
//
Point(1) = { 0, 0, 0};
Point(2) = { L, 0, 0};

Point(3) = { 0, H/3, 0};
Point(4) = { L, H/3, 0};

Point(5) = { 0, 2*H/3, 0};
Point(6) = { L, 2*H/3, 0};

Point(7) = { 0, H, 0};
Point(8) = { L, H, 0};

// Cerca cada pedaço com as suas linhas 
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Line(5) = {4, 6};
Line(6) = {6, 5};
Line(7) = {5, 3};

Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {7, 5};


// Divisões usando transfinito
nx = 10;
ny = 3;
Transfinite Line {1} = nx;
Transfinite Line {2} = ny;
Transfinite Line {3} = nx;
Transfinite Line {4} = ny;
Transfinite Line {5} = ny;
Transfinite Line {6} = nx;
Transfinite Line {7} = ny;
Transfinite Line {8} = ny;
Transfinite Line {9} = nx;
Transfinite Line {10} = ny;

// Agora define as superfícies
// Superfícies
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {-3, 5, 6, 7};
Curve Loop(3) = {-6, 8, 9, 10};

Plane Surface(100) = {1};
Plane Surface(200) = {2};
Plane Surface(300) = {3};


// E converte para transfinitas. Os dados são os pontos dos cantos
Transfinite Surface {100} = {1, 2, 4, 3};
Recombine Surface {100};

Transfinite Surface {200} = {3, 4, 6, 5};
Recombine Surface {200};

Transfinite Surface {300} = {5, 6, 8, 7};
Recombine Surface {300};

// Material
Physical Surface("Material,aço,1,210E9,0.3,7850.0") = {100,200,300};

// Força concentrada no meio do lado da direita, direção Y
Physical Curve("Ft,2,10000.0") = {5};

// Prende os gls x dos nós da face da esquerda
Physical Curve("U,1,0.0") = {4,7,10};

// Prende os gls y dos nós da face da esquerda
Physical Curve("U,2,0.0") = {4,7,10};

//
// Até aqui, geometria. Agora, geração de malha
//

// Geração da malha 
// https://gmsh.info/doc/texinfo/gmsh.html#Mesh-options

// Algoritmo para geração de malha
// 8 (OpenCASCADE):
// Leverages OpenCASCADE's meshing capabilities for 2D surfaces.
Mesh.Algorithm = 8;

// 0=none, 1=all quadrangles, 2=all hexahedra.
Mesh.SubdivisionAlgorithm = 1;

// Opção global para transformar os triângulos em retângulos
Mesh.RecombineAll = 1;

// Algoritmo de recombinação de malha
// 0: Simple:
//    This is a basic recombination algorithm that attempts to create quads or hexes from triangles. 
// 1: Blossom: 
//    This algorithm uses a minimum cost perfect matching algorithm to generate fully quadrilateral meshes from triangulations according to a paper published in 2011. 
// 2: Simple Full-Quad:
//    This option aims to produce a fully quadrilateral mesh with a simpler approach than Blossom. 
// 3: Blossom Full-Quad:
//    Combines the Blossom approach with the goal of generating a fully quadrilateral mesh. 
Mesh.RecombinationAlgorithm = 2;

// Cria a malha
Mesh 2;

// Smooth the mesh
// OptimizeMesh "Laplace2D";
// OptimizeMesh "Gmsh";

// Grava a malha
Save "8x5-distrib-meio_transfinite.msh";
