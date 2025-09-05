//
// Problema 8 x 5 engastado com força embaixo
//

// Comprimento da viga
L = 8.0;

// Altura da viga
H = 5.0;

// Espessura 
b = 0.01;

// tamanho do elemento
lc = 0.2;

// Cantos
Point(1) = { 0, 0, 0};
Point(2) = { L, 0, 0};
Point(3) = { L, H, 0};
Point(4) = { 0, H, 0};

// Lados
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Superfície
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Material
Physical Surface("Material,aço,1,210E9,0.3,7850.0") = {1};

// Força concentrada no canto superior direito 
// direção Y 
Physical Point("Fc,2,100.0") = {2};

// Prende os gls x dos nós da face da esquerda
Physical Curve("U,1,0.0") = {4};

// Prende os gls y dos nós da face da esquerda
Physical Curve("U,2,0.0") = {4};

//
// Até aqui, geometria. Agora, geração de malha
//

// Geração da malha 
// https://gmsh.info/doc/texinfo/gmsh.html#Mesh-options

// Algoritmo para geração de malha
// 8 (OpenCASCADE):
// Leverages OpenCASCADE's meshing capabilities for 2D surfaces.
Mesh.Algorithm = 8;

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
Mesh.RecombinationAlgorithm = 0;

// Tamanho da malha
Mesh.CharacteristicLengthFactor = lc;

// Cria a malha
Mesh 2;

// Smooth the mesh
// OptimizeMesh "Laplace2D";
OptimizeMesh "Gmsh";

// Grava a malha
Save "8x5.msh";
