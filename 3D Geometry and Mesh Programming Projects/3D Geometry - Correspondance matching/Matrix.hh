#ifndef __MATRIX_HEADER__
#define __MATRIX_HEADER__

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <iostream>
#include "MeshViewer.hh"


/* THis class defines general matrix and matrix operations */

typedef double Matrix3x3[3][3];
typedef double* pMatrix3x3[3][3];
typedef double Matrix9x4[9][4];
typedef double Matrix9x1[9];

void MatrixInverse(Matrix3x3 A, Matrix3x3 Ainverse);

void MatrixProduct(Matrix3x3 A, Matrix3x3 B, Matrix3x3 C);


#endif