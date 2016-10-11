#include "Matrix.hh"

void MatrixInverse(Matrix3x3 A, Matrix3x3 Ainverse)
{
	double a00 = A[0][0];
	double a01 = A[0][1];
	double a02 = A[0][2]; 
    double a10 = A[1][0];
	double a11 = A[1][1];
	double a12 = A[1][2];
	double a20 = A[2][0];
	double a21 = A[2][1];
	double a22 = A[2][2];

	/* Reciprocal of determinant */
    double dinverse = 1.0 / ( 
        a00 * a11 * a22  -  a00 * a12 * a21  -  a01 * a10 * a22 + 
        a01 * a12 * a20 + a02 * a10 * a21  -  a02 * a11 * a20);


    /*  calculate inverse matrix  */
    Ainverse[0][0] = dinverse  *  (a11 * a22 - a12 * a21);
    Ainverse[0][1] = dinverse  *  (a02 * a21 - a01 * a22);
    Ainverse[0][2] = dinverse  *  (a01 * a12 - a02 * a11);
    
    Ainverse[1][0] = dinverse  *  (a12 * a20 - a10 * a22);
    Ainverse[1][1] = dinverse  *  (a00 * a22 - a02 * a20);
    Ainverse[1][2] = dinverse  *  (a02 * a10 - a00 * a12);
    
    Ainverse[2][0] = dinverse  *  (a10 * a21 - a11 * a20);
    Ainverse[2][1] = dinverse  *  (a01 * a20 - a00 * a21);
    Ainverse[2][2] = dinverse  *  (a00 * a11 - a01 * a10);
}

void MatrixProduct(Matrix3x3 A, Matrix3x3 B, Matrix3x3 result)
{
	result[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
    result[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
    result[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

    result[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
    result[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
    result[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

    result[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
    result[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
    result[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
}