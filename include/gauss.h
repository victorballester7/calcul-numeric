#ifndef GAUSS_H
#define GAUSS_H

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

// ----------------------------------------------
// // Gauss-Seidel method
// // ----------------------------------------------
// // Purpose:
// //  Solve a system of linear equations Ax = b using the Gauss-Seidel method
// //
// // Arguments:
// //  n: number of equations
// //  A: matrix of coefficients
// //  b: vector of constants
// //  x: initial guess and is where the solution is stored
// //  tol: tolerance
// //  maxIter: maximum number of iterations
// //
// // Returns:
// //  0 if converged, 1 if not converged
// // ----------------------------------------------
// int GaussSeidel(int n, double** A, double* b, double* x, double tol, int maxIter);

// ----------------------------------------------
// Gauss elimination method
// ----------------------------------------------
// Purpose:
//  Solve a system of linear equations Ax = b using the Gauss elimination method
//
// Arguments:
//  n: number of equations
//  A: matrix of coefficients (as a vector)
//  b: vector of constants and is where the solution is stored at the end
//  tol: tolerance for the pivot to be considered zero
//
// Returns:
//  0 if the system has a unique solution, 1 if the matrix is singular
int gauss(int n, double *A, double *b, double tol);
int resgauss(int n, double *a, double *b, double *x);
#endif  // GAUSS_H
