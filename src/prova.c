#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int gauss(int n, double* A, double* b, double tol) {
  int maxRow;
  // Copy the matrix A to A_aux
  double* A_aux = (double*)malloc(n * n * sizeof(double));
  memcpy(A_aux, A, n * n * sizeof(double));
  double temp, maxElem;
// Swap the ith row with the row that has the maximum element in the ith column
#define A(i, j) A[(i)*n + (j)]
  for (int i = 0; i < n; i++) {
    maxRow = i;
    maxElem = fabs(A(i, i));
    for (int j = i + 1; j < n; j++) {
      if (fabs(A(j, i)) > maxElem) {
        maxRow = j;
        maxElem = fabs(A(j, i));
      }
    }

    // If the maximum element in the ith column is too small, the matrix is singular
    if (maxElem < tol) return 1;
    // Swap the ith row with the row that has the maximum element in the ith column
    if (maxRow != i) {
      for (int j = 0; j < n; j++) {
        temp = A(i, j);
        A(i, j) = A(maxRow, j);
        A(maxRow, j) = temp;
      }
      temp = b[i];
      b[i] = b[maxRow];
      b[maxRow] = temp;
    }

    // Eliminate the ith column of the matrix
    for (int j = i + 1; j < n; j++) {
      temp = A(j, i) / A(i, i);
      A(j, i) = 0;
      for (int k = i + 1; k < n; k++)
        A(j, k) -= temp * A(i, k);
      b[j] -= temp * b[i];
    }
  }
  // From here the matrix is upper triangular
  // Back-substitution to find the solution x (stored in b)
  for (int i = n - 1; i >= 0; i--) {
    for (int j = i + 1; j < n; j++)
      b[i] -= A(i, j) * b[j];
    b[i] /= A(i, i);
  }
#undef A

  // Copy the matrix A_aux back to A
  memcpy(A, A_aux, n * n * sizeof(double));
  free(A_aux);
  return 0;
}

int main(int argc, char const* argv[]) {
  const int n = 4;
  double tol = 0.000001;
  double* A = (double*)malloc(n * n * sizeof(double*));
  A[0] = 10;
  A[1] = -1;
  A[2] = 2;
  A[3] = 0;
  A[4] = -1;
  A[5] = 11;
  A[6] = -1;
  A[7] = 3;
  A[8] = 2;
  A[9] = -1;
  A[10] = 10;
  A[11] = -1;
  A[12] = 0;
  A[13] = 3;
  A[14] = -1;
  A[15] = 8;

  // A[0][0] = 10;
  // A[0][1] = -1;
  // A[0][2] = 2;
  // A[0][3] = 0;
  // A[1][0] = -1;
  // A[1][1] = 11;
  // A[1][2] = -1;
  // A[1][3] = 3;
  // A[2][0] = 2;
  // A[2][1] = -1;
  // A[2][2] = 10;
  // A[2][3] = -1;
  // A[3][0] = 0;
  // A[3][1] = 3;
  // A[3][2] = -1;
  // A[3][3] = 8;
  // A = {
  //     {0, 3, -1, 8},
  //     {-1, 11, -1, 3},
  //     {2, -1, 10, -1},
  //     {10, -1, 2, 0}};
  double* b = (double*)malloc(n * sizeof(double));
  b[0] = 6;
  b[1] = 25;
  b[2] = -11;
  b[3] = 15;
  // b= {15, 25, -11, 6};

  if (gauss(n, A, b, tol)) {
    printf("Maximum number of iterations exceeded.\n");
    return 1;
  }

  for (int i = 0; i < n; i++) {
    printf("A[%d] = [", i);
    for (int j = 0; j < n; j++)
      printf("%lf ", A[i * n + j]);
    printf("]\n");
  }

  printf("x = [");
  for (int i = 0; i < n; i++) {
    printf("%lf ", b[i]);
  }
  printf("]\n");

  return 0;
}
