#include "../include/gauss.h"

// int GaussSeidel(int n, double** A, double* b, double* x, double tol, int maxIter) {
//   // rearrange the rows of A to ensure that A(i,i) != 0
//   double *tmpA, tmpb;
//   int index[n];
//   for (int i = 0; i < n; i++) {
//     if (A[i][i] == 0) {
//       for (int j = i + 1; j < n; j++) {
//         if (A(j,i) != 0) {
//           tmpA = A[i];
//           A[i] = A[j];
//           A[j] = tmpA;
//           tmpb = b[i];
//           b[i] = b[j];
//           b[j] = tmpb;
//           index[i] = j;
//           break;
//         }
//       }
//     } else
//       index[i] = i;
//   }
//   int iter = 0, sum;
//   double* x0 = (double*)malloc(n * sizeof(double));
//   memcpy(x0, x, n * sizeof(double));
//   double error = tol + 1;
//   // start the iteration
//   while (error > tol && iter < maxIter) {
//     for (int i = 0; i < n; i++) {
//       sum = 0;
//       for (int j = 0; j < n; j++) {
//         if (j != i) sum += A(i,j) * x[j];
//       }
//       x[i] = (b[i] - sum) / A[i][i];
//     }
//     error = 0;
//     for (int i = 0; i < n; i++) {
//       error += fabs(x[i] - x0[i]);
//       x0[i] = x[i];
//     }
//     iter++;
//   }
//   free(x0);
//   if (iter == maxIter) {
//     // printf("Maximum number of iterations exceeded.\n");
//     return 1;
//   }
//   // rearrange the rows of A, b and x to their original order
//   double aux_A[n][n], aux_b[n];
//   for (int i = 0; i < n; i++) {
//     memcpy(aux_A[index[i]], A[i], n * sizeof(double));
//     aux_b[index[i]] = b[i];
//     x0[index[i]] = x[i];
//   }
//   memcpy(A, aux_A, n * n * sizeof(double));
//   memcpy(b, aux_b, n * sizeof(double));
//   memcpy(x, x0, n * sizeof(double));

//   return 0;
// }

// Gaussian elimination
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

int resgauss(int n, double* a, double* b, double* x) {
#define A(j, i) a[((j)-1) * n + ((i)-1)]

  int i, j, k;
  float pivot = 0.0;
  float factor = 0.0;
  float sum = 0.0;

  //    for(i=1;i<=n;i++){
  //        for(j=1;j<=n;j++)
  //            printf("\t%10.4f",A(i,j));
  //        printf("\t%10.4f",b[i-1]);
  //        printf("\n\n");
  //    }

  for (k = 1; k <= n - 1; k++) {
    if (A(k, k) == 0.0) {
      return (1);
    } else {
      pivot = A(k, k);
      for (j = k; j <= n; j++)
        A(k, j) = A(k, j) / pivot;
      b[k - 1] = b[k - 1] / pivot;

      for (i = k + 1; i <= n; i++) {
        factor = A(i, k);

        for (j = k; j <= n; j++)
          A(i, j) = A(i, j) - factor * A(k, j);
        b[i - 1] = b[i - 1] - factor * b[k - 1];
      }
    }

    if (A(n, n) == 0)
      return (1);

    else {
      x[n - 1] = b[n - 1] / A(n, n);

      for (i = n - 1; i >= 1; i--) {
        sum = 0.0;
        for (j = i + 1; j <= n; j++)
          sum = sum + A(i, j) * x[j - 1];
        x[i - 1] = (b[i - 1] - sum) / A(i, i);
      }
    }
  }

  return (0);

#undef A
}
