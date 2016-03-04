#pragma once

namespace Covariance {

   /* Recursively compute the determinant of the NxN matrix A
    */
   template<typename T> T recurse_determinant(T* A, unsigned int N) {
      if(N == 1) {
         return A[0];
      } else if(N == 2) {
         return A[0]*A[3]-A[1]*A[2];
      } else {
         // Size of a submatrix
         unsigned int NN = (N-1)*(N-1);

         // Calculate the cofactor for the first line
         T det = 0.0f;
         for(unsigned int i=0; i<N; ++i) {
            T B[NN];
            unsigned int U=0;
            for(unsigned int u=0; u<N; ++u) {
               if(u == i)
                  continue;

               for(unsigned int v=1; v<N; ++v) {

                  B[U*(N-1)+v-1] = A[u*N+v];
               }

               U++;
            }

            det += pow(-1, i) * A[i*N] * recurse_determinant(B, N-1);
         }

         return det;
      }
   }

   /* Compute the cofactor for the given NxN matrix A
    */
   template<typename T> T cofactor(T* A, int n, int i, int j) {
      // Create the cofactor
      T C[(n-1)*(n-1)];
      int I = 0;
      int J = 0;
      for(int ii=0; ii<n; ++ii) {
         if(ii == i) continue;

         J = 0;
         for(int jj=0; jj<n; ++jj) {
            if(jj == j) continue;

            C[I*(n-1)+J] = A[ii*n+jj];

            ++J;
         }
         ++I;
      }

      return recurse_determinant(C, n-1);
   }

   /* Calculate the inverse of the matrix A using the cofactor and the inverse
    * determinant. Return 'true' if the determinant can be calculated and put
    * the resulting matrix in place 'A'.
    */
   template<typename T> bool Inverse(T* A, int size) {
      T B[size*size];
      T det = recurse_determinant<T>(A, size);

      if(det < T(0.0)) {
         return false;
      }

      for(int i=0; i<size; ++i)
         for(int j=0; j<size; ++j) {
            B[i*size + j] = pow(-1, i+j+2) * (cofactor<T>(A, size, j, i) / det);
         }

      for(int i=0; i<size; ++i)
         for(int j=0; j<size; ++j) {
            A[i*size +j] = B[i*size + j];
         }

      return true;
   }

   /* Calculate the determinant of the matrix A using the recursive method.
    * Since the covariance is a symmetric positive matrix, the determinant
    * should always be positive.
    */
   template<typename T>
   T Determinant(T* A, int size) {
      return recurse_determinant<T>(A, size);
   }
}
