/*
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/

#ifndef TNT_ARRAY2D_LINALG_H
#define TNT_ARRAY2D_LINALG_H

// local includes
#include "tnt_array2d.h"

namespace TNT
{
  /**
 * @brief Returns an Array1D containing the sum of each column of the Array2D
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T Array1D containing the sum of each column of the Array2D
 * 
 * Example:
 * @code{.cpp}
 * double a_values[3][3] = {{2, 7, 6},
                           {9, 5, 1},
                           {4, 3, 8}};
 * TNT::Array2D<double> A(3, 3, *a_values);
 *
 * TNT::Array1D<double> sum_col = TNT::sum(A); // 15 x 3
 * @endcode
 */
  template <class T>
  Array1D<T> sum(const Array2D<T> &A)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (n < 1 || m < 1)
      return Array1D<T>();

    Array1D<T> S(n, 0.);
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < m; i++)
        S[j] += A[i][j];
    }

    return S;
  }

  /**
 * @brief Returns the 1-norm of matrix A. 
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T  the maximum absolute column sum of the matrix.
 * 
 * Example:
 * @code{.cpp}
  double a_values[3][3] = {{5, -4, 2},
                           {-1, 2, 3},
                           {-2, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  double norm_A = TNT::norm_1(A); // 8.0
 * @endcode
 */
  template <class T>
  T norm_1(const Array2D<T> &A)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (n < 1 || m < 1)
      return 0.0;

    double f = 0;
    for (int j = 0; j < n; j++)
    {
      double s = 0;
      for (int i = 0; i < m; i++)
      {
        s += std::abs(A[i][j]);
      }
      f = std::max(f, s);
    }
    return f;
  }

  /**
 * @brief Returns the inf-norm of the matrix A.
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T the maximum absolute row sum of the matrix A.
 * 
 * Example:
 * @code{.cpp}
  double a_values[3][3] = {{5, -4, 2},
                           {-1, 2, 3},
                           {-2, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  double norm_A = TNT::norm_inf(A); // 11.0
 * @endcode
 */
  template <class T>
  T norm_inf(const Array2D<T> &A)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (n < 1 || m < 1)
      return 0.0;

    double f = 0;
    for (int i = 0; i < m; i++)
    {
      double s = 0;
      for (int j = 0; j < n; j++)
      {
        s += std::abs(A[i][j]);
      }
      f = std::max(f, s);
    }
    return f;
  }

  /**
 * @brief Returns the frobenius norm of the matrix A.
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T the maximum absolute row sum of the matrix A.
 * 
 * Example:
 * @code{.cpp}
  double a_values[3][3] = {{5, -4, 2},
                           {-1, 2, 3},
                           {-2, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  double norm_A = TNT::norm_fro(A); // 80.0
 * @endcode
 */
  template <class T>
  T norm_fro(const Array2D<T> &A)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (n < 1 || m < 1)
      return 0.0;

    T sum = 0.;
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < m; i++)
      {
        sum += A[i][j] * A[i][j];
      }
    }
    return std::sqrt(sum);
  }

  /**
 * @brief Returns an Array1D containing the sum in a given direction 
 * of the Array2D
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @param i direction 1: rows, 2: columns.
 * @return T Sum: the sum of the elements in the direction 1
 * 
 * Example:
 * @code{.cpp}
 * double a_values[3][3] = {{2, 7, 6},
                           {9, 5, 1},
                           {4, 3, 8}};
 * TNT::Array2D<double> A(3, 3, *a_values);
 *
 * TNT::Array1D<double> sum_col = TNT::sum(A,2); // 15 x 3
 * TNT::Array1D<double> sum_row = TNT::sum(A,1); // 15 x 3
 * @endcode
 */
  template <class T>
  Array1D<T> sum(const Array2D<T> &A, int i)
  {
    switch (i)
    {
    case 1:
      return sum(transpose(A));
      break;
    case 2:
      return sum(A);
      break;
    default:
      return Array1D<T>();
      break;
    }
  }


  /**
 * @brief Returns the trace of an Array1D
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T Array1D containing the sum of diagonal fo the array
 * 
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(); // + Initialization 
 * double res = TNT::trace(A);
 * @endcode
 */
  template <class T>
  T trace(const Array2D<T> &A)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (n < 1 || m < 1)
      return 0.;

    T sum = 0.;

    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (i == j)
          sum += A[i][j];
      }
    }

    return sum;
  }

  /**
  * @brief Matrix Multiply:  compute C = A*B, where C[i][j] is the dot-product 
  * of row i of A and column j of B.
  *
  * @param A an (m x n) array
  * @param B an (n x k) array
  * @return the (m x k) array A*B, or a null array (0x0) if the matrices are
  * non-conformant (i.e. the number of columns of A are 
  * different than the number of rows of B.)
    * Example:
  * @code{.cpp}
  * TNT::Array2D<double>
  * C = TNT::matmult(A, B);
  * @endcode
  */
  template <class T>
  Array2D<T> matmult(const Array2D<T> &A, const Array2D<T> &B)
  {
    if (A.dim2() != B.dim1())
      return Array2D<T>();

    int M = A.dim1();
    int N = A.dim2();
    int K = B.dim2();

    Array2D<T> C(M, K);

    T sum = 0;

    for (int i = 0; i < M; i++)
      for (int j = 0; j < K; j++)
      {
        for (int k = 0; k < N; k++)
        {
          sum += A[i][k] * B[k][j];
        }
        C[i][j] = sum;
        sum = 0;
      }

    return C;
  }

  /**
 * @brief Transpose: return a new array2d B, where B(i,j) is A(j,i).
 * 
 * @tparam T data type
 * @param A matrix MxN
 * @return TNT::Array2D<T> of size N x M, where each (i,j) is  A(j,i).
 *
 * Example:
 * @code{.cpp}
 * double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
   TNT::Array2D<double> A(3, 3, *a_values);
   TNT::Array2D<double> B = TNT::transpose(A);
 * @endcode
 */
  template <class T>
  TNT::Array2D<T> transpose(const TNT::Array2D<T> &A)
  {
    int M = A.dim1();
    int N = A.dim2();

    TNT::Array2D<T> B(N, M);

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        B[j][i] = A[i][j];

    return B;
  }

  /**
 * @brief Matrix-Matrix tranpose multiplication, i.e. compute tranpose(A)*B.
 * 
 * @tparam T data type
 * @param A matrix: size M x N.
 * @param B matrix: size M x K.
 * @return Array2D<T> a new matrix of size N x K.
 *
 * @note this is more efficient than computing the tranpose(A) explicitly,
 * and then multiplying, as the tranpose of A is never really constructed.
 */
  template <class T>
  Array2D<T> transpose_mult(const Array2D<T> &A, const Array2D<T> &B)
  {

    if (A.dim1() != B.dim1())
      return Array2D<T>();

    int M = A.dim1();
    int N = A.dim2();
    int K = B.dim2();

    Array2D<T> C(M, K);
    T sum;

    for (int i = 0; i < N; i++)
      for (int k = 0; k < K; k++)
      {
        sum = 0;
        for (int j = 0; j < M; j++)
          sum = sum + A[j][i] * B[j][k];

        C[i][k] = sum;
      }

    return C;
  }

  /**
 * @brief Tranpose Matrix-Matrix, i.e. compute transpose(A)*B.

 * 
 * @tparam T 
 * @param A matrix: size M x N.
 * @param B matrix: size K x N.
 * @return Array2D<T> a new matrix of size N x K.
 */
  template <class T>
  Array2D<T> mult_transpose(const Array2D<T> &A, const Array2D<T> &B)
  {

    if (A.dim2() != B.dim2())
      return Array2D<T>();

    int M = A.dim1();
    int N = A.dim2();
    int K = B.dim1();

    Array2D<T> C(M, K);
    T sum;

    for (int i = 0; i < M; i++)
      for (int k = 0; k < K; k++)
      {
        sum = 0;
        for (int j = 0; j < N; j++)
          sum = sum + A[i][j] * B[k][j];

        C[i][k] = sum;
      }

    return C;
  }

  /**
 * @brief Matrix-Vector multiplication, i.e. compute A*b.
 * 
 * @tparam T 
 * @param A Matrix: size M x N.
 * @param b Vector: size M.
 * @return Array1D<T> a new vector of size N.
 */
  template <class T>
  Array1D<T> matmult(const Array2D<T> &A, const Array1D<T> &b)
  {

    if (A.dim2() != b.dim1())
      return Array1D<T>();

    int M = A.dim1();
    int N = A.dim2();

    Array1D<T> tmp(M);

    for (int i = 0; i < M; i++)
    {
      T sum = 0;
      for (int j = 0; j < N; j++)
        sum = sum + A[i][j] * b[j];

      tmp[i] = sum;
    }

    return tmp;
  }

  /**
 * @brief Matrix-Vector transpose multiplication, i.e. compute tranpose(A)*b.
 * 
 * @tparam T 
 * @param A Matrix: size M x N.
 * @param b Vector: size M.
 * @return Array1D<T> a new vector of size N.
 */
  template <class T>
  Array1D<T> transpose_mult(const Array2D<T> &A, const Array1D<T> &b)
  {

    if (A.dim1() != b.dim1())
      return Array1D<T>();

    int M = A.dim2();
    int N = A.dim1();

    Array1D<T> tmp(M);

    for (int i = 0; i < M; i++)
    {
      T sum = 0;
      for (int j = 0; j < N; j++)
        sum = sum + A[j][i] * b[j];

      tmp[i] = sum;
    }

    return tmp;
  }

  /**
 * @brief Returns the invert of a square matrix. if the dimension is <3
 * the inverse is given explicitlely, if not the LU method is used.
 * 
 * @tparam T data type
 * @param M Input MxM matrix
 * @return TNT::Array2D<T> Inverse MxM Matrix 
 */
  template <class T>
  TNT::Array2D<T> invert(const TNT::Array2D<T> &M)
  {
    // square matrices only
    assert(M.dim1() == M.dim2());
    int m = M.dim1();

    T det_m;
    TNT::Array2D<T> inv(m, m);

    T C11, C12, C13;
    T C21, C22, C23;
    T C31, C32, C33;

    switch (m)
    {
    case 1:
      if (M[0][0] == 0)
        return Array2D<T>();
      inv[0][0] = 1. / M[0][0];
      break;
    case 2:
      det_m = M[0][0] * M[1][1] - M[0][1] * M[1][0];
      if (det_m == 0)
        return Array2D<T>();
      inv[0][0] = M[1][1] / det_m;
      inv[1][1] = M[0][0] / det_m;
      inv[0][1] = -M[0][1] / det_m;
      inv[1][0] = -M[1][0] / det_m;
      break;
    case 3:
      C11 = M[0][0];
      C12 = M[0][1];
      C13 = M[0][2];
      C21 = M[1][0];
      C22 = M[1][1];
      C23 = M[1][2];
      C31 = M[2][0];
      C32 = M[2][1];
      C33 = M[2][2];

      det_m = C11 * (C22 * C33 - C23 * C32) - C12 * (C21 * C33 - C23 * C31) + C13 * (C21 * C32 - C22 * C31);

      if (det_m == 0)
        return Array2D<T>();

      inv[0][0] = (C22 * C33 - C32 * C23) / det_m;
      inv[1][0] = (C23 * C31 - C33 * C21) / det_m;
      inv[2][0] = (C21 * C32 - C31 * C22) / det_m;
      inv[0][1] = (C32 * C13 - C12 * C33) / det_m;
      inv[1][1] = (C33 * C11 - C13 * C31) / det_m;
      inv[2][1] = (C31 * C12 - C11 * C32) / det_m;
      inv[0][2] = (C12 * C23 - C22 * C13) / det_m;
      inv[1][2] = (C13 * C21 - C23 * C11) / det_m;
      inv[2][2] = (C11 * C22 - C21 * C12) / det_m;

      break;
    default:
      /* @todo
      // solve for inverse with LU decomposition
      JAMA::LU<T> lu(M);
      // create identity matrix
      TNT::Array2D<T> id(M.dim1(), M.dim2(), (T)0);
      for (int i = 0; i < M.dim1(); i++)
        id[i][i] = 1;
      // solves A * A_inv = Identity
      inv = lu.solve(id);
      
      */
      break;
    }
    return inv;
  }

  /**
 * @brief Returns the determinant of a square matrix. if the dimension is <3
 * the inverse is given explicitlely, if not the LU method is used.
 * 
 * @tparam T data type
 * @param M Input MxM matrix
 * @return T the determinant of the matrix
 */
  template <class T>
  T det(const TNT::Array2D<T> &M)
  {
    // square matrices only
    assert(M.dim1() == M.dim2());
    int m = M.dim1();

    T det_m;

    T C11, C12, C13;
    T C21, C22, C23;
    T C31, C32, C33;

    switch (m)
    {
    case 1:
      det_m = M[0][0];
      break;
    case 2:
      det_m = M[0][0] * M[1][1] - M[0][1] * M[1][0];
      break;
    case 3:
      C11 = M[0][0];
      C12 = M[0][1];
      C13 = M[0][2];
      C21 = M[1][0];
      C22 = M[1][1];
      C23 = M[1][2];
      C31 = M[2][0];
      C32 = M[2][1];
      C33 = M[2][2];

      det_m = C11 * (C22 * C33 - C23 * C32) - C12 * (C21 * C33 - C23 * C31) + C13 * (C21 * C32 - C22 * C31);
      break;
    default:
      /* @todo
      // solve for inverse with LU decomposition
      JAMA::LU<T> lu(M);
      // create identity matrix
      TNT::Array2D<T> id(M.dim1(), M.dim2(), (T)0);
      for (int i = 0; i < M.dim1(); i++)
        id[i][i] = 1;
      // solves A * A_inv = Identity
      inv = lu.solve(id);
      
      */
      break;
    }
    return det_m;
  }

  /**
   * @brief Verify if the all the values of A are equal to the values of B
   * 
   * @tparam T data type
   * @param A First array to compare
   * @param B Second Array to compare
   * @return true if A==B
   * @return false if at least one value is different 
   */
  template <class T>
  bool operator==(const Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() != m || B.dim2() != n)
      return false;

    else
    {
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          if (A[i][j] != B[i][j])
            return false;
      }
      return true;
    }
  }

  /**
   * @brief Verify if the all the values of A are near to the values of B
   * 
   * @tparam T data type
   * @param A First array to compare
   * @param B Second Array to compare
   * @return true if A==B
   * @return false if at least one value is different 
   */
  template <class T>
  bool near(const Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    T eps = std::numeric_limits<T>::epsilon();

    if (B.dim1() != m || B.dim2() != n)
      return false;

    T norm_A = norm_1(A);
    T norm_B = norm_1(B);

    if (norm_A == 0. && norm_B < 10. * eps)
      return false;
    if (norm_B == 0. && norm_A < 10. * eps)
      return false;

    if (TNT::norm_1(A - B) > 100. * eps * std::max(norm_A, norm_B))
      return false;
    else
      return true;
  }

  /**
   * @brief Verify if the all the values of A are near to the values of, value per value 
   * given a tolerance
   * 
   * @tparam T data type
   * @param A First array to compare
   * @param B Second Array to compare
   * @param tol tolerance for the comparaison
   * @return true if abs(A-B) > tol 
   * @return false if at least one value is different 
   */
  template <class T>
  bool near(const Array2D<T> &A, const Array2D<T> &B, T tol)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() != m || B.dim2() != n)
      return false;

    for (int j = 0; j < n; j++)
      for (int i = 0; i < m; i++)
        if (std::abs(A[i][j]-B[i][j]) > tol)
          return false;

    return true;
  }

} // namespace TNT

#endif
