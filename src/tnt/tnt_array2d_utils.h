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

#ifndef TNT_ARRAY2D_UTILS_H
#define TNT_ARRAY2D_UTILS_H

// local includes
#include "tnt_array2d.h"

namespace TNT
{

  /**
 * @brief Overload the ouput operator.
 * Print the array on the ostream
 * 
 * @tparam T Type of data
 * @param s outstream
 * @param A input Array to print (this)
 * @return std::ostream& 
 */
  template <class T>
  std::ostream &operator<<(std::ostream &s, const Array2D<T> &A)
  {
    int M = A.dim1();
    int N = A.dim2();

    s << "Dimensions: " << M << " " << N << "\n";
    s << "Components: \n";
    for (int i = 0; i < M; i++)
    {
      for (int j = 0; j < N; j++)
      {
        s << A[i][j] << " ";
      }
      s << "\n";
    }
    s << "\n";

    return s;
  }

  /**
 * @brief Overload the extraction operator. This beauty allows 
 to initialize the array2d from an istream. The first two values must be 
 integers an define the size of the array. The following values 
 initialize the array, if not value is given, value is set to zero. 
 * 
 * @tparam T Type of data
 * @param s input istream
 * @param A Array2D to be created from the istream
 * @return std::istream& input istream
 *
 * Example:
 * @code{.cpp}
   TNT::Array2D<double> A;
   std::istringstream is("3 3 0 1 2 1 2 3 2 3 4");
   is >> A;
 * @endcode
 * 
 */
  template <class T>
  std::istream &operator>>(std::istream &s, Array2D<T> &A)
  {

    int M, N;

    s >> M >> N;

    Array2D<T> B(M, N);

    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
      {
        s >> B[i][j];
      }

    A = B;
    return s;
  }

  /**
 * @brief  Adds two Array2D A and B and creates a *new* Array2D C.
 * 
 * @tparam T 
 * @param A First array to be added
 * @param B Second array to be added
 * @return Array2D<T> C = A + B , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A + B;
 * @endcode
 */
  template <class T>
  Array2D<T> operator+(const Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() != m || B.dim2() != n)
      return Array2D<T>();

    else
    {
      Array2D<T> C(m, n);

      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          C[i][j] = A[i][j] + B[i][j];
      }
      return C;
    }
  }

  /**
 * @brief  Substracts two Array2D A and B and creates a new Array2D C.
 * 
 * @tparam T 
 * @param A First array 
 * @param B Second array 
 * @return Array2D<T> C = A - B , a new array.
  *
 * Example:
 * @code{.cpp}
 * C = A - B;
 * @endcode
 */
  template <class T>
  Array2D<T> operator-(const Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() != m || B.dim2() != n)
      return Array2D<T>();

    else
    {
      Array2D<T> C(m, n);

      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          C[i][j] = A[i][j] - B[i][j];
      }
      return C;
    }
  }

  /**
 * @brief  Multiplies two Array2D A and B, element by element 
 * and creates a new Array2D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param B Second array 
 * @return Array2D<T> C = A * B , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A * B;
 * @endcode
 */
  template <class T>
  Array2D<T> operator*(const Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() != m || B.dim2() != n)
      return Array2D<T>();

    else
    {
      Array2D<T> C(m, n);

      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          C[i][j] = A[i][j] * B[i][j];
      }
      return C;
    }
  }

  /**
 * @brief  Multiplies an Array2D A by a scalar 
 * and creates a new Array2D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param b Scalar value 
 * @return Array2D<T> C = A * b , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A * b;
 * @endcode
 */
  template <class T>
  Array2D<T> operator*(const Array2D<T> &A, const double &b)
  {
    int m = A.dim1();
    int n = A.dim2();

    Array2D<T> C(m, n);

    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
        C[i][j] = A[i][j] * b;
    }
    return C;
  }

  /**
 * @brief  Multiplies a scalar by an Array2D A 
 * and creates a new Array2D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param b Scalar value 
 * @return Array2D<T> C = b * A , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = b * A;
 * @endcode
 */
  template <class T>
  inline Array2D<T> operator*(const T &b, Array2D<T> &A)
  {
    return A * b;
  }

  /**
 * @brief  Divides two Array2D A and B, element by element 
 * and creates a new Array2D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param B Second array 
 * @return Array2D<T> C = A / B , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A / B;
 * @endcode
 */
  template <class T>
  Array2D<T> operator/(const Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() != m || B.dim2() != n)
      return Array2D<T>();

    else
    {
      Array2D<T> C(m, n);

      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          C[i][j] = A[i][j] / B[i][j];
      }
      return C;
    }
  }

  /**
 * @brief  Divides an Array2D A by a scalar b, element by element 
 * and creates a new Array2D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param b Scalar 
 * @return Array1D<T> C = A / b , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A / b;
 * @endcode
 */
  template <class T>
  Array2D<T> operator/(const Array2D<T> &A, const T &b)
  {
    int m = A.dim1();
    int n = A.dim2();

    Array2D<T> C(m, n);

    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
        C[i][j] = A[i][j] / b;
    }
    return C;
  }

  /**
  * @brief  Adds  A and B, element by element 
  * and save the result into A. Avoid to construct a new Array C
  * 
  * @tparam T Type of data
  * @param A First array 
  * @param B Second array 
  * @return Array2D<T> A = A + B. 
  *
  * Example:
  * @code{.cpp}
  * A += B;
  * @endcode
  */
  template <class T>
  Array2D<T> &operator+=(Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() == m || B.dim2() == n)
    {
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          A[i][j] += B[i][j];
      }
    }
    return A;
  }

  /**
  * @brief  Adds  A and B, element by element 
  * and save the result into A. Avoid to construct a new Array C
  * 
  * @tparam T Type of data
  * @param A First array 
  * @param B Second array 
  * @return Array2D<T> A = A - B. 
  *
  * Example:
  * @code{.cpp}
  * A -= B;
  * @endcode
  */
  template <class T>
  Array2D<T> &operator-=(Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() == m || B.dim2() == n)
    {
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          A[i][j] -= B[i][j];
      }
    }
    return A;
  }

  /**
  * @brief  Multiplies A and B, element by element 
  * and save the result into A. Avoid to construct a new Array C
  * 
  * @tparam T Type of data
  * @param A First array 
  * @param B Second array 
  * @return Array2D<T> A = A * B. 
  *
  * Example:
  * @code{.cpp}
  * A *= B;
  * @endcode
  */
  template <class T>
  Array2D<T> &operator*=(Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() == m || B.dim2() == n)
    {
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          A[i][j] *= B[i][j];
      }
    }
    return A;
  }

  /**
   * @brief  Multiplies A by a sacalar, element by element 
  * and save the result into A. Avoid to construct a new Array C
  * 
  * @tparam T Type of data
  * @param A Array1D 
  * @param b Scalar 
  * @return Array1D<T> A = A * b 
  *
  * Example:
  * @code{.cpp}
  * A *= b;
  * @endcode
  */
  template <class T>
  Array2D<T> &operator*=(Array2D<T> &A, const T &b)
  {
    int m = A.dim1();
    int n = A.dim2();

    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
        A[i][j] *= b;
    }
    return A;
  }

  /**
  * @brief  Divides A and B, element by element 
  * and save the result into A. Avoid to construct a new Array C
  * 
  * @tparam T Type of data
  * @param A First array 
  * @param B Second array 
  * @return Array2D<T> A = A / B. 
  *
  * Example:
  * @code{.cpp}
  * A /= B;
  * @endcode
  */
  template <class T>
  Array2D<T> &operator/=(Array2D<T> &A, const Array2D<T> &B)
  {
    int m = A.dim1();
    int n = A.dim2();

    if (B.dim1() == m || B.dim2() == n)
    {
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          A[i][j] /= B[i][j];
      }
    }
    return A;
  }

  /**
  * @brief  Divides an Array1D A by a sacalar b (element by element)
  * and save the result into A. Avoid to construct a new Array C
  * 
  * @tparam T Type of data
  * @param A Array1D 
  * @param b Scalar 
  * @return Array2D<T> A = A/b 
  *
  * Example:
  * @code{.cpp}
  * A /= b;
  * @endcode
  */
  template <class T>
  Array2D<T> &operator/=(Array2D<T> &A, const T &b)
  {
    int m = A.dim1();
    int n = A.dim2();

    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
        A[i][j] /= b;
    }
    return A;
  }

  /* ........................ extended functions ......................*/
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
      if(M[0][0] == 0)
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

      det_m = C11 * (C22 * C33 - C23 * C32)
            - C12 * (C21 * C33 - C23 * C31)
            + C13 * (C21 * C32 - C22 * C31);

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

      det_m = C11 * (C22 * C33 - C23 * C32) 
            - C12 * (C21 * C33 - C23 * C31) 
            + C13 * (C21 * C32 - C22 * C31);
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

} // namespace TNT

#endif
