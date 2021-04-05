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

  

} // namespace TNT

#endif
