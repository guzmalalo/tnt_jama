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

#ifndef TNT_ARRAY1D_LINALG_H
#define TNT_ARRAY1D_LINALG_H

// local includes
#include "tnt_array1d.h"

namespace TNT
{

  /**
 * @brief Finds the max element in an array1D A.
 * 
 * @tparam T Data type
 * @param A an input Array1D
 * @return T max  = \f$ max(A) \f$
 *
 * Example:
 * @code{.cpp}
 * double a_values[6] = {-2, 0, 3, 4,+100,1e5};
 * TNT::Array1D<double> A(6, a_values);
 * double C = TNT::max(A);
 * @endcode
 */
  template <class T>
  T max(const Array1D<T> &A)
  {
    int n = A.dim1();

    if (n < 1)
      return 0.;

    T max_e = A[0];
    for (int i = 0; i < n; i++)
    {
      max_e = std::max(max_e, A[i]);
    }
    return max_e;
  }

  /**
 * @brief Finds the min element in an array1D A.
 * 
 * @tparam T Data type
 * @param A an input Array1D
 * @return T max  = \f$ min(A) \f$
 *
 * Example:
 * @code{.cpp}
 * double a_values[6] = {-2, 0, 3, 4,+100,1e5};
 * TNT::Array1D<double> A(6, a_values);
 * double C = TNT::min(A);
 * @endcode
 */
  template <class T>
  T min(const Array1D<T> &A)
  {
    int n = A.dim1();

    if (n < 1)
      return 0.;

    T min_e = A[0];
    for (int i = 0; i < n; i++)
    {
      min_e = std::min(min_e, A[i]);
    }
    return min_e;
  }

  /**
 * @brief  Gives the absolute absolue value of an Array1D 
 * in a new Array1D C.
 * 
 * @tparam T Type of data
 * @param A Input array 
 * @return Array1D<T> C[i] = |A[i]| , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = TNT::abs(A);
 * @endcode
 */
  template <class T>
  Array1D<T> abs(const Array1D<T> &A)
  {
    int n = A.dim1();
    Array1D<T> C(n);

    for (int i = 0; i < n; i++)
    {
      C[i] = std::abs(A[i]);
    }
    return C;
  }

   /**
 * @brief Sum of array elements
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T Sum: the sum of the elements. 
 * 
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(100,1.); 
 * double res = TNT::sum(A); // Expected equal to 100
 * @endcode
 */
  template <class T>
  T sum(const Array1D<T> &A)
  {
    int n = A.dim1();

    T sum = 0.;
    for (int i = 0; i < n; i++)
      sum += A[i] ;
    return sum;
  }

  /**
 * @brief Computes the infinity norm of an Array1d.
 * 
 *  \f$ \left\| \mathbf{x} \right\| _\infty := 
 *  \max \left( \left| x_1 \right| , \ldots , \left| x_n \right| \right) \f$
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T the infinity norm 
 * 
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(); // + Initialization 
 * double res = TNT::norm_inf(A);
 * @endcode
 */
  template <class T>
  T norm_inf(const Array1D<T> &A)
  {
    int n = A.dim1();

    if (n < 1)
      return 0.;

    T max_e = std::abs(A[0]);
    for (int i = 0; i < n; i++)
    {
      max_e = std::max(max_e, std::abs(A[i]));
    }
    return max_e;
  }

  /**
 * @brief Computes the l2​-norm of an Array1d (Euclidean norm)
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T the norm of the array 
 * 
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(); // + Initialization 
 * double res = TNT::norm(A);
 * @endcode
 */
  template <class T>
  T norm(const Array1D<T> &A)
  {
    int n = A.dim1();

    T sum = 0.;
    for (int i = 0; i < n; i++)
      sum += A[i] * A[i];
    return std::sqrt(sum);
  }

  /**
 * @brief Computes the l1-norm of an Array1d
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @return T the norm of the array 
 * 
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(); // + Initialization 
 * double res = TNT::norm_2(A);
 * @endcode
 */
  template <class T>
  T norm_1(const Array1D<T> &A)
  {
    int n = A.dim1();

    T sum = 0.;
    for (int i = 0; i < n; i++)
      sum += std::abs(A[i]);
    return sum;
  }

  /**
 * @brief Computes the generalized p-norm of an Array1d
 * 
 * @tparam T Type od data
 * @param A Input array 
 * @param p norm dimension
 * @return T the norm of the array 
 * 
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(); // + Initialization 
 * double res = TNT::norm_2(A);
 * @endcode
 */
  template <class T>
  T norm_p(const Array1D<T> &A, int p)
  {
    int n = A.dim1();

    T sum = 0.;
    for (int i = 0; i < n; i++)
      sum += std::pow(A[i],p);
    return std::pow(sum,1.0/p);
  }

  
  /**
 * @brief  Computes the cross product of A and B, this
 * operation is  only defined in a three-dimensional space.
 * 
 * @param A an Array1D
 * @param B an Array1D
 * @return Array1D<T> C  = \f$ A \times B \f$
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(); // + Initialization 
 * TNT::Array1D<double> B(); // + Initialization 
 * double res = TNT::cross(A,B);
 * @endcode
 */
  template <class T>
  Array1D<T> cross(const Array1D<T> &A, const Array1D<T> &B)
  {
    int n = A.dim1();
    int m = B.dim1();

    if (n != 3 || m != 3)
      return Array1D<T>();
    else
    {
      Array1D<T> C(n);
      C[0] = A[1] * B[2] - A[2] * B[1];
      C[1] = A[2] * B[0] - A[0] * B[2];
      C[2] = A[0] * B[1] - A[1] * B[0];
      return C;
    }
  }

  /**
 * @brief Computes the dot product of two Array1D A and B.
 * 
 * @tparam T Data type
 * @param A an input Array1D
 * @param B an input Array1D
 * @return T C  = \f$ A \cdot B \f$
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(); // + Initialization 
 * TNT::Array1D<double> B(); // + Initialization 
 * double C = TNT::dot_product(A, B);
 * @endcode
 */
  template <class T>
  T dot_product(const Array1D<T> &A, const Array1D<T> &B)
  {
    int n = A.dim1();
    int m = B.dim1();

    if (n != 3 || m != 3)
      return 0;
    else
    {
      T sum = 0;
      for (int i = 0; i < n; i++)
        sum += A[i] * B[i];

      return sum;
    }
  }



} // namespace TNT

#endif
