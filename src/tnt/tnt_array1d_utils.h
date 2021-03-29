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

#ifndef TNT_ARRAY1D_UTILS_H
#define TNT_ARRAY1D_UTILS_H

#include <cstdlib>
#include <cassert>

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
std::ostream& operator<<(std::ostream &s, const Array1D<T> &A)
{
    int N=A.dim1();

#ifdef TNT_DEBUG
	s << "addr: " << (void *) &A[0] << "\n";
#endif
    s << "Dimension: " << N << "\n";
	s << "Components: \n";
    for (int j=0; j<N; j++)
    {
       s << A[j] << "\n";
    }
    s << "\n";

    return s;
}

/**
 * @brief Overload the extraction operator. This beauty allows to initialize the array from an istream. The first value must be 
 an integer an defines the size of the array. The following values 
 initialize the array, if not value is given, value is set to zero. 
 * 
 * @tparam T Type of data
 * @param s input istream
 * @param A Array1D to be created from the istream
 * @return std::istream& input istream
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(a_size_);
 * std::istringstream is("5 1 2 3 4 5");
 * is >> A;
 * @endcode
 * 
 */
template <class T>
std::istream& operator>>(std::istream &s, Array1D<T> &A)
{
	int N;
	s >> N;

	Array1D<T> B(N);
	for (int i=0; i<N; i++)
		s >> B[i];
	A = B;
	return s;
}


/**
 * @brief  Adds two Array1D A and B and creates a new Array1D C.
 * 
 * @tparam T 
 * @param A First array to be added
 * @param B Second array to be added
 * @return Array1D<T> C = A + B , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A + B;
 * @endcode
 */
template <class T>
Array1D<T> operator+(const Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() != n )
		return Array1D<T>();

	else
	{
		Array1D<T> C(n);

		for (int i=0; i<n; i++)
		{
			C[i] = A[i] + B[i];
		}
		return C;
	}
}


/**
 * @brief  Substracts two Array1D A and B and creates a new Array1D C.
 * 
 * @tparam T 
 * @param A First array 
 * @param B Second array 
 * @return Array1D<T> C = A - B , a new array.
  *
 * Example:
 * @code{.cpp}
 * C = A - B;
 * @endcode
 */
template <class T>
Array1D<T> operator-(const Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() != n )
		return Array1D<T>();

	else
	{
		Array1D<T> C(n);

		for (int i=0; i<n; i++)
		{
			C[i] = A[i] - B[i];
		}
		return C;
	}
}


/**
 * @brief  Multiplies two Array1D A and B, element by element 
 * and creates a new Array1D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param B Second array 
 * @return Array1D<T> C = A * B , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A * B;
 * @endcode
 */
template <class T>
Array1D<T> operator*(const Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() != n )
		return Array1D<T>();

	else
	{
		Array1D<T> C(n);

		for (int i=0; i<n; i++)
		{
			C[i] = A[i] * B[i];
		}
		return C;
	}
}

/**
 * @brief  Multiplies an Array1D A by a scalar 
 * and creates a new Array1D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param b Scalar value 
 * @return Array1D<T> C = A * b , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A * b;
 * @endcode
 */
template <class T>
Array1D<T>  operator*(Array1D<T> &A, const T &b)
{
	int n = A.dim1();
    Array1D<T> C(n);

	for (int i=0; i<n; i++)
	{
		C[i] = A[i]*b;
	}
	return C;
}

/**
 * @brief  Multiplies a scalar by an Array1D
 * and creates a new Array1D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param b Scalar value 
 * @return Array1D<T> C = b * A , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = b * A;
 * @endcode
 */
template <class T>
inline Array1D<T>  operator*(const T &b, Array1D<T> &A)
{
	return A*b;
}

/**
 * @brief  Divides two Array1D A and B, element by element 
 * and creates a new Array1D C.
 * 
 * @tparam T Type of data
 * @param A First array 
 * @param B Second array 
 * @return Array1D<T> C = A / B , a new array.
 *
 * Example:
 * @code{.cpp}
 * C = A / B;
 * @endcode
 */
template <class T>
Array1D<T> operator/(const Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() != n )
		return Array1D<T>();

	else
	{
		Array1D<T> C(n);

		for (int i=0; i<n; i++)
		{
			C[i] = A[i] / B[i];
		}
		return C;
	}
}

/**
 * @brief  Divides an Array1D A by a scalar b, element by element 
 * and creates a new Array1D C.
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
Array1D<T>  operator/(const Array1D<T> &A, const T &b)
{
	int n = A.dim1();
    Array1D<T> C(n);
	for (int i=0; i<n; i++)
	{
		C[i] = A[i]/b;
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
 * @return Array1D<T> A = A + B. 
 *
 * Example:
 * @code{.cpp}
 * A += B;
 * @endcode
 */
template <class T>
Array1D<T>&  operator+=(Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() == n)
	{
		for (int i=0; i<n; i++)
		{
				A[i] += B[i];
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
 * @return Array1D<T> A = A - B. 
 *
 * Example:
 * @code{.cpp}
 * A -= B;
 * @endcode
 */
template <class T>
Array1D<T>&  operator-=(Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() == n)
	{
		for (int i=0; i<n; i++)
		{
				A[i] -= B[i];
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
 * @return Array1D<T> A = A * B. 
 *
 * Example:
 * @code{.cpp}
 * A *= B;
 * @endcode
 */
template <class T>
Array1D<T>&  operator*=(Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() == n)
	{
		for (int i=0; i<n; i++)
		{
				A[i] *= B[i];
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
Array1D<T>&  operator*=(Array1D<T> &A, const T &b)
{
	int n = A.dim1();
		for (int i=0; i<n; i++)
		{
				A[i] *= b;
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
 * @return Array1D<T> A = A / B. 
 *
 * Example:
 * @code{.cpp}
 * A /= B;
 * @endcode
 */
template <class T>
Array1D<T>&  operator/=(Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() == n)
	{
		for (int i=0; i<n; i++)
		{
				A[i] /= B[i];
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
 * @return Array1D<T> A = A/b 
 *
 * Example:
 * @code{.cpp}
 * A /= b;
 * @endcode
 */
template <class T>
Array1D<T>&  operator/=(Array1D<T> &A, const T &b)
{
	int n = A.dim1();
		for (int i=0; i<n; i++)
		{
				A[i] /= b;
		}
	return A;
}

/* ........................ extended functions ......................*/

/**
 * @brief Compute the norm
 * 
 * @tparam T 
 * @param A 
 * @return double 
 */
template <class T>
double norm(const Array1D<T> &A)
{
	int n = A.dim1();

	double sum = 0.0;
	for (int i=0; i<n; i++)
		sum +=  abs(A[i])*abs(A[i]);
	return sqrt(sum);
}

/**
 * @brief  Computes A^B
 * 
 * @param A an Array1D
 * @param B an Array1D
 * @return Array1D<T>  A^B
 */
template <class T>
Array1D<T>  cross(const Array1D<T> &A, const Array1D<T> &B)
{
    int n = A.dim1();
    int m = B.dim1();
    
    if (n != 3 ||  m != 3 )
        return Array1D<T>();
    else
    {
        Array1D<T> C(n);
        C[0] = A[1]*B[2] - A[2]*B[1];
        C[1] = A[2]*B[0] - A[0]*B[2];
        C[2] = A[0]*B[1] - A[1]*B[0];
            return C;
    }
}


/**
 * @brief  Computes the dot product A.B
 * 
 * @param A an Array1D
 * @param B an Array1D
 * @return <T>  sum(Ai*Bi)
 */
template <class T>
T dot_product(const Array1D<T> &A, const Array1D<T> &B)
{
    int n = A.dim1();
    int m = B.dim1();
    
    if (n != 3 ||  m != 3 )
        return Array1D<T>();
    else
    {
		T sum = 0; 
		for (int i=0; i<n; i++)
			sum += A[i] * B[i];

        return sum;
    }
}

} // namespace TNT

#endif
