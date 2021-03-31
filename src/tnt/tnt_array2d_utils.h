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

#include <iostream>

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
std::ostream& operator<<(std::ostream &s, const Array2D<T> &A)
{
    int M=A.dim1();
    int N=A.dim2();

    s << "Dimensions: " << M << " " << N << "\n";
	s << "Components: \n";
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
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
std::istream& operator>>(std::istream &s, Array2D<T> &A)
{

    int M, N;

    s >> M >> N;

	Array2D<T> B(M,N);

    for (int i=0; i<M; i++)
        for (int j=0; j<N; j++)
        {
            s >>  B[i][j];
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

	if (B.dim1() != m ||  B.dim2() != n )
		return Array2D<T>();

	else
	{
		Array2D<T> C(m,n);

		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
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

	if (B.dim1() != m ||  B.dim2() != n )
		return Array2D<T>();

	else
	{
		Array2D<T> C(m,n);

		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
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

	if (B.dim1() != m ||  B.dim2() != n )
		return Array2D<T>();

	else
	{
		Array2D<T> C(m,n);

		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
				C[i][j] = A[i][j] * B[i][j];
		}
		return C;
	}
}




template <class T>
Array2D<T> operator/(const Array2D<T> &A, const Array2D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();

	if (B.dim1() != m ||  B.dim2() != n )
		return Array2D<T>();

	else
	{
		Array2D<T> C(m,n);

		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
				C[i][j] = A[i][j] / B[i][j];
		}
		return C;
	}
}





template <class T>
Array2D<T>&  operator+=(Array2D<T> &A, const Array2D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();

	if (B.dim1() == m ||  B.dim2() == n )
	{
		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
				A[i][j] += B[i][j];
		}
	}
	return A;
}



template <class T>
Array2D<T>&  operator-=(Array2D<T> &A, const Array2D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();

	if (B.dim1() == m ||  B.dim2() == n )
	{
		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
				A[i][j] -= B[i][j];
		}
	}
	return A;
}



template <class T>
Array2D<T>&  operator*=(Array2D<T> &A, const Array2D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();

	if (B.dim1() == m ||  B.dim2() == n )
	{
		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
				A[i][j] *= B[i][j];
		}
	}
	return A;
}





template <class T>
Array2D<T>&  operator/=(Array2D<T> &A, const Array2D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();

	if (B.dim1() == m ||  B.dim2() == n )
	{
		for (int i=0; i<m; i++)
		{
			for (int j=0; j<n; j++)
				A[i][j] /= B[i][j];
		}
	}
	return A;
}

/**
    Matrix Multiply:  compute C = A*B, where C[i][j]
    is the dot-product of row i of A and column j of B.


    @param A an (m x n) array
    @param B an (n x k) array
    @return the (m x k) array A*B, or a null array (0x0)
        if the matrices are non-conformant (i.e. the number
        of columns of A are different than the number of rows of B.)


*/
template <class T>
Array2D<T> matmult(const Array2D<T> &A, const Array2D<T> &B)
{
    if (A.dim2() != B.dim1())
        return Array2D<T>();

    int M = A.dim1();
    int N = A.dim2();
    int K = B.dim2();

    Array2D<T> C(M,K);

    for (int i=0; i<M; i++)
        for (int j=0; j<K; j++)
        {
            T sum = 0;

            for (int k=0; k<N; k++)
                sum += A[i][k] * B [k][j];

            C[i][j] = sum;
        }

    return C;

}

} // namespace TNT

#endif
