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

#ifndef TNT_ARRAY1D_H
#define TNT_ARRAY1D_H

#include "tnt_i_refvec.h"

namespace TNT
{

	/**
 * @brief One-dimensional numerical array which
	looks like a conventional C multiarray. 
	Storage corresponds to C (row-major) ordering.
	Elements are accessed via A[i] notation for 0-based indexing,
	and A(i) for 1-based indexing.. 
	
	<p>
	Array assignment is by reference (i.e. shallow assignment).
	That is, B=A implies that the A and B point to the
	same array, so modifications to the elements of A
	will be reflected in B. If an independent copy
	is required, then B = A.copy() can be used.  Note
	that this facilitates returning arrays from functions
	without relying on compiler optimizations to eliminate
	extensive data copying.
 * 
 * @tparam T   Used to determined the data type of array entries.  This is most
  commonly used when requiring scalar temporaries in templated algorithms
  that have TNT arrays as input.  For example,
  @code{.cpp}
  template <class Array1D>
  void foo(Array1D &A)
  {
    A::value_type first_entry = A[0];
    ...
  }
  @endcode
*/
	template <class T>
	class Array1D
	{

	private:
		/* ... */
		i_refvec<T> v_;
		int n_;
		T *data_; /* this normally points to v_.begin(), but
                             * could also point to a portion (subvector)
							 * of v_.
                            */

		void copy_(T *p, const T *q, int len) const;
		void set_(T *begin, T *end, const T &val);

	public:
		/**
 * @brief Used to determined the data type of array entries.
 * 
 */
		typedef T value_type;

		Array1D();
		explicit Array1D(int n);
		Array1D(int n, const T &a);
		Array1D(int n, T *a);
		inline Array1D(const Array1D &A);
		inline operator T *();
		inline operator const T *();
		inline Array1D &operator=(const T &a);
		inline Array1D &operator=(const Array1D &A);
		inline Array1D &ref(const Array1D &A);
		Array1D copy() const;
		Array1D &inject(const Array1D &A);
		inline T &operator[](int i);
		inline T &operator()(int i);
		inline const T &operator[](int i) const;
		inline const T &operator()(int i) const;
		inline int dim1() const;
		inline int dim() const;
		~Array1D();

		/* ... extended interface ... */

		inline int ref_count() const;
		inline Array1D<T> subarray(int i0, int i1);
	};

	/**
 * @brief Construct a new Array 1D<T> without initialization.
 * 
 * @tparam T Type
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A;
 * @endcode
 * 
 */
	template <class T>
	Array1D<T>::Array1D() : v_(), n_(0), data_(0) {}

	/**
 * @brief Construct a new Array1D<T> A from another Array1D B. A and B point to the same data (clone operation).
 * 
 * @tparam T Type of data
 * @param A input Array1D<T> to be cloned
 *
 * Example :
 * @code{.cpp}
 * TNT::Array1D<double> A(42, 1.0);
 * TNT::Array1D<double> C(A);
 * @endcode
 */
	template <class T>
	Array1D<T>::Array1D(const Array1D<T> &A) : v_(A.v_), n_(A.n_),
											   data_(A.data_)
	{
#ifdef TNT_DEBUG
		std::cout << "Created Array1D(const Array1D<T> &A) \n";
#endif
	}

	/**
 * @brief Construct a new Array1D of size n
 * @tparam T Type of the array data
 * @param n dimension of the array
 *
 * Example :
 * @code{.cpp}
 * 	TNT::Array1D<double> A(42);
 * @endcode
 */
	template <class T>
	Array1D<T>::Array1D(int n) : v_(n), n_(n), data_(v_.begin())
	{
#ifdef TNT_DEBUG
		std::cout << "Created Array1D(int n) \n";
#endif
	}

	/**
 * @brief Construct a new Array1D of size n with an initial value of n  
 * @tparam T Type of the data
 * @param n size of the array
 * @param val initial value for all the members
 *
 * Example :
 * @code{.cpp}
 * TNT::Array1D<double> A(42, 1.0);
 * @endcode
 */
	template <class T>
	Array1D<T>::Array1D(int n, const T &val) : v_(n), n_(n), data_(v_.begin())
	{
#ifdef TNT_DEBUG
		std::cout << "Created Array1D(int n, const T& val) \n";
#endif
		set_(data_, data_ + n, val);
	}

	/**
 * @brief Construct a new Array1D from an raw array of doubles of size n
 * (cloning operation)  
 * 
 * @tparam T type of the input array
 * @param n size of the array
 * @param a Input raw array of size n an type T
 *
 * Example:
 * @code{.cpp}
 * double ra[42] = {1.0};
 * TNT::Array1D<double> A(42, *ra);
 * @endcode
 *
 * @rst
 .. note:: 
 	The storage for this pre-existing array will never be destroyed by TNT
	because the referece is set to 0, this information is propaged to other
	cloning operations. So reference is reference_count is always 0.
 * @endrst
 */
	template <class T>
	Array1D<T>::Array1D(int n, T *a) : v_(a), n_(n), data_(v_.begin())
	{
#ifdef TNT_DEBUG
		std::cout << "Created Array1D(int n, T* a) \n";
#endif
	}

	/**
 * @brief Conversion operator function to raw arrays,
 * data is cloned
 * 
 * @tparam T Type of data
 * @return the adresse to the first member of the array
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42, 1.0);
 * double *ra = A;
 * @endcode
 *
 * @rst
 .. warning:: 
 	The raw array information is erased after the Array1D destruction because
 	the cloning operation. Only for local operations. 
 * @endrst
 */
	template <class T>
	inline Array1D<T>::operator T *()
	{
		return &(v_[0]);
	}

	/**
 * @brief Conversion operator function to constant raw arrays,
 * data is cloned
 * 
 * @tparam T Type of data
 * @return the adresse to the first member of the array
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42, 1.0);
 * const double *ra = A;
 * @endcode
 */
	template <class T>
	inline Array1D<T>::operator const T *()
	{
		return &(v_[0]);
	}

	/**
 * @brief Overload the acces operator [] 
 * 
 * @tparam T Type fo data
 * @param i index value (must start at 0 and end at dim-1)
 * @return T& Value of type T at index i 
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42, 1.0);
 * double first_number = A[0];
 * double last_number  = A[41];
 * @endcode
 */
	template <class T>
	inline T &Array1D<T>::operator[](int i)
	{
#ifdef TNT_BOUNDS_CHECK
		assert(i >= 0);
		assert(i < n_);
#endif
		return data_[i];
	}

	/**
 * @brief Overload the acces operator [] 
 * 
 * @tparam T Type fo data
 * @param i index value (must start at 0 and end at dim-1)
 * @return T& Constant value of type T at index i 
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42, 1.0);
 * const double first_number = A[0];
 * const double last_number  = A[41];
 * @endcode
 */
	template <class T>
	inline const T &Array1D<T>::operator[](int i) const
	{
#ifdef TNT_BOUNDS_CHECK
		assert(i >= 0 && "Index is negative");
		assert(i < n_ && "Index is out of bonds");
#endif
		return data_[i];
	}

	/**
 * @brief Overload the acces operator () 
 * 
 * @tparam T Type fo data
 * @param i index value (must start at 1 and end at dim)
 * @return T& Value of type T at index i-1 
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42, 1.0);
 * double first_number = A(1);
 * double last_number  = A(42);
 * @endcode
 */
	template <class T>
	inline T &Array1D<T>::operator()(int i)
	{
#ifdef TNT_BOUNDS_CHECK
		assert(i > 0 && "Index is negative");
		assert(i <= n_ && "Index is out of bonds");
#endif
		return data_[i - 1];
	};

	/**
 * @brief Overload the acces operator () 
 * 
 * @tparam T Type fo data
 * @param i index value (must start at 1 and end at dim)
 * @return T& Constant value of type T at index i-1 
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42, 1.0);
 * const double first_number = A(1);
 * const double last_number  = A(42);
 * @endcode
 */
	template <class T>
	inline const T &Array1D<T>::operator()(int i) const
	{
#ifdef TNT_BOUNDS_CHECK
		assert(i > 0 && "Index is negative");
		assert(i <= n_ && "Index is out of bonds");
#endif
		return data_[i - 1];
	}

	/**
 * @brief Initialaize all members af the array to a given value. 
 * Overload the operator =
 * 
 * @tparam T Type fo data
 * @param a input value 
 * @return Array1D<T>&  this Array
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42);
 * A = 1.0;
 * @endcode
 */
	template <class T>
	Array1D<T> &Array1D<T>::operator=(const T &a)
	{
		set_(data_, data_ + n_, a);
		return *this;
	}

	/**
 * @brief Make a copy of the current array to a different one. Avoid a
 * cloning operation.
 * 
 * @tparam T Type of data
 * @return Array1D<T> A copy of the Array
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42,1.0);
 * TNT::Array1D<double> B(A.copy()); 
 * @endcode
 */
	template <class T>
	Array1D<T> Array1D<T>::copy() const
	{
		Array1D A(n_);
		copy_(A.data_, data_, n_);

		return A;
	}

	/**
 * @brief Inject the values of an array A into another array B
 * of the same size. The values are copied (not cloning).
 * 
 * @tparam T Data type of the array
 * @param A Input Array, must be the same size of *this
 * @return Array1D<T>& A pointer to the *this array 
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> A(42,1.0);
 * TNT::Array1D<double> B(42,0.0); 
 * B.inject(A);
 * @endcode
 */
	template <class T>
	Array1D<T> &Array1D<T>::inject(const Array1D &A)
	{
		if (A.n_ == n_)
			copy_(data_, A.data_, n_);

		return *this;
	}

	/**
 * @brief Change the reference of one array B to another array B,
 * if B has only one reference the data is destroyed, if not the
 * reference counter is decreased.
 * 
 * @tparam T Type of data
 * @param A A reference to an Array1D
 * @return Array1D<T>& *this with the same reference of A
 *
 * Examples:
 * @code{.cpp}
 *
 * TNT::Array1D<double> A(a_size_, init_value_);
 * TNT::Array1D<double> Z(a_size_, init_value_);
 * TNT::Array1D<double> B(13);
 * TNT::Array1D<double> C;
 * TNT::Array1D<double> *D = new TNT::Array1D<double>(1);
 *
 *  B.ref(A);  //  B is now a reference to A
 *  C.ref(A);  //  C is now a reference to A
 *  D->ref(A); //  D is now a reference to A
 * @endcode
 */
	template <class T>
	Array1D<T> &Array1D<T>::ref(const Array1D<T> &A)
	{
		if (this != &A)
		{
			v_ = A.v_; /* operator= handles the reference counting. */
			n_ = A.n_;
			data_ = A.data_;
		}
		return *this;
	}

	/**
 * @brief Overload the = operator. Use the method Array1D<T>::ref .
 * This means that "=" performs a clone operation not a copy operation.
 * 
 * @tparam T Type of the data
 * @param A A reference to an Array1D
 * @return Array1D<T>& Array1D<T>& *this with the same reference of A
 *
 * Examples:
 * @code{.cpp}
 *
 * TNT::Array1D<double> A(a_size_, init_value_);
 * TNT::Array1D<double> Z(a_size_, init_value_);
 * TNT::Array1D<double> B(13);
 * TNT::Array1D<double> C;
 * TNT::Array1D<double> *D = new TNT::Array1D<double>(1);
 *
 *  B = ref(A);  //  B is now a reference to A
 *  C = ref(A);  //  C is now a reference to A
 *  (*D) = ref(A); // D is now a reference to A
 * @endcode
 */
	template <class T>
	Array1D<T> &Array1D<T>::operator=(const Array1D<T> &A)
	{
		return ref(A);
	}

	/**
 * @brief Returns the dimension of the Array
 * 
 * @tparam T Type of data
 * @return int the size of the Array
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> B(42);
 * int size = B.dim() // equal to 42
 * @endcode
 */
	template <class T>
	inline int Array1D<T>::dim1() const { return n_; }

	/**
 * @brief Returns the dimension of the Array
 * 
 * @tparam T Type of data
 * @return int the size of the Array
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> B(42);
 * int size = B.dim1() // equal to 42
 * @endcode
 */
	template <class T>
	inline int Array1D<T>::dim() const { return n_; }

	/**
 * @brief Destroy the Array1D<T> object
 * 
 * @tparam T Type of data
 *
 * Example:
 * @code{.cpp}
 * TNT::Array1D<double> B(42);
 * int size = B.dim1() // equal to 42
 * @endcode
 */
	template <class T>
	Array1D<T>::~Array1D() {}

	/* ............................ extended interface ......................*/

	/**
 * @brief Return the number of reference to *this array
 * 
 * @tparam T Type of data
 * @return int The number of references
 *
 * Examples:
 * @code{.cpp}
 *
 * TNT::Array1D<double> A(a_size_, init_value_);
 * TNT::Array1D<double> Z(a_size_, init_value_);
 * TNT::Array1D<double> B(13);
 * TNT::Array1D<double> C;
 * TNT::Array1D<double> *D = new TNT::Array1D<double>(1);
 *
 * B = ref(A);  //  B is now a reference to A
 * C = ref(A);  //  C is now a reference to A
 * (*D) = ref(A); // D is now a reference to A
 * 
 * int num_ref = A.ref_count() // equal to 4
 * @endcode
 */
	template <class T>
	inline int Array1D<T>::ref_count() const
	{
		return v_.ref_count();
	}

	/**
 * @brief Creates an Array1B B composed by a REFERENCE to a subarray of this 
 * between the index i0 and i1
 * 
 * @tparam T Type of data
 * @param i0 start index
 * @param i1 end index
 * @return Array1D<T> A reference subarray.  
 *
 * Examples:
 * @code{.cpp}
 *
 * TNT::Array1D<double> A(a_size_, init_value_);
 *
 * TNT::Array1D<double> B = A.subarray(0, 9);
 * 
 * bool a = (&A[i]==&B[i]) // true 
 * @endcode
 */
	template <class T>
	inline Array1D<T> Array1D<T>::subarray(int i0, int i1)
	{
		if ((i0 >= 0) && (i1 < n_) && (i0 <= i1))
		{
			Array1D<T> X(*this); /* create a new instance of this array. */
			X.n_ = i1 - i0 + 1;
			X.data_ += i0;

			return X;
		}
		else
		{
			return Array1D<T>();
		}
	}

	/* private internal functions */

	template <class T>
	void Array1D<T>::set_(T *begin, T *end, const T &a)
	{
		for (T *p = begin; p < end; p++)
			*p = a;
	}

	template <class T>
	void Array1D<T>::copy_(T *p, const T *q, int len) const
	{
		T *end = p + len;
		while (p < end)
			*p++ = *q++;
	}

} /* namespace TNT */

#endif
/* TNT_ARRAY1D_H */
