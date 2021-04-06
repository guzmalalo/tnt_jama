#ifndef JAMA_CHOLESKY_H
#define JAMA_CHOLESKY_H

#include "tnt_math_utils.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"

namespace JAMA
{
  /**
 * @brief For a symmetric, positive definite matrix A, this function
 * computes the Cholesky factorization, i.e. it computes a lower
 * triangular matrix L such that A = L*L'.If the matrix is not symmetric
 * or positive definite, the functioncomputes only a partial decomposition.
 * This can be tested with the is_spd() flag.
 * 
 * <p>Typical usage looks like:
 * 
 * @code{cpp}
	Array2D<double> A(n,n);
	Array2D<double> L;
	Cholesky<double> chol(A);

	if (chol.is_spd())
		L = chol.getL();	
  else
		cout << "factorization was not complete.\n";
 * @endcode
 * 
 * @tparam T A real data type
 * 
 * @note (Adapted from JAMA, a Java Matrix Library, developed by jointly 
 * by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
 */
  template <class T>
  class Cholesky
  {
    ///< Array for internal storage of decomposition
    TNT::Array2D<T> L_; // lower triangular factor

    ///< Symmetric and positive definite flag.
    int isspd; // 1 if matrix to be factored was SPD

  public:
    // Default Cholesky constructor
    Cholesky();

    // Cholesky algorithm for symmetric and positive definite matrix.
    Cholesky(const TNT::Array2D<T> &A);

    // Triangular factor
    TNT::Array2D<T> getL() const;

    /// Solve A*X = B
    TNT::Array1D<T> solve(const TNT::Array1D<T> &B);
    TNT::Array2D<T> solve(const TNT::Array2D<T> &B);

    // is the matrix symmetric and positive definite?
    int is_spd() const;
  };

  /**
 * @brief Construct a new Cholesky<T> object, and initialize 
 * 
 */
  template <class T>
  Cholesky<T>::Cholesky() : L_(0, 0), isspd(0) {}

  /**
 * @brief Verify if original matrix to be factored was symmetric 
 * positive-definite (SPD).
 * 
 * @tparam T 
 * @return int 1, if true , 0 if false
 */
  template <class T>
  int Cholesky<T>::is_spd() const
  {
    return isspd;
  }

  /**
 * @brief Return triangular factor, such that L*L'=A.
 * 
 * @return An Array2D<T> with the values of L
 */
  template <class T>
  TNT::Array2D<T> Cholesky<T>::getL() const
  {
    return L_;
  }

  /**
 * @brief Constructs a lower triangular matrix L, such that L*L'= A. 
 * If A is not symmetric positive-definite (SPD), only a partial factorization
 * is performed. If is_spd()evaluate true (1) then the factorization was successful.
 *  
 * @param A Input Array2D 
 */
  template <class T>
  Cholesky<T>::Cholesky(const TNT::Array2D<T> &A)
  {
    int m = A.dim1();
    int n = A.dim2();

    isspd = (m == n);

    if (m != n)
    {
      L_ = TNT::Array2D<T>();
      return;
    }

    L_ = TNT::Array2D<T>(n, n);

    // Main loop.
    for (int j = 0; j < n; j++)
    {
      T d(0.0);
      for (int k = 0; k < j; k++)
      {
        T s(0.0);
        for (int i = 0; i < k; i++)
        {
          s += L_[k][i] * L_[j][i];
        }
        L_[j][k] = s = (A[j][k] - s) / L_[k][k];
        d = d + s * s;
        isspd = isspd && (A[k][j] == A[j][k]);
      }
      d = A[j][j] - d;
      isspd = (isspd && (d > 0.0));
      L_[j][j] = std::sqrt(d > 0.0 ? d : 0.0);
      for (int k = j + 1; k < n; k++)
      {
        L_[j][k] = 0.0;
      }
    }
  }

  /**
 * @brief Solve a linear system A*x = b, using the previously 
 * computed cholesky factorization of A: L*L'.
 * 
 * @tparam T 
 * @param b An Array1D of dimension n = A.dim1() (Number of  rows of A)
 * @return TNT::Array1D<T> x so that L*L'*x = b. If b is nonconformat, 
 * or if A was not symmetric positive definite, a null (0x0) array is 
 * returned.
 */
  template <class T>
  TNT::Array1D<T> Cholesky<T>::solve(const TNT::Array1D<T> &b)
  {
    int n = L_.dim1();
    if (b.dim1() != n)
      return TNT::Array1D<T>();

    TNT::Array1D<T> x = b.copy();

    // Solve L*y = b;
    for (int k = 0; k < n; k++)
    {
      for (int i = 0; i < k; i++)
        x[k] -= x[i] * L_[k][i];
      x[k] /= L_[k][k];
    }

    // Solve L'*X = Y;
    for (int k = n - 1; k >= 0; k--)
    {
      for (int i = k + 1; i < n; i++)
        x[k] -= x[i] * L_[i][k];
      x[k] /= L_[k][k];
    }

    return x;
  }

  /**
 * @brief Solve a linear system A*X = B, using the previously computed
 * cholesky factorization of A: L*L'.
 * 
 * @param B  A Matrix with as many rows as A and any number of columns.
 * @return TNT::Array2D<T> X so that L*L'*X = B.  If B is nonconformat,
 * or if A was not symmetric positive definite, a null (0x0)	array is
 * returned. 
 */
  template <class T>
  TNT::Array2D<T> Cholesky<T>::solve(const TNT::Array2D<T> &B)
  {
    int n = L_.dim1();
    if (B.dim1() != n)
      return TNT::Array2D<T>();

    TNT::Array2D<T> X = B.copy();
    int nx = B.dim2();

    // Solve L*y = b;
    for (int j = 0; j < nx; j++)
    {
      for (int k = 0; k < n; k++)
      {
        for (int i = 0; i < k; i++)
          X[k][j] -= X[i][j] * L_[k][i];
        X[k][j] /= L_[k][k];
      }
    }

    // Solve L'*X = Y;
    for (int j = 0; j < nx; j++)
    {
      for (int k = n - 1; k >= 0; k--)
      {
        for (int i = k + 1; i < n; i++)
          X[k][j] -= X[i][j] * L_[i][k];
        X[k][j] /= L_[k][k];
      }
    }

    return X;
  }

}
// namespace JAMA

#endif
// JAMA_CHOLESKY_H