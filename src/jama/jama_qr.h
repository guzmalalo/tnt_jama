#ifndef JAMA_QR_H
#define JAMA_QR_H

#include "tnt_math_utils.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"

namespace JAMA
{

  /**
 * @brief Classical QR Decomposition.
 * For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
 * orthogonal matrix Q and an n-by-n upper triangular matrix R so that
 * A = Q*R.
 * 
 * The QR decomposition always exists, even if the matrix does not have
 * full rank, so the constructor will never fail.  The primary use of the
 * QR decomposition is in the least squares solution of nonsquare systems
 * of simultaneous linear equations.  This will fail if isFullRank()
 * returns 0 (false).
 * 
 * The Q and R factors can be retrived via the getQ() and getR()
 * methods. Furthermore, a solve() method is provided to find the
 * least squares solution of Ax=b using the QR factors.  
 * 
 * @note : Adapted from JAMA, a Java Matrix Library, developed by jointly 
	by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama.
 * 
 * @tparam T data type 
 */
  template <class T>
  class QR
  {

    ///< Array for internal storage of decomposition.
    TNT::Array2D<T> QR_;

    ///< Row and column dimensions
    int m, n;

    ///< Array for internal storage of diagonal of R.
    TNT::Array1D<T> Rdiag;

  public:
    /**
 * @brief Create a QR factorization object for A.
 * 
 * @param A rectangular (m>=n) matrix.
 */
    QR(const TNT::Array2D<T> &A) /* constructor */
    {
      QR_ = A.copy();
      m = A.dim1();
      n = A.dim2();
      Rdiag = TNT::Array1D<T>(n);
      int i = 0, j = 0, k = 0;

      // Main loop.
      for (k = 0; k < n; k++)
      {
        // Compute 2-norm of k-th column without under/overflow.
        T nrm = 0;
        for (i = k; i < m; i++)
        {
          nrm = TNT::hypot(nrm, QR_[i][k]);
        }

        if (nrm != 0.0)
        {
          // Form k-th Householder vector.
          if (QR_[k][k] < 0)
          {
            nrm = -nrm;
          }
          for (i = k; i < m; i++)
          {
            QR_[i][k] /= nrm;
          }
          QR_[k][k] += 1.0;

          // Apply transformation to remaining columns.
          for (j = k + 1; j < n; j++)
          {
            T s = 0.0;
            for (i = k; i < m; i++)
            {
              s += QR_[i][k] * QR_[i][j];
            }
            s = -s / QR_[k][k];
            for (i = k; i < m; i++)
            {
              QR_[i][j] += s * QR_[i][k];
            }
          }
        }
        Rdiag[k] = -nrm;
      }
    }

    /**
 * @brief Flag to denote the matrix is of full rank.
 * 
 * @return int 1 if matrix is full rank, 0 otherwise.
 */
    int isFullRank() const
    {
      for (int j = 0; j < n; j++)
      {
        if (Rdiag[j] == 0)
          return 0;
      }
      return 1;
    }

    /**
 * @brief Retreive the Householder vectors from QR factorization
 * 
 * @return TNT::Array2D<T> lower trapezoidal matrix whose columns define the reflections
 */
    TNT::Array2D<T> getHouseholder(void) const
    {
      TNT::Array2D<T> H(m, n);

      /* note: H is completely filled in by algorithm, so
	     initialization of H is not necessary.
	  */
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
        {
          if (i >= j)
          {
            H[i][j] = QR_[i][j];
          }
          else
          {
            H[i][j] = 0.0;
          }
        }
      }
      return H;
    }

    /**
 * @brief Return the upper triangular factor, R, of the QR factorization
 * 
 * @return TNT::Array2D<T> R
 */
    TNT::Array2D<T> getR() const
    {
      TNT::Array2D<T> R(n, n);
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          if (i < j)
          {
            R[i][j] = QR_[i][j];
          }
          else if (i == j)
          {
            R[i][j] = Rdiag[i];
          }
          else
          {
            R[i][j] = 0.0;
          }
        }
      }
      return R;
    }

    /**
 * @brief Generate and return the (economy-sized) orthogonal factor
 * 
 * @return TNT::Array2D<T>  Q the (economy-sized) orthogonal factor (Q*R=A).
 */
    TNT::Array2D<T> getQ() const
    {
      int i = 0, j = 0, k = 0;

      TNT::Array2D<T> Q(m, n);
      for (k = n - 1; k >= 0; k--)
      {
        for (i = 0; i < m; i++)
        {
          Q[i][k] = 0.0;
        }
        Q[k][k] = 1.0;
        for (j = k; j < n; j++)
        {
          if (QR_[k][k] != 0)
          {
            T s = 0.0;
            for (i = k; i < m; i++)
            {
              s += QR_[i][k] * Q[i][j];
            }
            s = -s / QR_[k][k];
            for (i = k; i < m; i++)
            {
              Q[i][j] += s * QR_[i][k];
            }
          }
        }
      }
      return Q;
    }

    /**
 * @brief Least squares solution of A*x = b
 * 
 * @param b m-length array (vector).
 * @return TNT::Array1D<T> x n-length array (vector) that minimizes 
 * the two norm of Q*R*X-B.
 * If B is non-conformant, or if QR.isFullRank() is false,
 * the routine returns a null (0-length) vector.
 */
    TNT::Array1D<T> solve(const TNT::Array1D<T> &b) const
    {
      if (b.dim1() != m) /* arrays must be conformant */
        return TNT::Array1D<T>();

      if (!isFullRank()) /* matrix is rank deficient */
      {
        return TNT::Array1D<T>();
      }

      TNT::Array1D<T> x = b.copy();

      // Compute Y = transpose(Q)*b
      for (int k = 0; k < n; k++)
      {
        T s = 0.0;
        for (int i = k; i < m; i++)
        {
          s += QR_[i][k] * x[i];
        }
        s = -s / QR_[k][k];
        for (int i = k; i < m; i++)
        {
          x[i] += s * QR_[i][k];
        }
      }
      // Solve R*X = Y;
      for (int k = n - 1; k >= 0; k--)
      {
        x[k] /= Rdiag[k];
        for (int i = 0; i < k; i++)
        {
          x[i] -= x[k] * QR_[i][k];
        }
      }

      /* return n x nx portion of X */
      TNT::Array1D<T> x_(n);
      for (int i = 0; i < n; i++)
        x_[i] = x[i];

      return x_;
    }

    /**
 * @brief Least squares solution of A*X = B
 * 
 * @param B  m x k Array (must conform).
 * @return TNT::Array2D<T> X n x k Array that minimizes the two norm of Q*R*X-B.
 * If B is non-conformant, or if QR.isFullRank() is false,
 * the routine returns a null (0x0) array.
 */
    TNT::Array2D<T> solve(const TNT::Array2D<T> &B) const
    {
      if (B.dim1() != m) /* arrays must be conformant */
        return TNT::Array2D<T>(0, 0);

      if (!isFullRank()) /* matrix is rank deficient */
      {
        return TNT::Array2D<T>(0, 0);
      }

      int nx = B.dim2();
      TNT::Array2D<T> X = B.copy();
      int i = 0, j = 0, k = 0;

      // Compute Y = transpose(Q)*B
      for (k = 0; k < n; k++)
      {
        for (j = 0; j < nx; j++)
        {
          T s = 0.0;
          for (i = k; i < m; i++)
          {
            s += QR_[i][k] * X[i][j];
          }
          s = -s / QR_[k][k];
          for (i = k; i < m; i++)
          {
            X[i][j] += s * QR_[i][k];
          }
        }
      }
      // Solve R*X = Y;
      for (k = n - 1; k >= 0; k--)
      {
        for (j = 0; j < nx; j++)
        {
          X[k][j] /= Rdiag[k];
        }
        for (i = 0; i < k; i++)
        {
          for (j = 0; j < nx; j++)
          {
            X[i][j] -= X[k][j] * QR_[i][k];
          }
        }
      }

      /* return n x nx portion of X */
      TNT::Array2D<T> X_(n, nx);
      for (i = 0; i < n; i++)
        for (j = 0; j < nx; j++)
          X_[i][j] = X[i][j];

      return X_;
    }
  };

}
// namespace JAMA

#endif
// JAMA_QR__H
