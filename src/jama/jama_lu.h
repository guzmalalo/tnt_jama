#ifndef JAMA_LU_H
#define JAMA_LU_H

#include "tnt_math_utils.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"

namespace JAMA
{

  /**
 * @brief LU Decomposition.
   <P>
   For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
   unit lower triangular matrix L, an n-by-n upper triangular matrix U,
   and a permutation vector piv of length m so that A(piv,:) = L*U.
   If m < n, then L is m-by-m and U is m-by-n.
   <P>
   The LU decomposition with pivoting always exists, even if the matrix is
   singular, so the constructor will never fail.  The primary use of the
   LU decomposition is in the solution of square systems of simultaneous
   linear equations.  This will fail if isNonsingular() returns false.
 * 
 * @tparam T A real type data:  
 * - int
 * - double
 * - float
 * - etc ..
 */
  template <class T>
  class LU
  {
  private:
    ///< Array for internal storage of decomposition.
    TNT::Array2D<T> LU_;

    ///< m column dimension.
    ///< n row dimension.
    ///< p pivot sign.
    int m, n, pivsign;

    ///< Internal storage of pivot vector.
    TNT::Array1D<int> piv;

    TNT::Array2D<T> permute_copy(const TNT::Array2D<T> &A,
                                 const TNT::Array1D<int> &piv, int j0, int j1)
    {
      int piv_length = piv.dim();

      TNT::Array2D<T> X(piv_length, j1 - j0 + 1);

      for (int i = 0; i < piv_length; i++)
        for (int j = j0; j <= j1; j++)
          X[i][j - j0] = A[piv[i]][j];

      return X;
    }

    TNT::Array1D<T> permute_copy(const TNT::Array1D<T> &A,
                                 const TNT::Array1D<int> &piv)
    {
      int piv_length = piv.dim();
      if (piv_length != A.dim())
        return TNT::Array1D<T>();

      TNT::Array1D<T> x(piv_length);

      for (int i = 0; i < piv_length; i++)
        x[i] = A[piv[i]];

      return x;
    }

  public:
    /**
   * @brief LU decomposition constructor.
   * @param A Rectangular matrix
   * @return LU Decomposition object to access L, U and piv.
   */
    LU(const TNT::Array2D<T> &A) : LU_(A.copy()), m(A.dim1()), n(A.dim2()),
                                   piv(A.dim1())
    {

      // Use a "left-looking", dot-product, Crout/Doolittle algorithm.

      for (int i = 0; i < m; i++)
      {
        piv[i] = i;
      }
      pivsign = 1;
      T *LUrowi = 0;

      TNT::Array1D<T> LUcolj(m);

      // Outer loop.
      for (int j = 0; j < n; j++)
      {
        // Make a copy of the j-th column to localize references.
        for (int i = 0; i < m; i++)
        {
          LUcolj[i] = LU_[i][j];
        }

        // Apply previous transformations.
        for (int i = 0; i < m; i++)
        {
          LUrowi = LU_[i];

          // Most of the time is spent in the following dot product.
          int kmax = std::min(i, j);
          double s = 0.0;
          for (int k = 0; k < kmax; k++)
          {
            s += LUrowi[k] * LUcolj[k];
          }

          LUrowi[j] = LUcolj[i] -= s;
        }

        // Find pivot and exchange if necessary.
        int p = j;
        for (int i = j + 1; i < m; i++)
        {
          if (std::abs(LUcolj[i]) > std::abs(LUcolj[p]))
          {
            p = i;
          }
        }
        if (p != j)
        {
          int k = 0;
          for (k = 0; k < n; k++)
          {
            double t = LU_[p][k];
            LU_[p][k] = LU_[j][k];
            LU_[j][k] = t;
          }
          k = piv[p];
          piv[p] = piv[j];
          piv[j] = k;
          pivsign = -pivsign;
        }

        // Compute multipliers.
        if ((j < m) && (LU_[j][j] != 0.0))
        {
          for (int i = j + 1; i < m; i++)
          {
            LU_[i][j] /= LU_[j][j];
          }
        }
      }
    }

    /**
 * @brief Check is the matrix is non singular.
 * 
 * @return int 1 (true)  if upper triangular factor U (and hence A) 
 * is nonsingular, 0 otherwise.
 */
    int isNonsingular()
    {
      for (int j = 0; j < n; j++)
      {
        if (LU_[j][j] == 0)
          return 0;
      }
      return 1;
    };

    /**
 * @brief Return lower triangular factor
 * 
 * @return TNT::Array2D<T> L
 */
    TNT::Array2D<T> getL()
    {
      TNT::Array2D<T> L_(m, n);
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
        {
          if (i > j)
          {
            L_[i][j] = LU_[i][j];
          }
          else if (i == j)
          {
            L_[i][j] = 1.0;
          }
          else
          {
            L_[i][j] = 0.0;
          }
        }
      }
      return L_;
    }

    /**
 * @brief Return upper triangular factor
 * 
 * @return TNT::Array2D<T>  U portion of LU factorization.
 */
    TNT::Array2D<T> getU()
    {
      TNT::Array2D<T> U_(m, n);
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
        {
          if (i <= j)
          {
            U_[i][j] = LU_[i][j];
          }
          else
          {
            U_[i][j] = 0.0;
          }
        }
      }
      return U_;
    }

    /**
 * @brief Get the Pivot permutation vector
 * 
 * @return TNT::Array1D<int> piv
 */
    TNT::Array1D<int> getPivot()
    {
      return piv;
    }

    /**
 * @brief Compute determinant using LU factors.
 * 
 * @return T  determinant of A, or 0 if A is not square.
 */
    T det()
    {
      if (m != n)
      {
        return T(0);
      }
      T d = T(pivsign);
      for (int j = 0; j < n; j++)
      {
        d *= LU_[j][j];
      }
      return d;
    }

    /**
 * @brief Solve A*X = B
 * 
 * @param B A Matrix with as many rows as A and any number of columns.
 * @return TNT::Array2D<T> X so that L*U*X = B(piv,:), if B is nonconformant, 
 * 0x0 (null) array.
 */
    TNT::Array2D<T> solve(const TNT::Array2D<T> &B)
    {

      /* Dimensions: A is mxn, X is nxk, B is mxk */

      if (B.dim1() != m)
      {
        return TNT::Array2D<T>();
      }
      if (!isNonsingular())
      {
        return TNT::Array2D<T>();
      }

      // Copy right hand side with pivoting
      int nx = B.dim2();

      TNT::Array2D<T> X = permute_copy(B, piv, 0, nx - 1);

      // Solve L*Y = B(piv,:)
      for (int k = 0; k < n; k++)
      {
        for (int i = k + 1; i < n; i++)
        {
          for (int j = 0; j < nx; j++)
          {
            X[i][j] -= X[k][j] * LU_[i][k];
          }
        }
      }
      // Solve U*X = Y;
      for (int k = n - 1; k >= 0; k--)
      {
        for (int j = 0; j < nx; j++)
        {
          X[k][j] /= LU_[k][k];
        }
        for (int i = 0; i < k; i++)
        {
          for (int j = 0; j < nx; j++)
          {
            X[i][j] -= X[k][j] * LU_[i][k];
          }
        }
      }
      return X;
    }

    /**
 * @brief Solve A*x = b, where x and b are vectors of length equal	
 * to the number of rows in A.
 * 
 * @param b a vector (Array1D> of length equal to the first dimension
 * of A.
 * @return TNT::Array1D<T> x a vector (Array1D> so that L*U*x = b(piv), 
 * if B is nonconformant,	returns 0x0 (null) array.
 */
    TNT::Array1D<T> solve(const TNT::Array1D<T> &b)
    {

      /* Dimensions: A is mxn, X is nxk, B is mxk */

      if (b.dim1() != m)
      {
        return TNT::Array1D<T>();
      }
      if (!isNonsingular())
      {
        return TNT::Array1D<T>();
      }

      TNT::Array1D<T> x = permute_copy(b, piv);

      // Solve L*Y = B(piv)
      for (int k = 0; k < n; k++)
      {
        for (int i = k + 1; i < n; i++)
        {
          x[i] -= x[k] * LU_[i][k];
        }
      }

      // Solve U*X = Y;
      for (int k = n - 1; k >= 0; k--)
      {
        x[k] /= LU_[k][k];
        for (int i = 0; i < k; i++)
          x[i] -= x[k] * LU_[i][k];
      }

      return x;
    }

  }; /* class LU */

} /* namespace JAMA */

#endif
/* JAMA_LU_H */
