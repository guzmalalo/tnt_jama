#ifndef JAMA_SVD_H
#define JAMA_SVD_H

#include "tnt_math_utils.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"

namespace JAMA
{

  /** 
 * @brief Singular Value Decomposition.
 * 
 * For an m-by-n matrix A with m >= n, the singular value decomposition is
 * an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
 * an n-by-n orthogonal matrix V so that A = U*S*V'.
 * 
 * The singular values, sigma[k] = S[k][k], are ordered so that
 * sigma[0] >= sigma[1] >= ... >= sigma[n-1]. i.e  singular values of matrix A
 * in descending order.
 * 
 * @tparam T data type (real)
 */
  template <class T>
  class SVD
  {

    ///< Arrays for internal storage of U and V.
    TNT::Array2D<T> U, V;

    ///< Array for internal storage of singular values.
    TNT::Array1D<T> s;

    ///< Row and column dimensions.
    int m, n;

  public:
    /**
    * @brief Construct the singular value decomposition
    * Structure to access U, S and V.
    * 
    * @param Arg Rectangular matrix
    */
    SVD(const TNT::Array2D<T> &Arg)
    {
      // Derived from LINPACK code.
      // Initialize.
      m = Arg.dim1();
      n = Arg.dim2();
      int nu = std::min(m, n);
      s = TNT::Array1D<T>(std::min(m + 1, n));
      U = TNT::Array2D<T>(m, nu, T(0));
      V = TNT::Array2D<T>(n, n);
      TNT::Array1D<T> e(n);
      TNT::Array1D<T> work(m);
      TNT::Array2D<T> A(Arg.copy());
      int wantu = 1; /* boolean */
      int wantv = 1; /* boolean */
      int i = 0, j = 0, k = 0;

      // Reduce A to bi-diagonal form, storing the diagonal elements
      // in s and the super-diagonal elements in e.

      int nct = std::min(m - 1, n);
      int nrt = std::max(0, std::min(n - 2, m));
      for (k = 0; k < std::max(nct, nrt); k++)
      {
        if (k < nct)
        {

          // Compute the transformation for the k-th column and
          // place the k-th diagonal in s[k].
          // Compute 2-norm of k-th column without under/overflow.
          s[k] = 0;
          for (i = k; i < m; i++)
          {
            s[k] = TNT::hypot(s[k], A[i][k]);
          }
          if (s[k] != 0.0)
          {
            if (A[k][k] < 0.0)
            {
              s[k] = -s[k];
            }
            for (i = k; i < m; i++)
            {
              A[i][k] /= s[k];
            }
            A[k][k] += 1.0;
          }
          s[k] = -s[k];
        }
        for (j = k + 1; j < n; j++)
        {
          if ((k < nct) && (s[k] != 0.0))
          {

            // Apply the transformation.

            T t(0.0);
            for (i = k; i < m; i++)
            {
              t += A[i][k] * A[i][j];
            }
            t = -t / A[k][k];
            for (i = k; i < m; i++)
            {
              A[i][j] += t * A[i][k];
            }
          }

          // Place the k-th row of A into e for the
          // subsequent calculation of the row transformation.

          e[j] = A[k][j];
        }
        if (wantu & (k < nct))
        {

          // Place the transformation in U for subsequent back
          // multiplication.

          for (i = k; i < m; i++)
          {
            U[i][k] = A[i][k];
          }
        }
        if (k < nrt)
        {

          // Compute the k-th row transformation and place the
          // k-th super-diagonal in e[k].
          // Compute 2-norm without under/overflow.
          e[k] = 0;
          for (i = k + 1; i < n; i++)
          {
            e[k] = TNT::hypot(e[k], e[i]);
          }
          if (e[k] != 0.0)
          {
            if (e[k + 1] < 0.0)
            {
              e[k] = -e[k];
            }
            for (i = k + 1; i < n; i++)
            {
              e[i] /= e[k];
            }
            e[k + 1] += 1.0;
          }
          e[k] = -e[k];
          if ((k + 1 < m) & (e[k] != 0.0))
          {

            // Apply the transformation.

            for (i = k + 1; i < m; i++)
            {
              work[i] = 0.0;
            }
            for (j = k + 1; j < n; j++)
            {
              for (i = k + 1; i < m; i++)
              {
                work[i] += e[j] * A[i][j];
              }
            }
            for (j = k + 1; j < n; j++)
            {
              T t(-e[j] / e[k + 1]);
              for (i = k + 1; i < m; i++)
              {
                A[i][j] += t * work[i];
              }
            }
          }
          if (wantv)
          {

            // Place the transformation in V for subsequent
            // back multiplication.

            for (i = k + 1; i < n; i++)
            {
              V[i][k] = e[i];
            }
          }
        }
      }

      // Set up the final bi-diagonal matrix or order p.

      int p = std::min(n, m + 1);
      if (nct < n)
      {
        s[nct] = A[nct][nct];
      }
      if (m < p)
      {
        s[p - 1] = 0.0;
      }
      if (nrt + 1 < p)
      {
        e[nrt] = A[nrt][p - 1];
      }
      e[p - 1] = 0.0;

      // If required, generate U.

      if (wantu)
      {
        for (j = nct; j < nu; j++)
        {
          for (i = 0; i < m; i++)
          {
            U[i][j] = 0.0;
          }
          U[j][j] = 1.0;
        }
        for (k = nct - 1; k >= 0; k--)
        {
          if (s[k] != 0.0)
          {
            for (j = k + 1; j < nu; j++)
            {
              T t(0.0);
              for (i = k; i < m; i++)
              {
                t += U[i][k] * U[i][j];
              }
              t = -t / U[k][k];
              for (i = k; i < m; i++)
              {
                U[i][j] += t * U[i][k];
              }
            }
            for (i = k; i < m; i++)
            {
              U[i][k] = -U[i][k];
            }
            U[k][k] = 1.0 + U[k][k];
            for (i = 0; i < k - 1; i++)
            {
              U[i][k] = 0.0;
            }
          }
          else
          {
            for (i = 0; i < m; i++)
            {
              U[i][k] = 0.0;
            }
            U[k][k] = 1.0;
          }
        }
      }

      // If required, generate V.

      if (wantv)
      {
        for (k = n - 1; k >= 0; k--)
        {
          if ((k < nrt) & (e[k] != 0.0))
          {
            for (j = k + 1; j < nu; j++)
            {
              T t(0.0);
              for (i = k + 1; i < n; i++)
              {
                t += V[i][k] * V[i][j];
              }
              t = -t / V[k + 1][k];
              for (i = k + 1; i < n; i++)
              {
                V[i][j] += t * V[i][k];
              }
            }
          }
          for (i = 0; i < n; i++)
          {
            V[i][k] = 0.0;
          }
          V[k][k] = 1.0;
        }
      }

      // Main iteration loop for the singular values.

      int pp = p - 1;
      int iter = 0;
      T eps = std::numeric_limits<T>::epsilon();
      while (p > 0)
      {
        int k = 0;
        int kase = 0;

        // Here is where a test for too many iterations would go.

        // This section of the program inspects for
        // negligible elements in the s and e arrays.  On
        // completion the variables kase and k are set as follows.

        // kase = 1     if s(p) and e[k-1] are negligible and k<p
        // kase = 2     if s(k) is negligible and k<p
        // kase = 3     if e[k-1] is negligible, k<p, and
        //              s(k), ..., s(p) are not negligible (qr step).
        // kase = 4     if e(p-1) is negligible (convergence).

        for (k = p - 2; k >= -1; k--)
        {
          if (k == -1)
          {
            break;
          }
          if (std::abs(e[k]) <= eps * (std::abs(s[k]) + std::abs(s[k + 1])))
          {
            e[k] = 0.0;
            break;
          }
        }
        if (k == p - 2)
        {
          kase = 4;
        }
        else
        {
          int ks;
          for (ks = p - 1; ks >= k; ks--)
          {
            if (ks == k)
            {
              break;
            }
            T t((ks != p ? std::abs(e[ks]) : 0.) +
                (ks != k + 1 ? std::abs(e[ks - 1]) : 0.));
            if (std::abs(s[ks]) <= eps * t)
            {
              s[ks] = 0.0;
              break;
            }
          }
          if (ks == k)
          {
            kase = 3;
          }
          else if (ks == p - 1)
          {
            kase = 1;
          }
          else
          {
            kase = 2;
            k = ks;
          }
        }
        k++;

        // Perform the task indicated by kase.

        switch (kase)
        {

          // Deflate negligible s(p).

        case 1:
        {
          T f(e[p - 2]);
          e[p - 2] = 0.0;
          for (j = p - 2; j >= k; j--)
          {
            T t(TNT::hypot(s[j], f));
            T cs(s[j] / t);
            T sn(f / t);
            s[j] = t;
            if (j != k)
            {
              f = -sn * e[j - 1];
              e[j - 1] = cs * e[j - 1];
            }
            if (wantv)
            {
              for (i = 0; i < n; i++)
              {
                t = cs * V[i][j] + sn * V[i][p - 1];
                V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
                V[i][j] = t;
              }
            }
          }
        }
        break;

          // Split at negligible s(k).

        case 2:
        {
          T f(e[k - 1]);
          e[k - 1] = 0.0;
          for (j = k; j < p; j++)
          {
            T t(TNT::hypot(s[j], f));
            T cs(s[j] / t);
            T sn(f / t);
            s[j] = t;
            f = -sn * e[j];
            e[j] = cs * e[j];
            if (wantu)
            {
              for (i = 0; i < m; i++)
              {
                t = cs * U[i][j] + sn * U[i][k - 1];
                U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
                U[i][j] = t;
              }
            }
          }
        }
        break;

          // Perform one qr step.

        case 3:
        {

          // Calculate the shift.

          T scale = std::max(
              std::max(
                  std::max(
                      std::max(std::abs(s[p - 1]), std::abs(s[p - 2])),
                      std::abs(e[p - 2])),
                  std::abs(s[k])),
              std::abs(e[k]));
          T sp = s[p - 1] / scale;
          T spm1 = s[p - 2] / scale;
          T epm1 = e[p - 2] / scale;
          T sk = s[k] / scale;
          T ek = e[k] / scale;
          T b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
          T c = (sp * epm1) * (sp * epm1);
          T shift = 0.0;
          if ((b != 0.0) || (c != 0.0))
          {
            shift = std::sqrt(b * b + c);
            if (b < 0.0)
            {
              shift = -shift;
            }
            shift = c / (b + shift);
          }
          T f = (sk + sp) * (sk - sp) + shift;
          T g = sk * ek;

          // Chase zeros.

          for (j = k; j < p - 1; j++)
          {
            T t = TNT::hypot(f, g);
            T cs = f / t;
            T sn = g / t;
            if (j != k)
            {
              e[j - 1] = t;
            }
            f = cs * s[j] + sn * e[j];
            e[j] = cs * e[j] - sn * s[j];
            g = sn * s[j + 1];
            s[j + 1] = cs * s[j + 1];
            if (wantv)
            {
              for (i = 0; i < n; i++)
              {
                t = cs * V[i][j] + sn * V[i][j + 1];
                V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
                V[i][j] = t;
              }
            }
            t = TNT::hypot(f, g);
            cs = f / t;
            sn = g / t;
            s[j] = t;
            f = cs * e[j] + sn * s[j + 1];
            s[j + 1] = -sn * e[j] + cs * s[j + 1];
            g = sn * e[j + 1];
            e[j + 1] = cs * e[j + 1];
            if (wantu && (j < m - 1))
            {
              for (i = 0; i < m; i++)
              {
                t = cs * U[i][j] + sn * U[i][j + 1];
                U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
                U[i][j] = t;
              }
            }
          }
          e[p - 2] = f;
          iter = iter + 1;
        }
        break;

          // Convergence.

        case 4:
        {

          // Make the singular values positive.

          if (s[k] <= 0.0)
          {
            s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
            if (wantv)
            {
              for (i = 0; i <= pp; i++)
              {
                V[i][k] = -V[i][k];
              }
            }
          }

          // Order the singular values.

          while (k < pp)
          {
            if (s[k] >= s[k + 1])
            {
              break;
            }
            T t = s[k];
            s[k] = s[k + 1];
            s[k + 1] = t;
            if (wantv && (k < n - 1))
            {
              for (i = 0; i < n; i++)
              {
                t = V[i][k + 1];
                V[i][k + 1] = V[i][k];
                V[i][k] = t;
              }
            }
            if (wantu && (k < m - 1))
            {
              for (i = 0; i < m; i++)
              {
                t = U[i][k + 1];
                U[i][k + 1] = U[i][k];
                U[i][k] = t;
              }
            }
            k++;
          }
          iter = 0;
          p--;
        }
        break;
        }
      }
    }

    /**
 * @brief Return the left singular vectors
 * 
 * @param A 
 */
    TNT::Array2D<T> getU()
    {
      int minm = std::min(m + 1, n);

      TNT::Array2D<T> A(m, minm);

      for (int i = 0; i < m; i++)
        for (int j = 0; j < minm; j++)
          A[i][j] = U[i][j];

      return A;
    };

    /**
 * @brief Return the right singular vectors
 * 
 * @param V 
 */
    TNT::Array2D<T> getV()
    {
      return V.copy();
    }

    /**
 * @brief Return the one-dimensional array of singular values
 * 
 * @param x  diagonal of S.
 */
    TNT::Array1D<T> getSingularValues()
    {
      return s.copy();
    }

    /**
   * @brief Return the diagonal matrix of singular values
   * 
   * @param A S
   */
    TNT::Array2D<T> getS()
    {
       TNT::Array2D<T> A(n, n);
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          A[i][j] = 0.0;
        }
        A[i][i] = s[i];
      }
      return A;
    }

    /**
 * @brief  Two norm  
 * 
 * @return T (max(S))
 */
    T norm2()
    {
      return s[0];
    }

    /**
 * @brief Two norm condition number
 * 
 * @return T max(S)/min(S)
 */
    T cond()
    {
      return s[0] / s[std::min(m, n) - 1];
    }

    /**
 * @brief Effective numerical matrix rank
 * 
 * @return int  Number of nonnegligible singular values.
 */
    int rank()
    {
      T eps = std::pow(2.0, -52.0);
      T tol = std::max(m, n) * s[0] * eps;
      int r = 0;
      for (int i = 0; i < s.dim(); i++)
      {
        if (s[i] > tol)
        {
          r++;
        }
      }
      return r;
    }
  };

}
#endif
// JAMA_SVD_H
