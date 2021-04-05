#include "Jama_Struct.hpp"
/// @defgroup all_tests
/// @defgroup Jama_tests
/// @ingroup all_tests

/**
* @brief Tests for Jama class 
* @addtogroup Jama_Cholesky 
* @ingroup Jama_tests
* @{
*/

/**
 * @brief Test Cholesky factorization using a 
 * 3x3 case
 * 
 */
TEST_F(Jama_test, cholesky_3x3)
{
  // Matrix test
  double a_values[3][3] = {{4., 1., 1.},
                           {1., 2., 3.},
                           {1., 3., 6.}};
  TNT::Array2D<double> A(3, 3, *a_values);

  // Cholesky factorization
  JAMA::Cholesky<double> B(A);

  // Get the lower matrix
  TNT::Array2D<double> L = B.getL();

  // Verify the result 
  EXPECT_TRUE(A == TNT::mult_transpose(L, L))
    << "A is not equal to L*transpose(L)";
}

/**
 * @brief Test matrix inversion using
 * Cholesky factorization using a 3x3 case
 * 
 */
TEST_F(Jama_test, cholesky_inverse)
{
  // Matrix test
  double a_values[3][3] = {{4., 1., 1.},
                           {1., 2., 3.},
                           {1., 3., 6.}};
  TNT::Array2D<double> A(3, 3, *a_values);

  // Cholesky factorization
  JAMA::Cholesky<double> B(A);

  // Solve the system A*X= Identity
  TNT::Array2D<double> A_inv = B.solve(TNT::Tools<double>::eye(3));

  // X must be the inverse
  EXPECT_TRUE(TNT::near(A_inv,TNT::invert(A)))
      << "A_inv using cholesky is not correct";
}

/**
 * @brief Hilbert matrix example 
 */
TEST_F(Jama_test, hilbert_chol)
{
  /// H holds the n-by-n Hilbert matrix 
  TNT::Array2D<double> H(n_, n_);

  // Fill in the entries in H
  for (int i = 0; i < n_; i++)
  {
    for (int j = 0; j < n_; j++)
    {
      H[i][j] = 1. / (i + j + 1);
    }
  }

  // Use the Cholesky class to find a matrix C so that C*Ct = H
  JAMA::Cholesky<double> Hchol(H);
  TNT::Array2D<double> C = Hchol.getL();

  // Construct C transpose
  TNT::Array2D<double> Ctrans(n_, n_);
  for (int i = 0; i < n_; i++)
  {
    for (int j = 0; j < n_; j++)
    {
      Ctrans[i][j] = C[j][i];
    }
  }

  // Check that C*Ct gives H
  EXPECT_TRUE(H == TNT::matmult(C, Ctrans))
      << "A_inv using cholesky is not correct";
}

/**
* @}
*/