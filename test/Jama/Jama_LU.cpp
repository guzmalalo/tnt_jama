#include "Jama_Struct.hpp"
/// @defgroup all_tests
/// @defgroup Jama_tests
/// @ingroup all_tests

/**
* @brief Tests for Jama class 
* @addtogroup Jama_LU 
* @ingroup Jama_tests
* @{
*/

/**
* @brief Test LU factorization 2x2
*/
TEST_F(Jama_test, LU_2x2)
{
  double a_values[2][2] = {{4, 3},
                           {6, 3}};
  TNT::Array2D<double> A(2, 2, *a_values);

  JAMA::LU<double> B(A);

  // lower
  TNT::Array2D<double> L(B.getL());
  TNT::Array2D<double> U(B.getU());
  TNT::Array1D<int> P(B.getPivot());

  // Verify the result
  EXPECT_TRUE(A.subarray(P,0,1) == TNT::matmult(L,U))
          << "A is not equal to L*U";
}

/**
* @brief Test LU det
*/
TEST_F(Jama_test, LU_det)
{
  double a_values[3][3] = {{1, 1, 1},
                           {1, 2, 3},
                           {0, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  JAMA::LU<double> B(A);

  // Verify the result
  EXPECT_EQ(B.det(),TNT::det(A))
      << "det A LU is not equal det(A)";
}

/**
* @brief Test LU inverse
*/
TEST_F(Jama_test, LU_inverse)
{
  double a_values[3][3] = {{1, 1, 1},
                           {1, 2, 3},
                           {0, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  // LU object
  JAMA::LU<double> B(A);

  // Set up identity matrix
  TNT::Array2D<double> eye(3, 3, 0.);
  for (int i = 0; i < 3; i++)
    eye[i][i] = 1;

  // Compute inverse using LU
  TNT::Array2D<double> A_inv = B.solve(eye);

  // Verify the result
  EXPECT_TRUE(A_inv == TNT::invert(A))
      << "A inv is not equal invert(A)";
}

/**
* @brief Test LU solve for Ax=b
*/
TEST_F(Jama_test, LU_solve)
{
  double a_values[3][3] = {{7, -2, 1},
                           {14, -7, -3},
                           {-7, 11, 18}};
  TNT::Array2D<double> M(3, 3, *a_values);

  // LU object
  JAMA::LU<double> A(M);

  // Set up identity matrix
  TNT::Array1D<double> b(3);
  b[0] = 12.;
  b[1] = 17.;
  b[2] = 5.;

  // Compute Ax = b, to find x
  TNT::Array1D<double> x = A.solve(b);

  // Verify the result
  EXPECT_NEAR(x[0],3.,1e-7)
      << "x is not equal ref";
  EXPECT_NEAR(x[1], 4.,1e-7)
      << "x is not equal ref";
  EXPECT_NEAR(x[2], -1.,1e-7)
      << "x is not equal ref";
}

/**
* @brief Test LU inverse of hilbert matrix
*/
TEST_F(Jama_test, LU_hilbert)
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
  // Set up identity matrix
  TNT::Array2D<double> eye(n_, n_, 0.);
  for (int i = 0; i < n_; i++)
    eye[i][i] = 1;

  JAMA::LU<double> HLU(H);
  TNT::Array2D<double> Hinv = HLU.solve(eye);

  // Verify the result
  EXPECT_TRUE(TNT::near(Hinv,TNT::invert(H)))
      << "H inv is not equal invert(H)";

  // Verify the result
  EXPECT_NEAR(HLU.det(), TNT::det(H),1e-12)
      << "H inv is not equal invert(H)";
}

/**
* @}
*/