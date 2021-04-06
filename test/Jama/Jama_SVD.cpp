#include "Jama_Struct.hpp"
/// @defgroup all_tests
/// @defgroup Jama_tests
/// @ingroup all_tests

/**
* @brief Tests for Jama class 
* @addtogroup Jama_SVD 
* @ingroup Jama_tests
* @{
*/

/**
* @brief Test SVD factorization 3x3
*/
TEST_F(Jama_test, SVD_3x3)
{
  double a_values[3][3] = {{1, 0, 1},
                           {-1, -2, 0},
                           {0, 1, -1}};
  TNT::Array2D<double> A(3, 3, *a_values);

  JAMA::SVD<double> B(A);

  // Singular values
  TNT::Array1D<double> s(B.getSingularValues());
 
  EXPECT_NEAR(s[0], 2.4605,1e-4);
  EXPECT_NEAR(s[1], 1.6996,1e-4);
  EXPECT_NEAR(s[2], 0.2391,1e-4);
}

/**
* @brief Test SVD factorization 3x3
*/
TEST_F(Jama_test, SVD_rectangular)
{
  double columnwise[12]= {1., 2., 3., 4., 5., 6., 7., 8.};
  TNT::Array2D<double> A(4, 2, columnwise);

  // SVD object
  JAMA::SVD<double> B(A);

  // Singular right vectors 
  TNT::Array2D<double> V(B.getV());

  // Singular left vectors
  TNT::Array2D<double> U(B.getU());

  // Singular values matrix
  TNT::Array2D<double> S(B.getS());

  // Reconstruction 
  TNT::Array2D<double> USV= TNT::matmult(U, TNT::mult_transpose(S, V));

  // Verify the result
  EXPECT_TRUE(TNT::near(A,USV));
}

/**
* @brief Test SVD rank mat
*/
TEST_F(Jama_test, SVD_rank)
{
  double a_values[4][4] = {{1, 0, 2, 3},
                           {2, 0, 4, 6},
                           {0, 2, 2, 0},
                           {1, 2, 4, 3}};
  TNT::Array2D<double> A(4, 4, *a_values);

  // SVD object
  JAMA::SVD<double> B(A);

  // Verify the result
  EXPECT_EQ(B.rank(),2);
}


/**
* @brief Test SVD cond
*/
TEST_F(Jama_test, SVD_cond)
{
  TNT::Array2D<double> H = TNT::Tools<double>::hilbert(3);

  // SVD object
  JAMA::SVD<double> B(H);

  EXPECT_NEAR(B.cond(), 524.0567775860,1e-7);
}

/**
* @brief Test SVD cond with jama test
*/
TEST_F(Jama_test, SVD_cond_jama)
{
  double a_values[3][3] = {{1., 3.},
                           {7., 9.}};

  TNT::Array2D<double> A(3, 3, *a_values);

  JAMA::SVD<double> B(A);

  // Singular values
  TNT::Array1D<double> s(B.getSingularValues());

  double c = s[0]/s[std::min(A.dim1(),A.dim2())-1];

  EXPECT_EQ(B.cond(),c);
}

/**
* @}
*/