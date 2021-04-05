#include "Jama_Struct.hpp"
/// @defgroup all_tests
/// @defgroup Jama_tests
/// @ingroup all_tests

/**
* @brief Tests for Jama class 
* @addtogroup Jama_Eigen 
* @ingroup Jama_tests
* @{
*/

/**
* @brief Test Eigen symmetric positive matrix
*/
TEST_F(Jama_test, eigen_sym_pos)
{
  // epsilon standards double 
  double eps = std::pow(2.0, -53.0);
  // Symmetric positive Matrix test
  double a_values[3][3] = {{4., 1., 1.},
                           {1., 2., 3.},
                           {1., 3., 6.}};
  TNT::Array2D<double> A(3, 3, *a_values);
  JAMA::Eigenvalue<double> B(A);

  TNT::Array2D<double> D;
  B.getD(D);

  TNT::Array2D<double> V;
  B.getV(V);

  TNT::Array2D<double> AV = TNT::matmult(A, V);
  TNT::Array2D<double> VD = TNT::matmult(V, D);

  // Verify the result
  EXPECT_TRUE(TNT::near(AV, VD));
}

/**
* @brief Test Eigen 2x2
*/
TEST_F(Jama_test, eigen_2x2)
{
  double a_values[2][2] = {{0, 1},
                           {-2, -3}};
  TNT::Array2D<double> A(2, 2, *a_values);

  JAMA::Eigenvalue<double> B(A);

  TNT::Array2D<double> D;
  B.getD(D);

  TNT::Array2D<double> V;
  B.getV(V);

  TNT::Array2D<double> AV = TNT::matmult(A, V);
  TNT::Array2D<double> VD = TNT::matmult(V, D);
  // Verify the result
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
    {
      EXPECT_DOUBLE_EQ(AV[i][j], VD[i][j])
          << "AV is not equal to VD";
    }
}

/**
 * @brief Test using Magic squares,
 * max_eig = maximum eigenvalue of (A + A')/2, should equal
 * trace.
 * 
 */
TEST_F(Jama_test, max_eig)
{
  // Create a magic matrix of dimension 5
  TNT::Array2D<double> M = TNT::Tools<double>::magic(5);

  // Create an eigenvalue object
  JAMA::Eigenvalue<double>
      B((M + TNT::transpose(M)) * 0.5);

  // Get the real eigenvalues 
  TNT::Array1D<double> d;
  B.getTEigenvalues(d);

  // Find the max eigen value
  double max_value = TNT::max(d);

  // Compute the trace of the magi matrix
  double tr = TNT::trace(M);

  // Evaluate the result
  EXPECT_DOUBLE_EQ(max_value, tr)
      << "Max eigenvalue of (A + A')/2 != tr(A) equal";
}
/**
* @}
*/