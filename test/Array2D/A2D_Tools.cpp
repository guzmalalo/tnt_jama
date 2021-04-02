#include "A2D_Struct.hpp"

/**
* @ingroup Array2D_tests
* @{
*/

/**
* @brief Testing the eye tool
*/
TEST_F(Array2D_test, eye)
{
  TNT::Array2D<double> A = TNT::Tools<double>::eye(n_);

  // Check the matrix
  for (int i = 0; i < n_; i++)
  {
    for (int j = 0; j < n_; j++)
    {
      if (i == j)
        EXPECT_EQ(A[i][j], 1.0);
      else
        EXPECT_EQ(A[i][j], 0.0);
    }
  }
}

/**
* @brief Testing the magic matrix
*/
TEST_F(Array2D_test, magic)
{
  TNT::Array2D<double> A = TNT::Tools<double>::magic(n_);

  double sum_row;
  double sum_col;
  double sum_dia = 0.;
  // Verify results rows
  for (int i = 0; i < n_; i++)
  {
    sum_row = 0.;
    for (int j = 0; j < n_ ; j++)
    {
      sum_row += A[i][j];
    }
    EXPECT_EQ(sum_row, (std::pow(n_, 3) + n_) / 2.0);
  }

  for (int j = 0; j < n_; j++)
  {
    sum_col = 0.;
    for (int i = 0; i < n_; i++)
    {
      sum_col += A[i][j];
    }
    EXPECT_EQ(sum_col, (std::pow(n_, 3) + n_) / 2.0);
  }

  // Verify results rows
  for (int i = 0; i < n_; i++)
  {
      sum_dia += A[i][i];
  }

  EXPECT_EQ(sum_dia, (std::pow(n_,3)+n_)/2.0);
}

/**
* @}
*/