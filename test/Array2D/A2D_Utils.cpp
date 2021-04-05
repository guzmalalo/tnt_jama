#include "A2D_Struct.hpp"

/**
* @ingroup Array2D_tests
* @{
*/

/**
* @brief Testing the stream extraction operator >>
*/
TEST_F(Array2D_test, extraction_operator)
{
  TNT::Array2D<double> A;
  std::istringstream is("3 3 0 1 2 1 2 3 2 3 4");
  is >> A;

  EXPECT_EQ(A.dim1(), 3);
  EXPECT_EQ(A.dim2(), 3);

  // Verify the address

  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(A[i][j], i + j)
          << "A is not equal to init value []";
    }
}

/**
* @brief Testing the addition operation
*/
TEST_F(Array2D_test, addition_operator)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, init_value_);
  TNT::Array2D<double> C = A + B;

  // Verify the sum
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(C[i][j], (i + j) + init_value_)
          << "C is not equal to A+B";
    }
}

/**
* @brief Testing the substraction operation
*/
TEST_F(Array2D_test, subtraction_operator)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, init_value_);
  TNT::Array2D<double> C = A - B;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(C[i][j], (i + j) - init_value_)
          << "C is not equal to A-B";
    }
}

/**
* @brief Testing the multiplication operation
* element by element.
*/
TEST_F(Array2D_test, multiplication_operator)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, scalar_);
  TNT::Array2D<double> C = A * B;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(C[i][j], (i + j) * scalar_)
          << "C is not equal to A*B";
    }
}

/**
* @brief Testing the division operation
*/
TEST_F(Array2D_test, division_operator)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, scalar_);
  TNT::Array2D<double> C = A / B;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(C[i][j], (i + j) / scalar_)
          << "C is not equal to A/B";
    }
}

/**
* @brief Testing the addition assignement
*/
TEST_F(Array2D_test, addition_assignement)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, scalar_);

  A += B;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(A[i][j], (i + j) + scalar_)
          << "A is not equal to A+B";
    }
}

/**
* @brief Testing the substraction assignement
*/
TEST_F(Array2D_test, substraction_assignement)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, scalar_);

  A -= B;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(A[i][j], (i + j) - scalar_)
          << "A is not equal to A-B";
    }
}

/**
* @brief Testing the multiplication assignement
*/
TEST_F(Array2D_test, multiplication_assignement)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, scalar_);

  A *= B;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(A[i][j], (i + j) * scalar_)
          << "A is not equal to A*B";
    }
}

/**
* @brief Testing the multiplication assignement scalar
*/
TEST_F(Array2D_test, multiplication_assignement_scalar)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, *a_values);
  TNT::Array2D<double> C(3, 3, *a_values);

  // Addition assignement
  B = A * scalar_;
  C = scalar_ * A;
  A *= scalar_;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(A[i][j], (i + j) * scalar_)
          << "A is not equal to A*b";
      EXPECT_EQ(A[i][j], (i + j) * scalar_)
          << "B is not equal to A*b";
      EXPECT_EQ(A[i][j], (i + j) * scalar_)
          << "C is not equal to A*b";
    }
}

/**
* @brief Testing the division assignement 
*/
TEST_F(Array2D_test, division_assignement)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, scalar_);

  A /= B;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(A[i][j], (i + j) / scalar_)
          << "A is not equal to A/B";
    }
}

/**
* @brief Testing the division assignement scalar
*/
TEST_F(Array2D_test, division_assignement_scalar)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B(3, 3, *a_values);

  // Addition assignement
  B = A / scalar_;
  A /= scalar_;

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(A[i][j], (i + j) / scalar_)
          << "A is not equal to A/b";
      EXPECT_EQ(B[i][j], (i + j) / scalar_)
          << "B is not equal to A/b";
    }
}

/**
* @brief Testing the comparaison operator
*/
TEST_F(Array2D_test, comparaison)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);

  double b_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> B(3, 3, *b_values);

  EXPECT_TRUE(A == B)
      << "A is not equal to B";
}

/**
* @}
*/