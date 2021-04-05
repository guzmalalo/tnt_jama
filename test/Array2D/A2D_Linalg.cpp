#include "A2D_Struct.hpp"

/**
* @ingroup Array2D_tests
* @{
*/

/**
* @brief Testing the sum method
*/
TEST_F(Array2D_test, sum)
{
  double a_values[3][3] = {{2, 7, 6},
                           {9, 5, 1},
                           {4, 3, 8}};
  TNT::Array2D<double> A(3, 3, *a_values);

  TNT::Array1D<double> sum = TNT::sum(A);
  TNT::Array1D<double> sum_row = TNT::sum(A,1);
  TNT::Array1D<double> sum_col = TNT::sum(A,2);

  // Verify the result
  for (int j = 0; j < 3; j++)
  {
    EXPECT_EQ(sum_col[j], 15.0)
        << "Sum of magic matrix 3x3 is not 15";
    EXPECT_EQ(sum_col[j], 15.0)
        << "Sum of magic matrix col 3x3 is not 15";
    EXPECT_EQ(sum_row[j], 15.0)
        << "Sum of magic matrix row 3x3 is not 15";
  }
}

/**
* @brief Testing the norm-1 method 
*/
TEST_F(Array2D_test, norm_1)
{
  double a_values[3][3] = {{5, -4, 2},
                           {-1, 2, 3},
                           {-2, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  double norm_A = TNT::norm_1(A);

    EXPECT_EQ(norm_A, 8.0)
        << "Norm-1 for matrix A is not correct";
}

/**
* @brief Testing the infinity norm method for arrays2D
*/
TEST_F(Array2D_test, norm_inf)
{
  double a_values[3][3] = {{5, -4, 2},
                           {-1, 2, 3},
                           {-2, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  double norm_A = TNT::norm_inf(A);

  EXPECT_EQ(norm_A, 11.0)
      << "Infinity norm for matrix A is not correct";
}

/**
* @brief Testing the frobenius norm method for arrays2D
*/
TEST_F(Array2D_test, norm_fro)
{
  double a_values[3][3] = {{5, -4, 2},
                           {-1, 2, 3},
                           {-2, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  double norm_A = TNT::norm_fro(A);

  EXPECT_EQ(norm_A, 8.0)
      << "Infinity norm for matrix A is not correct";
}

/**
* @brief Testing the sum method
*/
TEST_F(Array2D_test, trace)
{
  double a_values[3][3] = {{2, 7, 6},
                           {9, 5, 1},
                           {4, 3, 8}};
  TNT::Array2D<double> A(3, 3, *a_values);

  double tr = TNT::trace(A);

  // Verify the result
    EXPECT_EQ(tr, 15.0)
        << "Trace of magic matrix 3x3 is not 15";

}


/**
* @brief Testing the matrix product
*/
TEST_F(Array2D_test, matrix_product)
{
  double a_values[5][5] = {
      {+0.061198, +0.201990, +0.019678, -0.493936, -0.126745},
      {+0.437242, +0.058956, -0.149362, -0.045465, +0.296153},
      {-0.492474, -0.031309, +0.314156, +0.419733, +0.068317},
      {+0.336352, +0.411541, +0.458476, -0.393139, -0.135040},
      {+0.239585, -0.428913, -0.406953, -0.291020, -0.353768}};
  TNT::Array2D<double> A(5, 5, *a_values);

  double b_values[5][5] = {
      {+0.051408, -0.126745, -0.493936, +0.019678, +0.201990},
      {+0.035437, +0.296153, -0.045465, -0.149362, +0.058956},
      {-0.454499, +0.068317, +0.419733, +0.314156, -0.031309},
      {+0.373833, -0.135040, -0.393139, +0.458476, +0.411541},
      {+0.258704, -0.353768, -0.291020, -0.406953, -0.428913}};
  TNT::Array2D<double> B(5, 5, *b_values);

  double res[5][5] = {
      {-0.2160787, +0.1649472, +0.1999190, -0.1976620, -0.1252585},
      {+0.1520715, -0.1467921, -0.3496545, -0.1884897, -0.0492639},
      {+0.0053737, -0.0062406, +0.1916407, +0.2583152, +0.0322787},
      {-0.3584056, +0.2114322, +0.2014480, -0.0361067, -0.0260243},
      {-0.0182371, -0.0207407, -0.0522859, -0.0485276, +0.0678171}};

  TNT::Array2D<double>
      C = TNT::matmult(A, B);

  // Verify the result
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
    {
      EXPECT_NEAR(C[i][j], res[i][j], 1e-5) << "C is not equal to A*B";
    }
}

/**
* @brief Testing the transpose method
*/
TEST_F(Array2D_test, transpose)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);
  TNT::Array2D<double> B = TNT::transpose(A);

  // Verify the result
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++)
    {
      EXPECT_EQ(A[i][j], B[j][i])
          << "B is not equal At";
    }
}

/**
* @brief Testing the mult_tranpose method
*/
TEST_F(Array2D_test, mult_transpose)
{
  double a_values[3][3] = {{1, -2, 2},
                           {2, -1, -2},
                           {2, 2, 1}};

  // An orthogonal matrix
  TNT::Array2D<double> A(3, 3, *a_values);
  A *= 1.0 / 3.0;

  // Long version
  TNT::Array2D<double> B = TNT::matmult(A, TNT::transpose(A));

  // Optimized
  TNT::Array2D<double> C = TNT::mult_transpose(A, A);

  for (int i = 0; i < m_; i++)
    EXPECT_DOUBLE_EQ(B[i][i], 1.0)
        << "B is not equal to identity";

  EXPECT_TRUE(B == C)
      << "B is not equal C";
}

/**
* @brief Testing the tranpose_mult method
*/
TEST_F(Array2D_test, transpose_mult)
{
  double a_values[3][3] = {{1, -2, 2},
                           {2, -1, -2},
                           {2, 2, 1}};

  // An orthogonal matrix
  TNT::Array2D<double> A(3, 3, *a_values);
  A *= 1.0 / 3.0;

  // Long version
  TNT::Array2D<double> B = TNT::matmult(TNT::transpose(A), A);

  // Optimized
  TNT::Array2D<double> C = TNT::transpose_mult(A, A);

  for (int i = 0; i < m_; i++)
    EXPECT_DOUBLE_EQ(B[i][i], 1.0)
        << "B is not equal to identity";

  // Verify the result
  EXPECT_TRUE(B == C)
      << "B is not equal C";
}

/**
* @brief Testing the Array2D-Array1D multiplication
*/
TEST_F(Array2D_test, matrix_vector_mult)
{
  double a_values[2][3] = {{1, -1, 2},
                           {0, -3, 1}};
  TNT::Array2D<double> A(2, 3, *a_values);

  double b_values[3] = {2, 1, 0};

  TNT::Array1D<double> b(3, b_values);

  TNT::Array1D<double> c = TNT::matmult(A, b);

  EXPECT_EQ(c.dim(), 2);
  EXPECT_EQ(c(1), 1);
  EXPECT_EQ(c(2), -3);
}

/**
* @brief Testing the transpose Array2D-Array1D 
* multiplication
*/
TEST_F(Array2D_test, matrix_vector_transmult)
{
  double a_values[3][2] = {{+1, +0},
                           {-1, -3},
                           {+2, +1}};
  TNT::Array2D<double> A(3, 2, *a_values);

  double b_values[3] = {2, 1, 0};
  TNT::Array1D<double> b(3, b_values);

  TNT::Array1D<double> c = TNT::transpose_mult(A, b);

  EXPECT_EQ(c.dim(), 2);
  EXPECT_EQ(c(1), 1);
  EXPECT_EQ(c(2), -3);
}
/**
* @brief Testing the inverse matrix method
*/
TEST_F(Array2D_test, inverse_matrix_2x2)
{
  double a_values[2][2] = {{1, 2},
                           {3, 4}};
  TNT::Array2D<double> A(2, 2, *a_values);

  TNT::Array2D<double> B = TNT::invert(A);

  double res[2][2] = {{4, -2},
                      {-3, 1}};
  TNT::Array2D<double> RES(2, 2, *res);
  RES *= -0.5;

  // Verify the result
  EXPECT_TRUE(B == RES)
      << "B is not equal inv A";
}
/**
* @brief Testing the inverse matrix method
*/
TEST_F(Array2D_test, inverse_matrix_3x3)
{
  double a_values[3][3] = {{1, 1, 1},
                           {1, 2, 3},
                           {0, 1, 0}};
  TNT::Array2D<double> A(3, 3, *a_values);

  TNT::Array2D<double> B = TNT::invert(A);

  double res[3][3] = {{-3, 1, 1},
                      {0, 0, -2},
                      {1, -1, 1}};
  TNT::Array2D<double> RES(3, 3, *res);
  RES *= -0.5;

  // Verify the result
  EXPECT_TRUE(B == RES)
      << "B is not equal inv A";
}

/**
* @brief Testing the inverse matrix method 
* with an orthogonal matrix
*/
TEST_F(Array2D_test, inverse_matrix_ortho)
{
  double a_values[3][3] = {{1, -2, 2},
                           {2, -1, -2},
                           {2, 2, 1}};

  // An orthogonal matrix
  TNT::Array2D<double> A(3, 3, *a_values);
  A *= 1.0 / 3.0;

  // Inverted matrix
  TNT::Array2D<double> B = TNT::invert(A);

  // transposed matrix
  TNT::Array2D<double> C = TNT::transpose(A);

  // Verify the result
  EXPECT_TRUE(B == C)
      << "inv A is not equal trans A";
}
/**
* @}
*/