#include "A1D_Struct.hpp"

/**
* @ingroup Array1D_tests
* @{
*/

/**
* @brief Testing the sum method 
*/
TEST_F(Array1D_test, sum)
{
  int n = 100;
  double a_values[n];

  for (int i = 0; i < n; i++)
  {
    a_values[i] = i+1;
  }

  TNT::Array1D<double> A(n, a_values);

  // Testing double values
  EXPECT_EQ(TNT::sum(A), n*(n+1)/2.0);
}

/**
* @brief Testing the abs method 
*/
TEST_F(Array1D_test, abs)
{
  int n = 100;
  double a_values[n];

  for (int i = 0; i < n; i++)
  {
    a_values[i] = (i + 1.0)*std::pow(-1.,i);
  }

  TNT::Array1D<double> A(n, a_values);

  A = TNT::abs(A);
  // Testing double values
  EXPECT_EQ(TNT::sum(A), n * (n + 1) / 2.0);
}

/**
* @brief Testing the l1-norm of an Array1D 
*/
TEST_F(Array1D_test, l1_norm)
{
  int n = 100;
  double a_values[n];

  for (int i = 0; i < n; i++)
  {
    a_values[i] = (i + 1) * std::pow(-1., i);
  }

  TNT::Array1D<double> A(n, a_values);

  // Testing double values
  EXPECT_EQ(TNT::norm_1(A), n * (n + 1) / 2.0);
}

/**
* @brief Testing the l2 norm of an Array1D 
*/
TEST_F(Array1D_test, l2_norm)
{
  double a_values[5] = {1.1, 2.1, 3.1, 4.1, 5.1};
  double res = 0.0;
  TNT::Array1D<double> A(5, a_values);
  for (int i = 0; i < 5; i++)
    res += a_values[i] * a_values[i];

  // Testing double values
  EXPECT_EQ(TNT::norm(A), std::sqrt(res));

  int b_values[5] = {1, 2, 3, 4, 5};
  TNT::Array1D<int> B(5, b_values);

  // Testing int values
  EXPECT_EQ(TNT::norm(B), static_cast<int>(std::sqrt(55)));
}

/**
* @brief Testing the lp norm of an Array1D 
*/
TEST_F(Array1D_test, lp_norm)
{
  double a_values[5] = {1.1, 2.1, 3.1, 4.1, 5.1};
  double res = 0.0;
  TNT::Array1D<double> A(5, a_values);
  for (int i = 0; i < 5; i++)
    res += std::pow(a_values[i],5);

  // Testing double values
  EXPECT_EQ(TNT::norm_p(A,5), std::pow(res,1.0/5.0));

}

/**
* @brief Testing the cross product of two Array1D 
*/
TEST_F(Array1D_test, cross_product)
{
  double a_values[3] = {0.1, 2.3, 4.5};
  double b_values[3] = {6.7, 8.9, 10.0};
  TNT::Array1D<double> A(3, a_values);
  TNT::Array1D<double> B(3, b_values);

  TNT::Array1D<double> C = TNT::cross(A, B);

  // Testing int values
  EXPECT_DOUBLE_EQ(C(1), -17.050);
  EXPECT_DOUBLE_EQ(C(2), 29.150);
  EXPECT_DOUBLE_EQ(C(3), -14.520);
}

/**
* @brief Testing the dot product of two Array1D 
*/
TEST_F(Array1D_test, dot_product)
{
  double a_values[3] = {-2, 3, 4};
  double b_values[3] = {4, 6, 7};
  TNT::Array1D<double> A(3, a_values);
  TNT::Array1D<double> B(3, b_values);

  double C = TNT::dot_product(A, B);

  // Testing int values
  EXPECT_EQ(C, -8 + 18 + 28);
}

/**
* @brief Testing the max element in array
*/
TEST_F(Array1D_test, max)
{
  double a_values[6] = {-2, 0, 3, 4,+100,1e5};
  TNT::Array1D<double> A(6, a_values);
  double C = TNT::max(A);

  // Testing int values
  EXPECT_EQ(C, 1e5);
}

/**
* @brief Testing the max element in array
*/
TEST_F(Array1D_test, min)
{
  double a_values[6] = {-2, 0, 3, 4, +100, 1e5};
  TNT::Array1D<double> A(6, a_values);
  double C = TNT::min(A);

  // Testing int values
  EXPECT_EQ(C, -2);
}

/**
* @}
*/