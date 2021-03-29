#include "A1D_Struct.hpp"

/**
* @ingroup Array1D_tests
* @{
*/

/**
* @brief Testing the stream extraction operator >>
*/
TEST_F(Array1D_test, extraction_operator)
{
    TNT::Array1D<double> A(a_size_);
    std::istringstream is("5 1 2 3 4 5");
    is >> A;

    EXPECT_EQ(A.dim(), 5);

    // Verify the adress
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(A[i], i + 1)
            << "A is not initialized from istream";
    }
}

/**
* @brief Testing the addition operation
*/
TEST_F(Array1D_test, addition_operator)
{
    double a_values[5] = {1, 2, 3, 4, 5};
    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, scalar_);
    TNT::Array1D<double> C = A + B;

    ASSERT_TRUE(C.dim() > 0) << "C is null";

    // Verify the address
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(C[i], (i + 1) + scalar_)
            << "C is not equal to A+B";
    }
}

/**
* @brief Testing the substraction operation
*/
TEST_F(Array1D_test, subtraction_operator)
{
    double a_values[5] = {1, 2, 3, 4, 5};

    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, scalar_);
    TNT::Array1D<double> C = A - B;

    // Verify the adress
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(C[i], (i + 1) - scalar_)
            << "C is not equal to A+B";
    }
}

/**
* @brief Testing the multiplication operation
*/
TEST_F(Array1D_test, multiplication_operator)
{
    double a_values[5] = {1, 2, 3, 4, 5};

    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, scalar_);
    TNT::Array1D<double> C = A * B;

    // Verify the adress
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(C[i], (i + 1.0) * scalar_)
            << "C is not equal to A*B";
    }
}

/**
* @brief Testing the division operation
*/
TEST_F(Array1D_test, division_operator)
{
    double a_values[5] = {1., 2., 3., 4., 5.};

    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, scalar_);
    TNT::Array1D<double> C = A / B;

    // Verify the adress
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(C[i], (i + 1.0) / scalar_)
            << "C is not equal to A*B";
    }
}

/**
* @brief Testing the addition assignement
*/
TEST_F(Array1D_test, addition_assignement)
{
    double a_values[5] = {1, 2, 3, 4, 5};

    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, scalar_);

    // Addition assignement
    A += B;

    // Verify the adress
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(A[i], (i + 1) + scalar_)
            << "A+=B is not equal to A = A + B";
    }
}

/**
* @brief Testing the substraction assignement
*/
TEST_F(Array1D_test, substraction_assignement)
{
    double a_values[5] = {1, 2, 3, 4, 5};

    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, scalar_);

    // Addition assignement
    A -= B;

    // Verify the adress
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(A[i], (i + 1) - scalar_)
            << "A-=B is not equal to A = A - B";
    }
}

/**
* @brief Testing the multiplication assignement
*/
TEST_F(Array1D_test, multiplication_assignement)
{
    double a_values[5] = {1, 2, 3, 4, 5};
    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, scalar_);

    // Addition assignement
    A *= B;

    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(A[i], (i + 1) * scalar_)
            << "A is not equal to A*B";
    }
}

/**
* @brief Testing the multiplication assignement scalar
*/
TEST_F(Array1D_test, multiplication_assignement_scalar)
{
    double a_values[5] = {1, 2, 3, 4, 5};
    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, a_values);
    TNT::Array1D<double> C(5, a_values);

    // Addition assignement
    B = A * scalar_;
    C = scalar_ * A;
    A *= scalar_;

    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(A[i], (i + 1) * scalar_)
            << "A is not equal to A=*scalar";
        EXPECT_EQ(B[i], (i + 1) * scalar_)
            << "B is not equal to A*scalar";
        EXPECT_EQ(C[i], (i + 1) * scalar_)
            << "C is not equal to scalar*A";
    }
}

/**
* @brief Testing the division assignement 
*/
TEST_F(Array1D_test, division_assignement)
{
    double a_values[5] = {1., 2., 3., 4., 5.};

    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, scalar_);

    A /= B;

    // Verify the adress
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(A[i], (i + 1.0) / scalar_)
            << "A is not equal to A*B";
    }
}

/**
* @brief Testing the division assignement scalar
*/
TEST_F(Array1D_test, division_assignement_scalar)
{
    double a_values[5] = {1, 2, 3, 4, 5};
    TNT::Array1D<double> A(5, a_values);
    TNT::Array1D<double> B(5, a_values);
    TNT::Array1D<double> C(5, a_values);

    // Addition assignement
    B = A / scalar_;
    A /= scalar_;

    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(A[i], (i + 1) / scalar_)
            << "A is not equal to A=/scalar";
        EXPECT_EQ(B[i], (i + 1) / scalar_)
            << "B is not equal to A/scalar";
    }
}

/**
* @brief Testing the l2 norm of an Array1D 
*/
TEST_F(Array1D_test, l2_norm)
{
    double a_values[5] = {1.1, 2.1, 3.1, 4.1, 5.1};
    double res = 0.0;
    TNT::Array1D<double> A(5, a_values);
	for (int i=0; i<5; i++)
		res +=  a_values[i] * a_values[i];

    // Testing double values
    EXPECT_EQ(TNT::norm(A), std::sqrt(res));

    int b_values[5] = {1, 2, 3, 4, 5};
    TNT::Array1D<int> B(5, b_values);

    // Testing int values 
    EXPECT_EQ(TNT::norm(B), static_cast<int>(std::sqrt(55)));

}

/**
* @}
*/