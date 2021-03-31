#include "A2D_Struct.hpp"
/// @defgroup all_tests
/// @defgroup Array2D_tests
/// @ingroup all_tests

/**
* @brief Tests for Array 1D basic operations
* @addtogroup Basic_Operations 
* @ingroup Array2D_tests
* @{
*/

/**
* @brief A simple unit test to test size initialization 
* (rows)
*/
TEST_F(Array2D_test, dim1)
{
    TNT::Array2D<double> A;
    TNT::Array2D<double> B(m_, n_);
    TNT::Array2D<double> C(m_, n_, 0.);
    TNT::Array2D<double> D = TNT::Array2D<double>(m_, n_);

    EXPECT_EQ(A.dim1(), 0)
        << "A Size is false using dim1";
    EXPECT_EQ(B.dim1(), m_)
        << "B Size is false using dim1";
    EXPECT_EQ(C.dim1(), m_)
        << "C Size is false using dim1";
    EXPECT_EQ(D.dim1(), m_)
        << "D Size is false using dim1";
}

/**
* @brief A simple unit test to test size initialization 
* (columns)
*/
TEST_F(Array2D_test, dim2)
{
    TNT::Array2D<double> A;
    TNT::Array2D<double> B(m_, n_);
    TNT::Array2D<double> C(m_, n_, 0.);
    TNT::Array2D<double> D = TNT::Array2D<double>(m_, n_);

    EXPECT_EQ(A.dim2(), 0)
        << "A Size is false using dim2";
    EXPECT_EQ(B.dim2(), n_)
        << "B Size is false using dim2";
    EXPECT_EQ(C.dim2(), n_)
        << "C Size is false using dim2";
    EXPECT_EQ(D.dim2(), n_)
        << "D Size is false using dim2";
}

/**
* @brief Testing initialization with constant value
*/
TEST_F(Array2D_test, init_constant_values)
{
    // Using constant values
    TNT::Array2D<double> A(m_, n_, init_value_);

    for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
        {
            EXPECT_EQ(A[i][j], init_value_)
                << "A is not equal to init value []";
        }
}

/**
* @brief Testing initialization by assignation
*/
TEST_F(Array2D_test, init_assignation_values)
{
    // Using constant values
    TNT::Array2D<double> A(m_, n_);

    A = 1.0;

    for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
        {
            EXPECT_EQ(A[i][j], init_value_)
                << "A is not equal to init value []";
        }
}

/**
* @brief Testing initialization using a preexisting 
* bi-dimensional array
*/
TEST_F(Array2D_test, init_raw_array)
{
    // Using a pre existing array
    double ad[3][2] = {{0, 1}, {2, 3}, {4, 5}};
    TNT::Array2D<double> C(3, 2, *ad);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            EXPECT_EQ(C[i][j], ad[i][j])
                << "C is not equal to init values [][]";

    EXPECT_TRUE(&ad[0][0] == &C[0][0])
        << "C point to the same address of raw array";
}
/**
* @}
*/