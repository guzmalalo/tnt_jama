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
    TNT::Array2D<double> B(3, 3, 5.5);
    TNT::Array2D<double> C = A * B; 

    // Verify the result
    for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
        {
            EXPECT_EQ(C[i][j], (i + j) *5.5)
                << "C is not equal to A*B";
        }
}
/**
* @}
*/