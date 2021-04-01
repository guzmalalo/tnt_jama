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

    double a = C[2][1];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            EXPECT_EQ(C[i][j], ad[i][j])
                << "C is not equal to init values [][]";

    EXPECT_TRUE(&ad[0][0] == &C[0][0])
        << "C point to the same address of raw array";

    // If the Array2d is destroyed the raw array is preserved
    // conserved
}

/**
* @brief Testing the copy constructor
*/
TEST_F(Array2D_test, copy_constructor)
{
    // Reference array
    TNT::Array2D<double> A(m_, n_, init_value_);
    // Copy constructor (a reference to A)
    TNT::Array2D<double> B(A);

    for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
            EXPECT_EQ(A[i][j], B[i][j])
                << "A is not equal B [][]";

    // Checking the number of references to A
    EXPECT_EQ(A.ref_count(), 2)
        << "A has only one reference instead of two";
}

/**
* @brief Testing the  operator []
*/
TEST_F(Array2D_test, bracket_operator)
{
    // Reference array
    TNT::Array2D<double> A(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A[i][j] = i + j;

    // Extracting the rows
    TNT::Array1D<double> row_0(n_, A[0]);
    TNT::Array1D<double> row_1(n_, A[1]);
    TNT::Array1D<double> row_2(n_, A[2]);

    // or using raw arrays
    double *p_row0 = A[0];
    double *p_row1 = A[1];
    double *p_row2 = A[2];

    // Comparing the values
    for (int i = 0; i < 3; i++)
    {
        EXPECT_EQ(row_0[i], A[0][i])
            << "Row 0 is not A[0]";
        EXPECT_EQ(row_1[i], A[1][i])
            << "Row 1 is not A[0]";
        EXPECT_EQ(row_2[i], A[2][i])
            << "Row 2 is not A[0]";

        EXPECT_EQ(p_row0[i], A[0][i])
            << "PRow 0 is not A[0]";
        EXPECT_EQ(p_row1[i], A[1][i])
            << "PRow 1 is not A[0]";
        EXPECT_EQ(p_row2[i], A[2][i])
            << "PRow 2 is not A[0]";
    }

    // cheking the pointers; both must point to the first  element
    // of the i row
    EXPECT_EQ(&row_0[0], A[0]) << "Row 0 is not A[0]";
}

/**
* @brief Testing the  operator [] const
*/
TEST_F(Array2D_test, bracket_operator_const)
{
    // Reference array
    TNT::Array2D<double> A(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A[i][j] = i + j;

    // Extracting the rows
    const TNT::Array1D<double> row_0(n_, A[0]);
    const TNT::Array1D<double> row_1(n_, A[1]);
    const TNT::Array1D<double> row_2(n_, A[2]);

    // or using raw arrays
    const double *p_row0 = A[0];
    const double *p_row1 = A[1];
    const double *p_row2 = A[2];

    // Comparing the values
    for (int i = 0; i < 3; i++)
    {
        EXPECT_EQ(row_0[i], A[0][i])
            << "Row 0 is not A[0]";
        EXPECT_EQ(row_1[i], A[1][i])
            << "Row 1 is not A[0]";
        EXPECT_EQ(row_2[i], A[2][i])
            << "Row 2 is not A[0]";

        EXPECT_EQ(p_row0[i], A[0][i])
            << "PRow 0 is not A[0]";
        EXPECT_EQ(p_row1[i], A[1][i])
            << "PRow 1 is not A[0]";
        EXPECT_EQ(p_row2[i], A[2][i])
            << "PRow 2 is not A[0]";
    }

    // cheking the pointers; both must point to the first  element
    // of the i row
    EXPECT_EQ(&row_0[0], A[0]) << "Row 0 is not A[0]";
}

/**
* @brief Testing the copy method
*/
TEST_F(Array2D_test, copy)
{
    // Reference array
    TNT::Array2D<double> A(m_, n_, init_value_);
    TNT::Array2D<double> B(A.copy());

    // Check if B is equal to A
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            EXPECT_EQ(A[i][j], B[i][j])
                << "B is not equal to A";

    // Check the pointed address for data
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            EXPECT_NE(&A[i][j], &B[i][j])
                << "B address equal to A";
}

/**
* @brief Test the inject method
*/
TEST_F(Array2D_test, inject_method)
{
    // Reference array
    TNT::Array2D<double> A(m_, n_, init_value_);
    TNT::Array2D<double> B(m_, n_);

    B.inject(A);

    // Check if B is equal to A
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            EXPECT_EQ(A[i][j], B[i][j])
                << "B is not equal to A";

    // Check the pointed address for data
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            EXPECT_NE(&A[i][j], &B[i][j])
                << "B address equal to A";
}

/**
* @brief Test the ref method
*/
TEST_F(Array2D_test, ref_method)
{
    // Definition des arrays
    TNT::Array2D<double> A(m_, n_, init_value_);
    TNT::Array2D<double> B(4, 4, init_value_);
    TNT::Array2D<double> C;
    TNT::Array2D<double> *D = new TNT::Array2D<double>(2, 2);

    // B and C must be now a reference to A.
    B.ref(A);
    C.ref(A);
    D->ref(A);

    // Check dimension
    EXPECT_EQ(A.dim1(), B.dim1());
    EXPECT_EQ(A.dim1(), C.dim1());
    EXPECT_EQ(A.dim1(), D->dim1());

    EXPECT_EQ(A.dim2(), B.dim2());
    EXPECT_EQ(A.dim2(), C.dim2());
    EXPECT_EQ(A.dim2(), D->dim2());

    // There must be 4 references to A
    EXPECT_EQ(A.ref_count(), 4);

    // Verify the values
    for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
        {
            EXPECT_EQ(A[i][j], B[i][j])
                << "B is not equal to A";
            EXPECT_EQ(A[i][j], C[i][j])
                << "C is not equal to A";
            EXPECT_EQ(A[i][j], (*D)[i][j])
                << "D is not equal to A";
        }

    // Delete a pointer
    delete D;
    D = nullptr;

    // There must be 3 references to A
    EXPECT_EQ(A.ref_count(), 3)
        << "Ref_count is not correct";
}

/**
* @brief Test the = operator
*/
TEST_F(Array2D_test, equal_operator)
{
    // Definition des arrays
    TNT::Array2D<double> A(m_, n_, init_value_);
    TNT::Array2D<double> Z(m_, n_, init_value_);
    TNT::Array2D<double> B(4, 4, init_value_);
    TNT::Array2D<double> C;
    TNT::Array2D<double> *D = new TNT::Array2D<double>(2, 2);

    // B and C must be now a reference to A.
    B = A;
    C = A;
    *D = A;

    // Check dimension
    EXPECT_EQ(A.dim1(), B.dim1());
    EXPECT_EQ(A.dim1(), C.dim1());
    EXPECT_EQ(A.dim1(), D->dim1());

    EXPECT_EQ(A.dim2(), B.dim2());
    EXPECT_EQ(A.dim2(), C.dim2());
    EXPECT_EQ(A.dim2(), D->dim2());

    // There must be 4 references to A
    EXPECT_EQ(A.ref_count(), 4);

    // Verify the values
    for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
        {
            EXPECT_EQ(A[i][j], B[i][j])
                << "B is not equal to A";
            EXPECT_EQ(A[i][j], C[i][j])
                << "C is not equal to A";
            EXPECT_EQ(A[i][j], (*D)[i][j])
                << "D is not equal to A";
        }

    // Delete a pointer
    delete D;
    D = nullptr;

    // There must be 3 references to A
    EXPECT_EQ(A.ref_count(), 3)
        << "Ref_count is not correct";

    // Change the reference of B to C
    B = C;

    // There must be 3 references to A
    EXPECT_EQ(A.ref_count(), 3)
        << "Ref_count is not correct";

    // Change the reference of C to Z
    C = Z;

    // There must be 2 references to A
    EXPECT_EQ(A.ref_count(), 2)
        << "Ref_count to A is not correct";

    // There must be 2 references to Z
    EXPECT_EQ(Z.ref_count(), 2)
        << "Ref_count to Z is not correct";
}

/**
* @brief Conversion of Array2D to a  regular 
         multidimensional C pointer.
*/
TEST_F(Array2D_test, conversion_operator)
{
    TNT::Array2D<double> A(m_, n_, init_value_);

    // Conversion of A to ad
    // Information is cloned
    double **ad = A;

    // Check if ad is equal to the init value
    for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
        {
            EXPECT_EQ(ad[i][j], init_value_)
                << "Raw array is not equal to init value";
            EXPECT_EQ(&ad[i][j], &A[i][j])
                << "Raw array is not equal to A";
        }

    EXPECT_EQ(A.ref_count(), 1) << "There are more than one reference ";
}

/**
* @brief Test the ref_count method
*/
TEST_F(Array2D_test, ref_count_method)
{

    // Definition des arrays
    TNT::Array2D<double> A(m_, n_, init_value_);
    TNT::Array2D<double> Z(m_, n_, init_value_);
    TNT::Array2D<double> B(4, 4, init_value_);
    TNT::Array2D<double> C;
    TNT::Array2D<double> *D = new TNT::Array2D<double>(2, 2);

    // B and C must be now a reference to A.
    B.ref(A);
    C.ref(A);
    D->ref(A);

    // There must be 4 references to A
    EXPECT_EQ(A.ref_count(), 4)
        << "Ref_count is not correct";

    // Delete a pointer
    delete D;
    D = nullptr;

    // There must be 3 references to A
    EXPECT_EQ(A.ref_count(), 3)
        << "Ref_count is not correct";

    // Change the reference of B to C
    B.ref(C);

    // There must be 3 references to A
    EXPECT_EQ(A.ref_count(), 3)
        << "Ref_count is not correct";

    // Change the reference of C to Z
    C.ref(Z);

    // There must be 2 references to A
    EXPECT_EQ(A.ref_count(), 2)
        << "Ref_count to A is not correct";

    // There must be 2 references to Z
    EXPECT_EQ(Z.ref_count(), 2)
        << "Ref_count to Z is not correct";
}

/**
* @brief Test the subarray method
*/
TEST_F(Array2D_test, subarray_method)
{
    TNT::Array2D<double> A(m_, n_, init_value_);

    // Fill the array
    for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
        {
            A[i][j] = i + j;
        }

    // Create an array from a subarray A
    TNT::Array2D<double> B = A.subarray(0, 1, 0, 1);

    EXPECT_EQ(B.dim1(), 2);

    // Verify the address

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            EXPECT_EQ(&A[i][j], &B[i][j])
                << "B is not equal to A";
        }

    // Invalid subarray
    TNT::Array2D<double> C = A.subarray(3, 0, 2, 0);
    EXPECT_EQ(C.dim1(), 0);

    // Invalid subarray
    TNT::Array2D<double> D = A.subarray(-1, 2, -1, 2);
    EXPECT_EQ(D.dim1(), 0);

    // Invalid subarray
    TNT::Array2D<double> E = A.subarray(0, 4, 0, 4);
    EXPECT_EQ(E.dim1(), 0);

    // Scalar array
    TNT::Array2D<double> F = A.subarray(1, 1, 1, 1);
    EXPECT_EQ(F.dim1(), 1);

    // Row array
    TNT::Array2D<double> G = A.subarray(0, 0, 0, 2);
    EXPECT_EQ(G.dim1(), 1);
    EXPECT_EQ(G.dim2(), 3);

    // Column array
    TNT::Array2D<double> H = A.subarray(0, 2, 0, 0);
    EXPECT_EQ(H.dim1(), 3);
    EXPECT_EQ(H.dim2(), 1);


 
}

/**
* @brief Test the subarray method with vector
*/
TEST_F(Array2D_test, subarray_method_rows)
{
  double a_values[3][3] = {{0, 1, 2},
                           {1, 2, 3},
                           {2, 3, 4}};
  TNT::Array2D<double> A(3, 3, *a_values);

  // Define the rows index to use in the new subarray
  TNT::Array1D<int> rows(2);
  rows[0] = 0;
  rows[1] = 2;

  // Extracting the values from the first and last line 
  // and from the second to the third colum 
  TNT::Array2D<double> B = A.subarray(rows, 1, 2);
  EXPECT_EQ(B.dim1(), 2);
  EXPECT_EQ(B.dim2(), 2);
  EXPECT_EQ(B[0][0], 1);
  EXPECT_EQ(B[0][1], 2);
  EXPECT_EQ(B[1][0], 3);
  EXPECT_EQ(B[1][1], 4);

  // Extracting all the values
  TNT::Array1D<int> rows2(3);
  rows2[0] = 0;
  rows2[1] = 1;
  rows2[2] = 2;

  B = A.subarray(rows2, 0, 2);

  EXPECT_TRUE(A == B);

  // Extract the last value of A
  // Extracting all the values
  TNT::Array1D<int> rows3(1);
  rows3[0] = 2;
  B = A.subarray(rows3, 2, 2);
  EXPECT_EQ(B.dim1(), 1);
  EXPECT_EQ(B.dim2(), 1);

  EXPECT_TRUE(A[2][2] == B[0][0]);
}

/**
* @}
*/