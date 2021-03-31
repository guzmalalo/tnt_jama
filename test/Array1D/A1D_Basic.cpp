#include "A1D_Struct.hpp"
/// @defgroup all_tests
/// @defgroup Array1D_tests
/// @ingroup all_tests

/**
* @brief Tests for Array 1D basic operations
* @addtogroup Basic_Operations 
* @ingroup Array1D_tests
* @{
*/

/**
* @brief A simple unit test to test size initialization
*/
TEST_F(Array1D_test, dim1)
{
    TNT::Array1D<double> A;
    TNT::Array1D<double> B(a_size_);
    TNT::Array1D<double> C(a_size_, 0.);
    TNT::Array1D<double> D = TNT::Array1D<double>(a_size_);

    EXPECT_EQ(A.dim1(), 0)
        << "A Size is false using dim1";
    EXPECT_EQ(B.dim1(), a_size_)
        << "B Size is false using dim1";
    EXPECT_EQ(C.dim1(), a_size_)
        << "C Size is false using dim1";
    EXPECT_EQ(D.dim1(), a_size_)
        << "D Size is false using dim1";
}

/**
* @brief A simple unit test to test values initialization
*/
TEST_F(Array1D_test, init_values)
{
    // An array of doubles
    double ad[a_size_];
    for (int i = 0; i < a_size_; ++i)
    {
        ad[i] = i;
    }

    // Initialize with a constat value
    TNT::Array1D<double> A(a_size_, init_value_);
    // Initialize with out a value
    TNT::Array1D<double> B(a_size_);
    // Initialize using an array as an entry
    TNT::Array1D<double> C(a_size_, ad);

    B = 1.0;

    for (int i = 0; i < a_size_; i++)
    {
        EXPECT_EQ(A[i], init_value_)
            << "A is not equal to init value []";
        EXPECT_EQ(A[i], B[i])
            << "A is not equal to B []";
        EXPECT_EQ(C[i], ad[i])
            << "C is not equal to raw array []";
    }
    for (int i = 1; i <= A.dim(); i++)
    {
        EXPECT_EQ(A(i), init_value_)
            << "A is not equal to init value ()";
        EXPECT_EQ(A(i), B(i))
            << "A is not equal to B ()";
        EXPECT_EQ(C(i), ad[i - 1])
            << "C is not equal to raw array ()";
    }
}

/**
* @brief A simple unit test for clone operations
*/
TEST_F(Array1D_test, clone)
{
    TNT::Array1D<double> A(a_size_, init_value_);
    TNT::Array1D<double> B(a_size_);

    // Clone using constructor
    TNT::Array1D<double> C(A);

    // Clone using = operator
    B = A;

    // A and B are different objets pointing to the same
    // data
    EXPECT_FALSE(&A == &B)
        << "B has the same address of A ";

    // B and C are different objets pointing to the same
    // data
    EXPECT_FALSE(&B == &C)
        << "B has the same address of C ";

    // Checking if A, B anc C point to the same data
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(&A[i], &B[i])
            << "B is not a clone of A []";
        EXPECT_EQ(&C[i], &A[i])
            << "C is not a clone of A []";
        EXPECT_EQ(&C[i], &B[i])
            << "C is not a clone of B []";
    }

    // Check if C has been initialized to 1.
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(C[i], init_value_)
            << "C is not equal to init value";
    }

    // Modify all C values to 0.0
    C = 0.0;

    // Check if A and B has been modified to 0.
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(A[i], 0.0)
            << "A is not zero";
        EXPECT_EQ(B[i], 0.0)
            << "B is not zero";
    }
}

/**
* @brief A simple unit test for clone operations
    using an array
*/
TEST_F(Array1D_test, clone_raw_array)
{
    // An array of doubles
    double ad[a_size_];
    for (int i = 0; i < a_size_; ++i)
    {
        ad[i] = i;
    }
    { // Scope of Array1D
        // Initialize with a constat value
        TNT::Array1D<double> A(a_size_, ad);

        // Checking if A, B anc C point to the same data
        for (int i = 0; i < a_size_; i++)
        {
            EXPECT_EQ(&A[i], &ad[i])
                << "A is not a clone of ad []";
        }
    } // Array1D is destroyed

    // Check if the raw array is still here
    for (int i = 0; i < a_size_; i++)
    {
        EXPECT_EQ(ad[i], i)
            << "ad has been destroyed ";
    }
}

/**
* @brief Conversion of Array1D to a raw array 
*/
TEST_F(Array1D_test, conversion_operator)
{
    TNT::Array1D<double> A(a_size_, init_value_);

    // Conversion of A to ad
    // Information is cloned
    double *ad = A;

    // Check if ad is equal to the init value
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(ad[i], init_value_)
            << "Raw array is not equal to A";
        EXPECT_EQ(&ad[i], &A[i])
            << "Raw array is not equal to A";
    }
}

/**
* @brief Conversion of Array1D to a const raw array 
*/
TEST_F(Array1D_test, conversion_operator_const)
{
    TNT::Array1D<double> A(a_size_, init_value_);

    // Conversion of A to ad
    // Information is cloned
    const double *ad = A;

    // Check if ad is equal to the init value
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(ad[i], init_value_)
            << "Raw array is not equal to A";
        EXPECT_EQ(&ad[i], &A[i])
            << "Raw array is not equal to A";
    }
}

/**
* @brief Test the copy method
*/
TEST_F(Array1D_test, copy_method)
{
    TNT::Array1D<double> A(a_size_, init_value_);
    TNT::Array1D<double> B(A.copy());

    // Check if ad is equal to the init value
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(A[i], B[i])
            << "B is not equal to A";
    }

    // Check the pointed address for data
    // Check if ad is equal to the init value
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_FALSE(&A[i] == &B[i])
            << "B address equal to A";
    }
}

/**
* @brief Test the inject method
*/
TEST_F(Array1D_test, inject_method)
{
    TNT::Array1D<double> A(a_size_, init_value_);
    TNT::Array1D<double> B(a_size_, 0.0);

    B.inject(A);

    // Check if B is equal to the init value given by A
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(A[i], B[i])
            << "B is not equal to A";
    }

    // Check the pointed address for data
    // Check if ad is equal to the init value
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_FALSE(&A[i] == &B[i])
            << "B address equal to A";
    }
}

/**
* @brief Test the ref method
*/
TEST_F(Array1D_test, ref_method)
{
    // Definition des arrays
    TNT::Array1D<double> A(a_size_, init_value_);
    TNT::Array1D<double> B(13);
    TNT::Array1D<double> C;
    TNT::Array1D<double> *D = new TNT::Array1D<double>(1);

    // B and C must be now a reference to A.
    B.ref(A);
    C.ref(A);
    D->ref(A);

    // Check dimension
    EXPECT_EQ(A.dim(), B.dim());
    EXPECT_EQ(A.dim(), C.dim());
    EXPECT_EQ(A.dim(), D->dim());

    // There must be 4 references to A
    EXPECT_EQ(A.ref_count(), 4);

    // Verify the values
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(A[i], B[i])
            << "B is not equal to A";
        EXPECT_EQ(A[i], C[i])
            << "C is not equal to A";
        EXPECT_EQ(A[i], D->operator[](i))
            << "D is not equal to A";
    }

    // Verify the adress
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(&A[i], &B[i])
            << "B is not equal to A";
        EXPECT_EQ(&A[i], &C[i])
            << "C is not equal to A";
        EXPECT_EQ(&A[i], &D->operator[](i))
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
* @brief Test the ref_count method
*/
TEST_F(Array1D_test, ref_count_method)
{
    // Definition des arrays
    TNT::Array1D<double> A(a_size_, init_value_);
    TNT::Array1D<double> Z(a_size_, init_value_);
    TNT::Array1D<double> B(13);
    TNT::Array1D<double> C;
    TNT::Array1D<double> *D = new TNT::Array1D<double>(1);

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
* @brief Test the = operator
*/
TEST_F(Array1D_test, equal_operator)
{
    // Definition des arrays
    TNT::Array1D<double> A(a_size_, init_value_);
    TNT::Array1D<double> Z(a_size_, init_value_);
    TNT::Array1D<double> B(13);
    TNT::Array1D<double> C;
    TNT::Array1D<double> *D = new TNT::Array1D<double>(1);

    // B and C must be now a reference to A.
    B = A;
    C = A;
    *D = A;

    // There must be 4 references to A
    EXPECT_EQ(A.ref_count(), 4)
        << "Ref_count is not correct";

    // Verify the values
    for (int i = 0; i < A.dim(); i++)
    {
        EXPECT_EQ(A[i], B[i])
            << "B is not equal to A";
        EXPECT_EQ(A[i], C[i])
            << "C is not equal to A";
        EXPECT_EQ(A[i], D->operator[](i))
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
* @brief Test the subarray method
*/
TEST_F(Array1D_test, subarray_method)
{
    TNT::Array1D<double> A(a_size_, init_value_);

    // Fill the array
    for (int i = 0; i < A.dim(); i++)
    {
        A[i] = i;
    }

    // Create an array from a subarray A
    TNT::Array1D<double> B = A.subarray(0, 9);

    EXPECT_EQ(B.dim(), 10);

    // Verify the address
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(&A[i], &B[i])
            << "B is not equal to A";
    }

    TNT::Array1D<double> C = A.subarray(9, 0);
    EXPECT_EQ(C.dim(), 0);

    TNT::Array1D<double> D = A.subarray(-1, 40);
    EXPECT_EQ(D.dim(), 0);

    TNT::Array1D<double> E = A.subarray(0, 43);
    EXPECT_EQ(E.dim(), 0);

    TNT::Array1D<double> F = A.subarray(1, 1);
    EXPECT_EQ(F.dim(), 1);
}

/**
* @}
*/