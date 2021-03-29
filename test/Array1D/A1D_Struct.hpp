#ifndef TEST_ARRAY1D_H
#define TEST_ARRAY1D_H

#include <gtest/gtest.h>
#include "tnt.h"
#include <fstream>
#include <iostream>

/**
* @brief General Array 1D test structure
*/
struct Array1D_test : public ::testing::Test
{
    int a_size_;
    double init_value_;
    double scalar_;
    virtual void SetUp() override
    {
        a_size_ = 42;
        init_value_ = 1.0;
        scalar_ = 10.0;
    }
    virtual void TearDown() override
    {
    }
};

#endif