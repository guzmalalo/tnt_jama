#ifndef TEST_ARRAY2D_H
#define TEST_ARRAY2D_H

#include "tnt.h"
#include <gtest/gtest.h>
#include <fstream>
#include <iostream>

/**
* @brief General Array 1D test structure
*/
struct Array2D_test : public ::testing::Test
{
    int m_;
    int n_;
    double init_value_;
    double scalar_;
    virtual void SetUp() override
    {
        m_ = 3;
        n_ = 3;
        init_value_ = 1.0;
        scalar_ = 10.0;
    }
    virtual void TearDown() override
    {
    }
};

#endif