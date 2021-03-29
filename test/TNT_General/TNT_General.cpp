
#include "tnt.h"
#include <gtest/gtest.h>
#include <fstream>
#include <iostream>


/**
* @brief General Array 1D test structure
*/
struct TNT_General_test : public ::testing::Test
{

    virtual void SetUp() override
    {
    }
    virtual void TearDown() override
    {
    }
};

/**
* @brief Testing the hypot declared function
* this is tested because of the abs method in 
the function. if no std is given abs only works 
with integers. 
*/
TEST_F(TNT_General_test, hypot)
{
    double leg_x, leg_y, result;
    leg_x = 3.;
    leg_y = 4.;

    result = TNT::hypot(leg_x, leg_y);

    EXPECT_EQ(result, 5.0);

    leg_x = 3.5;
    leg_y = 4.5;

    result = TNT::hypot(leg_x, leg_y);

    EXPECT_NEAR(result, 5.700877125, 0.000001);
};
