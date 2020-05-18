#include "tst_test2.h"

#include <gtest/gtest.h>

#include <func.h>
#include <glnsvpos.h>
#include <rungekutta.h>

using namespace testing;

TEST(test_case_name1, test_name1)
{
    EXPECT_EQ (10, 10);
}

TEST(test_case_name1, test_name2)
{
    ASSERT_NE (10, 10);
}

TEST(test_case_name2, test_name3)
{
    EXPECT_LE (10, 10);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
