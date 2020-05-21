#define BOOST_TEST_MODULE Test
#include <boost/test/included/unit_test.hpp>

// **************** HOW TO START ****************
// 1. Build this project
// 2. Open Console
// 3. Move to path "..\glnephexercise\libglnsvpos\build-Boost_test-Desktop_Qt_5_14_2_MinGW_64_bit-Debug\debug"
// 4. Launch "Boost_test.exe" with key " --log_level=test_suite --run_test=+test_net"
// Example "D:\Repository\glnephexercise\libglnsvpos\build-Boost_test-Desktop_Qt_5_14_2_MinGW_64_bit-Debug\debug>Boost_test.exe --log_level=test_suite --run_test=+test_net"
// **************** HOW TO START ****************

#include <func.h>
#include <glnsvpos.h>
#include <rungekutta.h>

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(test_add)
{
    int a = 2;
    int b = 3;
    int c = add(a,b);

    BOOST_CHECK( c == a+b );
}

BOOST_AUTO_TEST_CASE(test_mult)
{
    int a = 2;
    int b = 3;
    int c = mult(a,b);

    BOOST_CHECK( c == a*b );
}
