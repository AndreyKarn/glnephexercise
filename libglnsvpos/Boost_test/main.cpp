#define BOOST_TEST_MODULE Test
#include <boost/test/included/unit_test.hpp>

// **************** HOW TO START ****************
// 1. Build this project
// 2. Open Console
// 3. Move to dir "..\glnephexercise\libglnsvpos\build-Boost_test-Desktop_Qt_5_14_2_MinGW_64_bit-Debug\debug"
// 4. Launch "Boost_test.exe" with key " --log_level=test_suite --run_test=+test_net"
// Example "D:\Repository\glnephexercise\libglnsvpos\build-Boost_test-Desktop_Qt_5_14_2_MinGW_64_bit-Debug\debug>Boost_test.exe --log_level=test_suite --run_test=+test_net"
// **************** HOW TO START ****************

#include <func.h>
#include <glnsvpos.h>
#include <rungekutta.h>

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(test_NT_calc)
{
    // In data
    uint16_t Time_year = 2020;
    uint8_t Time_month = 2;
    uint8_t Time_day = 25;

    // expected
    uint8_t N4exp = 7;
    uint16_t NTexp = 57;

    // actual
    uint8_t N4act;
    uint16_t NTact;

    // Function
    uint16_t year_idx = 0;
    N4act = ((Time_year - 1996) / 4) + 1;
    while (N4act > 31) { // Учет 5-битности N4
        N4act -= 31;
        year_idx++;
    }
    NTact = NT_calc(N4act, Time_year, year_idx, Time_month, Time_day);

    // Result
    BOOST_CHECK( N4exp == N4act );
    BOOST_CHECK( NTexp == NTact );
}

//BOOST_AUTO_TEST_CASE(test_add)
//{
//    int a = 2;
//    int b = 3;
//    int c = add(a,b);

//    BOOST_CHECK( c == a+b );
//}

//BOOST_AUTO_TEST_CASE(test_mult)
//{
//    int a = 2;
//    int b = 3;
//    int c = mult(a,b);

//    BOOST_CHECK( c == a*b );
//}
