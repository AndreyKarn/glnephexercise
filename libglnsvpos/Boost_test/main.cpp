#define BOOST_TEST_MODULE Test
#include <boost/test/included/unit_test.hpp>

// **************** HOW TO START ****************
// 1. Build this project
// 2. Open Console
// 3. Move to dir "..\glnephexercise\libglnsvpos\build-Boost_test-Desktop_Qt_5_14_2_MinGW_64_bit-Debug\debug"
// 4. Launch "Boost_test.exe" with key " --log_level=test_suite --run_test=+test_net"
// Example "D:\Repository\glnephexercise\libglnsvpos\build-Boost_test-Desktop_Qt_5_14_2_MinGW_64_bit-Debug\debug>Boost_test.exe --log_level=test_suite --run_test=+test_net"
// **************** HOW TO START ****************

#include <math.h>
#include <iostream>

#include <func.h>
#include <glnsvpos.h>
#include <rungekutta.h>
#include <structures.h>

namespace utf = boost::unit_test;
using namespace std;

//BOOST_AUTO_TEST_CASE(test_NT_calc)
//{
//    // In data
//    uint16_t Time_year = 2020;
//    uint8_t Time_month = 2;
//    uint8_t Time_day = 25;

//    // expected
//    uint8_t N4exp = 7;
//    uint16_t NTexp = 57;

//    // actual
//    uint8_t N4act;
//    uint16_t NTact;

//    // Function
//    uint16_t year_idx = 0;
//    N4act = ((Time_year - 1996) / 4) + 1;
//    while (N4act > 31) { // Учет 5-битности N4
//        N4act -= 31;
//        year_idx++;
//    }
//    NTact = NT_calc(N4act, Time_year, year_idx, Time_month, Time_day);

//    // Result
//    BOOST_CHECK( N4exp == N4act );
//    BOOST_CHECK( NTexp == NTact );
//}

//BOOST_AUTO_TEST_CASE(test_mem_leak)
//{
//    struct Y_s *Ydec;
//    Ydec = new struct Y_s[100];

//    BOOST_CHECK( true );
//}


BOOST_AUTO_TEST_CASE(test_glnsvpos)
{
    uint32_t i;
    uint64_t N;
    double h = 1;

    N = glnsvpos(1, h);

    struct Y_s *Y_model;
    Y_model = new struct Y_s[N];

    struct Y_s *Y_data;
    Y_data = new struct Y_s[N];

    struct Y_s *Y_delta;
    Y_delta = new struct Y_s[N];

    // Чтение из файла (матлаб)
    read_struct_Y(Y_model, N, "Matlab_data_for_h1.txt");

    // Чтение из файла (С++)
    ifstream file("data_out.txt");
    if (file.is_open()) { //Если открытие файла прошло успешно
        string line; //Строчка текста
        uint32_t i = 0;
        while (getline(file, line)) {
            istringstream iss(line);
            iss >> Y_data[i].X >> Y_data[i].Y >> Y_data[i].Z >> Y_data[i].VX >> Y_data[i].VY >> Y_data[i].VZ;
            i++;
        }
    }
    else printf("Error. File: %s, Line: %d\n", __FILE__, __LINE__);

    double deltaXmax = 0;
    double deltaYmax = 0;
    double deltaZmax = 0;

    for (i = 0; i <= N; i++) {
        Y_delta[i].X = Y_data[i].X - Y_model[i].X;
        Y_delta[i].Y = Y_data[i].Y - Y_model[i].Y;
        Y_delta[i].Z = Y_data[i].Z - Y_model[i].Z;

        if ( fabs(Y_delta[i].X) > fabs(deltaXmax) ) deltaXmax = Y_delta[i].X;
        if ( fabs(Y_delta[i].Y) > fabs(deltaYmax) ) deltaYmax = Y_delta[i].Y;
        if ( fabs(Y_delta[i].Z) > fabs(deltaZmax) ) deltaZmax = Y_delta[i].Z;
    }

    cout << "Max delta X = " << deltaXmax << endl;
    cout << "Max delta Y = " << deltaYmax << endl;
    cout << "Max delta Z = " << deltaZmax << endl;

    delete []Y_data;
    delete []Y_model;
    delete []Y_delta;

    double Epsilon = 0.2e-6;
    if ( deltaXmax <= Epsilon && deltaYmax <= Epsilon && deltaZmax <= Epsilon)
        BOOST_CHECK( true );
    else
        BOOST_CHECK( false );
}

BOOST_AUTO_TEST_CASE(test_RK)
{
    uint64_t N = 100000;
    double h = 1;

    struct Y_s *Y_data;
    Y_data = new struct Y_s[N];

    Y_data[0].X = (double)rand();
    Y_data[0].Y = (double)rand();
    Y_data[0].Z = (double)rand();
    Y_data[0].VX = (double)rand();
    Y_data[0].VY = (double)rand();
    Y_data[0].VZ = (double)rand();

    if ( !RK( N, h, Y_data) )
        BOOST_CHECK( true );
    else
        BOOST_CHECK( false );
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
