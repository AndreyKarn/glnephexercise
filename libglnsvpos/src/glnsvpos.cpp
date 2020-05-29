#include <func.h>
#include <glnsvpos.h>
#include <rungekutta.h>
#include <structures.h>

using namespace std;

uint64_t glnsvpos(bool RK_valid, double h) {

    // Class Ephemeris
    struct Ephemeris_s Eph;

    // Coordinates
    Eph.X = -8444572.27;
    Eph.Y = -8664957.52;
    Eph.Z = 22466454.10;
    // Velocity
    Eph.VX = 2983.60348;
    Eph.VY = -743.76965;
    Eph.VZ =  832.83615;
    // Acceleration
    Eph.AX = -0.0000028;
    Eph.AY =  0.0000019;
    Eph.AZ = -0.0000019;

    uint16_t Time_year = 2020;
    uint8_t Time_month = 2;
    uint8_t Time_day = 25;
    uint8_t Time_hour = 13;
    uint8_t Time_minutes = 45;
    uint8_t Time_seconds = 18;

    uint16_t year_idx = 0;

    // Time in Gln
    Eph.N4 = ((Time_year - 1996) / 4) + 1;
    while (Eph.N4 > 31) { // Учет 5-битности N4
        Eph.N4 -= 31;
        year_idx++;
    }

    Eph.NT = NT_calc( Eph.N4, Time_year, year_idx, Time_month, Time_day );

    Eph.tb = Time_seconds + Time_minutes*60 + Time_hour*60*60 + 10800;
    if (Eph.tb >= 24*60*60) {
        Eph.tb -= 24*60*60;
        Eph.NT++;
        if (Eph.NT >= 1462) Eph.N4++;
    }

    double GMST = GMST_calc( Eph.N4, Eph.NT );

    struct Ephemeris_s Eph0 = CrdTrnsf2Inertial( Eph, GMST );

    uint32_t tn = Eph.tb; // Текущее время
    uint32_t Toe = (12+3)*60*60; // Начальное время
    uint32_t Tof = (24+3)*60*60; // Конечное время

    uint64_t N2inc = (Tof - tn) / (double)h; // Количесвио отcчетов для времени большего текущего Eph.tb
    uint64_t N2dec = (tn - Toe) / (double)h; // Количесвио отcчетов для времени меньшего текущего Eph.tb
    uint64_t N = N2inc + N2dec; // Общее число отсчетов

    // Численное интегрирования для времени меньшего текущего Eph.tb
    struct Y_s *Ydec;
    Ydec = new struct Y_s[N2dec];
    Ydec[0].X = Eph0.X;
    Ydec[0].Y = Eph0.Y;
    Ydec[0].Z = Eph0.Z;
    Ydec[0].VX = Eph0.VX;
    Ydec[0].VY = Eph0.VY;
    Ydec[0].VZ = Eph0.VZ;

    if (RK_valid) RK( N2dec, -h, Ydec);

    // Численное интегрирования для времени большего текущего Eph.tb
    struct Y_s *Yinc;
    Yinc = new struct Y_s[N2inc];
    Yinc[0].X = Eph0.X;
    Yinc[0].Y = Eph0.Y;
    Yinc[0].Z = Eph0.Z;
    Yinc[0].VX = Eph0.VX;
    Yinc[0].VY = Eph0.VY;
    Yinc[0].VZ = Eph0.VZ;

    if (RK_valid) RK( N2inc, h, Yinc);

    // Формирование выходного массива
    struct Y_s *Yout;
    Yout = new struct Y_s[N];

    uint32_t i;

    // -----------------------------------------
    // M:  1   2   3   4   5    - adress
    //    x1  x2  x3  x4  x5    - value
    // C:  0   1   2   3   4    - adress
    //    x1  x2  x3  x4  x5    - value
    // if MATLAB array size N, then in C++ N-1
    // -----------------------------------------

    for (i = 0; i < N2dec; i++) {
        Yout[i] = Ydec[N2dec-i-1];
        //cout << " i = " << i << " N2dec-i = " << N2dec-i << endl;
    }
    for (i = 0; i < N2inc; i++) {
        Yout[i+N2dec] = Yinc[i];
        //cout << " i = " << i << " i+N2dec = " << i+N2dec << endl;
    }

    // Учет ускорений
    double tau = (double)Toe - (double)tn;
    for (i = 0; i < N; i++) {

        //cout << "tau = " << tau << endl;

        Yout[i].X += Eph0.AX * (tau * tau) / (double)2;
        Yout[i].Y += Eph0.AY * (tau * tau) / (double)2;
        Yout[i].Z += Eph0.AZ * (tau * tau) / (double)2;

        Yout[i].VX += Eph0.AX * tau;
        Yout[i].VY += Eph0.AY * tau;
        Yout[i].VZ += Eph0.AZ * tau;

        tau += (double)h;
    }

    write_struct_Y(Yout, N, "data_out.txt");

    // Очищение памяти
//    delete []Yout;
//    delete []Ydec;
//    delete []Yinc;

    return N;
}
