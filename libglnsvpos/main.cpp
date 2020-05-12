#include <QCoreApplication>
#include <iostream>

#include "include/libglnsvpos/glnsvpos.h"
#include "include/libglnsvpos/rungekutta.h"
#include "include/libglnsvpos/structures.h"

using namespace std;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

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
    uint16_t Time_month = 2;
    uint16_t Time_day = 25;
    uint16_t Time_hour = 13;
    uint16_t Time_minutes = 45;
    uint16_t Time_seconds = 18;

    // Time in Gln
    // TODO заменить на функцию
    Eph.N4 = ((Time_year - 1996) / 4) + 1;
    Eph.NT = 365*(Time_year-1996-4*(Eph.N4-1)) + 31 + Time_day + 1; // BUG Сделать нормально учет месяца
    Eph.tb = Time_seconds + Time_minutes*60 + Time_hour*60*60 + 10800;

    double GMST = GMST_calc( Eph.N4, Eph.NT);

    struct Ephemeris_s Eph0 = CrdTrnsf2Inertial(Eph, GMST);

    uint32_t tn = Eph.tb; // Текущее время
    double h = 1;  // Шаг
    uint32_t Toe = (12+3)*60*60; // Начальное время
    uint32_t Tof = (24+3)*60*60; // Конечное время

//    uint32_t Toe = tn-2; // Начальное время
//    uint32_t Tof = tn+2; // Конечное время

    uint32_t N2inc = (Tof - tn) / h; // Количесвио отcчетов для времяни большего текущего Eph.tb
    uint32_t N2dec = (tn - Toe) / h; // Количесвио отcчетов для времяни меньшего текущего Eph.tb
    uint32_t N = N2inc + N2dec + 1; // Общее число отсчетов
    // TODO еще один отcчет это текущее время, не забывать

    struct Y_s* Y;
    Y = new struct Y_s[N2inc];

    Y[0].F1 = Eph0.X;
    Y[0].F2 = Eph0.Y;
    Y[0].F3 = Eph0.Z;
    Y[0].F4 = Eph0.VX;
    Y[0].F5 = Eph0.VY;
    Y[0].F6 = Eph0.VZ;

    RK( N2dec, -h, Y);

    RK( N2inc, h, Y);

    cout << "add(2,2) = " << add(2,2) << "\n";
    cout << "mult(2,2) = " << mult(2,2) << "\n";

    return a.exec();
}
