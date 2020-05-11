#include <QCoreApplication>
#include <iostream>

#include "include/libglnsvpos/glnsvpos.h"
#include "include/libglnsvpos/rungekutta.h"
#include "include/libglnsvpos/students.h"

using namespace std;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // Class Ephemeris
    Ephemeris Eph;

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



    cout << "add(2,2) = " << add(2,2) << "\n";
    cout << "mult(2,2) = " << mult(2,2) << "\n";

    return a.exec();
}
