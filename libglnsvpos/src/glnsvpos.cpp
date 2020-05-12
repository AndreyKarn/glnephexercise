#include <include/libglnsvpos/glnsvpos.h>
#include <include/libglnsvpos/rungekutta.h>
#include <include/libglnsvpos/structures.h>

using namespace std;

double GMST_calc(uint8_t N4, uint16_t NT) {
    // Текущая Юлианская дата на 0 часов шкалы МДВ
    double JD0 = 1461 * (N4 - 1) + NT + 2450082.5 - (NT - 3) / 25;
    // Время от эпохи 2000 г 1 января 12 ч (UTC(SU))
    double T_delta = (JD0 - 2451545) / 36525;
    // Угол поворота Земли [рад]
    double ERA = 2 * M_PI * ( 0.7790572732640 + 1.00273781191135448 * (JD0 - 2451545));
    // Среднее звездное время по Гринвичу [рад]
    double GMST = ERA;
    GMST += 0.0000000703270726;
    GMST += 0.0223603658710194 * T_delta;
    GMST += 0.0000067465784654 * T_delta * T_delta;
    GMST -= 0.0000000000021332 * T_delta * T_delta * T_delta;
    GMST -= 0.0000000001452308 * T_delta * T_delta * T_delta * T_delta;
    GMST -= 0.0000000000001784 * T_delta * T_delta * T_delta * T_delta * T_delta;

    return GMST;
}

Ephemeris_s CrdTrnsf2Inertial(struct Ephemeris_s Eph, double GMST) {
    struct Ephemeris_s Eph0;

    double Omega_E = 7.2921151467e-5;

    //double Theta_Ge = GMST + Omega_E * (Eph.tb - 3 * 60 * 60);
    double Theta_Ge = 4.636890363514948e+04;

    // Координаты:
    Eph0.X = Eph.X * cos(Theta_Ge) - Eph.Y * sin(Theta_Ge);
    Eph0.Y = Eph.X * sin(Theta_Ge) + Eph.Y * cos(Theta_Ge);
    Eph0.Z = Eph.Z;

    // Скорости:
    Eph0.VX = Eph.VX * cos(Theta_Ge) - Eph.VY * sin(Theta_Ge) - Omega_E * Eph0.Y;
    Eph0.VY = Eph.VX * sin(Theta_Ge) + Eph.VY * cos(Theta_Ge) + Omega_E * Eph0.X;
    Eph0.VZ = Eph.VZ;

    // Ускорения:
    Eph0.AX = Eph.AX * cos(Theta_Ge) - Eph.AY * sin(Theta_Ge);
    Eph0.AY = Eph.AX * sin(Theta_Ge) + Eph.AY * cos(Theta_Ge);
    Eph0.AZ = Eph.AZ;

    Eph0.N4 = Eph.N4;
    Eph0.NT = Eph.NT;
    Eph0.tb = Eph.tb;

    return Eph0;
}

int add(int a, int b) {
    return mult(a,b) + b;
}




